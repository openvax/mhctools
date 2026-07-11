# Copyright (c) 2014-2017. Mount Sinai School of Medicine
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#       http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

from collections import defaultdict
import logging
from multiprocessing import cpu_count
from subprocess import check_output
import tempfile

from typechecks import require_string, require_integer, require_iterable_of
from .allele_normalization import (
    normalize_allele_name,
    normalize_allele_name_or_raw,
    AlleleParseError,
)

from .base_predictor import BasePredictor
from .unsupported_allele import UnsupportedAllele
from .process_helpers import run_command
from .cleanup_context import CleanupFiles
from .input_file_formats import create_input_peptides_files
from .process_helpers import run_multiple_commands_redirect_stdout
from .binding_prediction_collection import BindingPredictionCollection
from .pred import PeptideResult

logger = logging.getLogger(__name__)

# Upper bound on how many alleles the "auto" policy will pack into one
# predictor invocation. Batching amortizes the per-process startup cost
# (~100ms for netMHCpan) across alleles, but the marginal cost of an extra
# allele in a running process is small (~15ms), so the benefit saturates
# quickly. Past this point, larger batches mostly add downside: a single
# failing allele takes out the whole batch, one process holds more output in
# memory, and there are fewer processes to spread across cores. 20 keeps
# ~85% of the amortization while bounding those risks.
AUTO_MAX_ALLELES_PER_COMMAND = 20

# Caches keyed by (command, supported_allele_flag), memoized per process so a
# binary's supported-allele list is only read/parsed once no matter how many
# predictors are constructed.
#
# _supported_alleles_cache holds the RAW names exactly as the predictor prints
# them (e.g. `netMHCpan -listMHC`). Building it is cheap: split lines, no
# mhcgnomes. User alleles are validated against it with a fast path that only
# runs mhcgnomes on the user's handful of alleles, so the ~13k-name list is
# never normalized for the common (human HLA) case.
_supported_alleles_cache = {}

# _normalized_supported_cache holds the lazily-built {normalized_name: raw_name}
# mapping, computed only when the fast path misses (a supported allele whose
# name our normalizer rewrites into a spelling the predictor doesn't accept,
# e.g. many non-human BoLA/Mamu names). The raw_name is the predictor's own
# spelling, which is what must be passed back on the command line.
_normalized_supported_cache = {}


def _unique_in_order(values):
    return list(dict.fromkeys(values))


def _prefer_allele_spelling(new, existing):
    """
    When several raw -listMHC spellings normalize to the same allele, pick the
    one most likely accepted on `-a`, deterministically (independent of set
    iteration order):

      1. a colon-containing spelling (e.g. BoLA-1:00901, which netMHCpan
         accepts, over BoLA-100901, which it rejects);
      2. among those, the spelling whose colon sits earliest — i.e. right
         after the gene (Mamu-A1:00101 over Mamu-A1001:01);
      3. the lexicographically smaller string as a final stable tie-break.
    """
    def rank(name):
        colon = name.find(":")
        # higher is preferred: has-colon first, then earliest colon
        return (colon >= 0, -colon if colon >= 0 else 1)

    rank_new, rank_existing = rank(new), rank(existing)
    if rank_new != rank_existing:
        return new if rank_new > rank_existing else existing
    return min(new, existing)


class BaseCommandlinePredictor(BasePredictor):
    """
    Base class for MHC binding predictors that run a local external
    program and write their output to a local file.
    """
    def __init__(
            self,
            program_name,
            alleles,
            parse_output_fn,
            supported_alleles_flag,
            input_file_flag,
            length_flag,
            allele_flag,
            peptide_mode_flags=["-p"],
            tempdir_flag=None,
            extra_flags=[],
            max_peptides_per_file=10 ** 4,
            max_alleles_per_command=1,
            process_limit=-1,
            default_peptide_lengths=[9],
            group_peptides_by_length=False,
            min_peptide_length=8,
            max_peptide_length=None,
            parse_to_preds_fn=None,):
        """
        Parameters
        ----------
        program_name : str
            Name of prediction program to run
            (e.g. "netMHCcons" or "netMHCIIpan")

        alleles : list of str
            MHC alleles

        supported_alleles_flag : str
            Flag to pass to the predictor to get a list of supported alleles
            (e.g. "-A", "-list", "-listMHC")

        parse_output_fn : fn
            Takes the stdout string from the predictor and returns a collection
            of BindingPrediction objects

        input_file_flag : str
            How to specify the input FASTA file of source sequences (e.g. "-f")

        length_flag : str
            How to specify the desired predicted peptide length (e.g. "-length")

        allele_flag : str
            How to specify the allele we want predictions for (e.g. "-a")

        peptide_mode_flags : list of str
            How to switch from the default FASTA subsequences input mode to
            where peptides are explicitly given one per line of a text file.

        tempdir_flag : str, optional
            How to specify the predictor's temporary directory (e.g. "-tdir")

        extra_flags : list of str
            Extra flags to pass to the predictor

        max_peptides_per_file : int, optional
            Maximum number of lines per file when predicting peptides directly.

        max_alleles_per_command : int, "auto", or None, optional
            How many alleles to pass to a single invocation of the predictor
            via a comma-separated allele flag (e.g. ``-a A0201,B3502``). Only
            predictors whose allele flag accepts a comma-separated list (e.g.
            the netMHCpan family) should batch more than one allele.

              - ``1`` (default): one allele per command, i.e. one process per
                (input file, allele). Preserves the historical behavior.
              - ``"auto"``: batch alleles to amortize the per-process startup
                cost, while keeping enough parallel processes (input files x
                allele groups) to use the available cores and capping the
                group size at AUTO_MAX_ALLELES_PER_COMMAND. Speeds up
                many-allele runs without pessimizing small runs on
                otherwise-idle cores or over-batching into fragile mega-calls.
              - ``None`` or ``<= 0``: all alleles in a single command
                (unbounded batching).
              - ``k > 1``: at most ``k`` alleles per command.

        process_limit : int, optional
            Maximum number of parallel processes to start
            (0 for no limit, -1 for use all available processors)

        default_peptide_lengths : list of int, optional
            When making predictions across subsequences of protein sequences,
            what peptide lengths to predict for.

        group_peptides_by_length : bool
            Run commandline predictor on groups of peptides of equal length

        min_peptide_length : int
            Shortest peptide this predictor can handle

        max_peptide_length : int
            Longest peptide this predictor can handle
        """
        require_string(program_name, "Predictor program name")
        self.program_name = program_name

        if supported_alleles_flag is not None:
            require_string(supported_alleles_flag, "Supported alleles flag")
        self.supported_alleles_flag = supported_alleles_flag

        require_string(input_file_flag, "Input file flag")
        self.input_file_flag = input_file_flag

        require_string(length_flag, "Peptide length flag")
        self.length_flag = length_flag

        require_string(allele_flag, "Allele flag")
        self.allele_flag = allele_flag

        require_iterable_of(peptide_mode_flags, str)
        self.peptide_mode_flags = peptide_mode_flags

        if tempdir_flag is not None:
            require_string(tempdir_flag, "Temporary directory flag")
        self.tempdir_flag = tempdir_flag

        require_iterable_of(extra_flags, str)
        self.extra_flags = extra_flags

        require_integer(
            max_peptides_per_file,
            "Maximum number of lines in a peptides input file")
        self.max_peptides_per_file = max_peptides_per_file

        if max_alleles_per_command not in (None, "auto"):
            require_integer(
                max_alleles_per_command,
                "Maximum number of alleles per command")
        self.max_alleles_per_command = max_alleles_per_command

        require_integer(process_limit, "Maximum number of processes")
        self.process_limit = process_limit

        self.parse_output_fn = parse_output_fn
        self.parse_to_preds_fn = parse_to_preds_fn

        if isinstance(default_peptide_lengths, int):
            default_peptide_lengths = [default_peptide_lengths]

        self.group_peptides_by_length = group_peptides_by_length

        if self.supported_alleles_flag:
            # Raw names exactly as the predictor prints them; validation and
            # normalization happen lazily in _resolve_supported_alleles().
            self._supported_allele_names = self._determine_supported_alleles(
                self.program_name,
                self.supported_alleles_flag)
        else:
            # if we're not running the tool to determine supported alleles
            # then at least try running it by itself to determine if it's
            # it's present
            try:
                run_command([self.program_name])
            except Exception:
                raise SystemError("Failed to run %s" % self.program_name)
            self._supported_allele_names = None

        # Normalize (and dedupe homozygous) the requested alleles without
        # validating against the supported list — we validate below against the
        # raw list so we can hand the predictor back its own spelling.
        BasePredictor.__init__(
            self,
            alleles=alleles,
            valid_alleles=None,
            default_peptide_lengths=default_peptide_lengths,
            min_peptide_length=min_peptide_length,
            max_peptide_length=max_peptide_length,
            # Non-human alleles the normalizer can't parse (e.g. H-2-Qa1) are
            # kept verbatim and validated against the tool's raw -listMHC list.
            keep_unparseable_alleles=True)

        self._resolve_supported_alleles()

    @staticmethod
    def _determine_supported_alleles(command, supported_allele_flag):
        """
        Ask the commandline predictor (e.g. `netMHCpan -listMHC`) which alleles
        it supports and return the RAW names it prints, verbatim.

        No allele-name normalization happens here: parsing all ~13k names
        through mhcgnomes is the slow part, and it is deferred to
        _normalized_supported_alleles() which only runs when the fast
        validation path misses. Memoized per (command, supported_allele_flag).
        """
        cache_key = (command, supported_allele_flag)
        cached = _supported_alleles_cache.get(cache_key)
        if cached is not None:
            return cached
        try:
            # convert to str since Python3 returns a `bytes` object
            supported_alleles_output = check_output([
                command, supported_allele_flag
            ])
            supported_alleles_str = supported_alleles_output.decode("ascii", "ignore")
            assert len(supported_alleles_str) > 0, \
                '%s returned empty allele list' % command
            supported_alleles = {
                line.strip()
                for line in supported_alleles_str.split("\n")
                if line.strip() and not line.strip().startswith("#")
            }
            if len(supported_alleles) == 0:
                raise ValueError("Unable to determine supported alleles")

            # Sanity check that this really is an allele list and not, e.g., a
            # mismatched executable's usage/error text (which would otherwise
            # sail through as a set of "raw" names and only fail later as an
            # UnsupportedAllele). At least one name must parse as an allele.
            # Probe lazily and stop at the first success — O(1) for a real
            # list, so this does not reintroduce the up-front full parse.
            def parses(name):
                try:
                    normalize_allele_name(name)
                    return True
                except AlleleParseError:
                    return False
            if not any(parses(name) for name in supported_alleles):
                raise ValueError(
                    "No parseable alleles in output of %s %s "
                    "(possibly an incorrect executable version?)" % (
                        command, supported_allele_flag))

            _supported_alleles_cache[cache_key] = supported_alleles
            return supported_alleles
        except Exception as e:
            logger.exception(e)
            raise SystemError("Failed to run %s %s. Possibly an incorrect executable version?" % (
                command,
                supported_allele_flag))

    @staticmethod
    def _normalized_supported_alleles(command, supported_allele_flag, raw_names):
        """
        Lazily build (and memoize) a {normalized_name: raw_name} mapping over
        the predictor's supported alleles. This is where the expensive
        mhcgnomes normalization of the whole list happens; it is only called
        when a requested allele is not found verbatim, so the common (human
        HLA) case never pays for it.

        raw_name is the predictor's own spelling, which is what round-trips on
        the command line. When several raw spellings normalize to the same name
        (netMHCpan lists e.g. both `BoLA-1:00901` and `BoLA-100901`, or two
        colon forms for one Mamu allele), _prefer_allele_spelling picks one
        deterministically.
        """
        cache_key = (command, supported_allele_flag)
        cached = _normalized_supported_cache.get(cache_key)
        if cached is not None:
            return cached
        mapping = {}
        skipped = 0
        for name in raw_names:
            try:
                key = normalize_allele_name(name)
            except AlleleParseError:
                # Some non-human/edge-case names (e.g. BoLA-amani.1, H-2-Qa1)
                # can't be normalized; they remain reachable only by verbatim
                # match on the fast path.
                skipped += 1
                continue
            existing = mapping.get(key)
            if existing is None:
                mapping[key] = name
            else:
                mapping[key] = _prefer_allele_spelling(name, existing)
        if skipped:
            logger.debug(
                "%s %s: %d of %d supported allele names could not be normalized",
                command, supported_allele_flag, skipped, len(raw_names))
        _normalized_supported_cache[cache_key] = mapping
        return mapping

    def _resolve_supported_allele_cli_names(self, alleles):
        """
        Validate alleles against the predictor's supported-allele list and
        return the exact spelling to pass on the command line for each.

        Fast path: if prepare_allele_name(allele) is printed verbatim by the
        predictor, use it directly — no full-list normalization needed. Slow
        path (only on a miss): normalize the whole supported list once and look
        the allele up by its normalized name, using the predictor's own
        spelling on the command line so names our normalizer rewrites (many
        non-human alleles) still round-trip.
        """
        raw_names = self._supported_allele_names
        allele_cli_names = {}
        if raw_names is None:
            return allele_cli_names
        unsupported = []
        normalized_map = None
        for allele in _unique_in_order(alleles):
            # Fast path: the predictor's expected spelling is printed verbatim.
            # prepare_allele_name can itself reject an allele (e.g. netMHCIIpan
            # re-parses and raises on unexpected genes); if so, fall through to
            # the slow path, which resolves via the predictor's own -listMHC
            # spelling and never needs prepare_allele_name.
            try:
                cli_name = self.prepare_allele_name(allele)
            except Exception:
                cli_name = None
            if cli_name is not None and cli_name in raw_names:
                allele_cli_names[allele] = cli_name
                continue
            if normalized_map is None:
                normalized_map = self._normalized_supported_alleles(
                    self.program_name,
                    self.supported_alleles_flag,
                    raw_names)
            original = normalized_map.get(allele)
            if original is not None:
                allele_cli_names[allele] = original
            else:
                unsupported.append(allele)
        if unsupported:
            raise UnsupportedAllele(
                "Unsupported alleles for %s: %s\n"
                "Run command %s %s to see a list of valid alleles" % (
                    self.program_name,
                    unsupported,
                    self.program_name,
                    self.supported_alleles_flag))
        return allele_cli_names

    def _resolve_supported_alleles(self):
        """
        Validate self.alleles and record command-line spellings for them.
        """
        self._allele_cli_names = self._resolve_supported_allele_cli_names(
            self.alleles)

    def prepare_allele_name(self, allele_name):
        """
        How does the predictor expect to see allele names?
        """
        return allele_name.replace("*", "")

    def _cli_allele_name(self, allele_name, allele_cli_names=None):
        """
        The spelling to pass on the command line for a (normalized) allele.

        Prefers the predictor's own -listMHC spelling recorded during
        validation (so names like BoLA-1:00901 round-trip instead of being
        rewritten to a rejected form); falls back to prepare_allele_name for
        predictors that don't enumerate their supported alleles.
        """
        if allele_cli_names:
            resolved = allele_cli_names.get(allele_name)
            if resolved is not None:
                return resolved
        resolved = getattr(self, "_allele_cli_names", {}).get(allele_name)
        if resolved is not None:
            return resolved
        return self.prepare_allele_name(allele_name)

    def _auto_allele_group_size(self, n_alleles, n_input_files):
        """
        Group size for max_alleles_per_command="auto": batch alleles to
        amortize the per-process startup cost, but (a) keep enough parallel
        processes (n_input_files * groups_per_file) to occupy the cores the
        command runner will use, and (b) never exceed
        AUTO_MAX_ALLELES_PER_COMMAND, past which batching stops paying off.
        """
        if self.process_limit and self.process_limit > 0:
            target = self.process_limit
        else:
            target = cpu_count()
        n_input_files = max(1, n_input_files)
        # allele groups per file needed to reach `target` processes (ceil div)
        groups_per_file = max(1, -(-target // n_input_files))
        groups_per_file = min(groups_per_file, n_alleles)
        # group size = ceil(n_alleles / groups_per_file), capped
        group_size = max(1, -(-n_alleles // groups_per_file))
        return min(group_size, AUTO_MAX_ALLELES_PER_COMMAND)

    def _allele_groups(self, n_input_files=1, alleles=None):
        """
        Partition alleles into groups, one command (one predictor
        process) per group per input file. Group size is controlled by
        max_alleles_per_command (see __init__).

        Returns a list of lists of allele names. Empty if there are no
        alleles.
        """
        if alleles is None:
            alleles = self.alleles
        alleles = list(alleles)
        if not alleles:
            return []
        n = self.max_alleles_per_command
        if n == "auto":
            group_size = self._auto_allele_group_size(len(alleles), n_input_files)
        elif n is None or n <= 0 or n >= len(alleles):
            group_size = len(alleles)
        else:
            group_size = n
        return [alleles[i:i + group_size]
                for i in range(0, len(alleles), group_size)]

    def _build_command(
            self,
            input_filename,
            alleles,
            length=None,
            temp_dirname=None,
            peptide_mode=False,
            allele_cli_names=None):
        # accept either a single allele name or a list of them
        if isinstance(alleles, str):
            alleles = [alleles]
        args = [self.program_name]
        if peptide_mode:
            args.extend(self.peptide_mode_flags)
        allele_arg = ",".join(
            self._cli_allele_name(allele, allele_cli_names)
            for allele in alleles)
        args.extend([self.allele_flag, allele_arg])
        if length:
            args.extend([self.length_flag, str(length)])
        if self.tempdir_flag and temp_dirname:
            args.extend([self.tempdir_flag, temp_dirname])
        args.extend(self.extra_flags)
        if self.input_file_flag:
            args.extend([self.input_file_flag, input_filename])
        else:
            args.append(input_filename)
        return args

    def _run_commands_and_collect_predictions(
            self,
            commands,
            input_filenames,
            temp_dir_list,
            sequence_key_mapping=None):
        if sequence_key_mapping is None:
            sequence_key_mapping = defaultdict(lambda: "seq")
        binding_predictions = []

        # Cleanup either when finished or if an exception gets raised by
        # deleting the input and output files
        filenames_to_delete = list(input_filenames) + [
            f.name for f in commands.keys()]
        with CleanupFiles(
                filenames=filenames_to_delete,
                directories=temp_dir_list):
            run_multiple_commands_redirect_stdout(
                commands,
                print_commands=True,
                process_limit=self.process_limit)
            for output_file, command in commands.items():
                # subprocess wrote via its own fd; no flush/fsync needed
                output_file.seek(0)
                file_contents = output_file.read()
                output_file.close()
                binding_predictions.extend(
                    self.parse_output_fn(
                        stdout=file_contents,
                        sequence_key_mapping=sequence_key_mapping,
                        prediction_method_name=self.program_name))

        if len(binding_predictions) == 0:
            logger.warning("No binding predictions from %s" % self.program_name)
        return BindingPredictionCollection(binding_predictions)

    def _run_commands_and_collect_preds(
            self,
            commands,
            input_filenames,
            temp_dir_list,
            sequence_key_mapping=None):
        """Like _run_commands_and_collect_predictions but returns list[Pred]."""
        if sequence_key_mapping is None:
            sequence_key_mapping = defaultdict(lambda: "seq")
        all_preds = []

        filenames_to_delete = list(input_filenames) + [
            f.name for f in commands.keys()]
        with CleanupFiles(
                filenames=filenames_to_delete,
                directories=temp_dir_list):
            run_multiple_commands_redirect_stdout(
                commands,
                print_commands=True,
                process_limit=self.process_limit)
            for output_file, command in commands.items():
                # subprocess wrote via its own fd; no flush/fsync needed
                output_file.seek(0)
                file_contents = output_file.read()
                output_file.close()
                all_preds.extend(
                    self.parse_to_preds_fn(
                        stdout=file_contents,
                        sequence_key_mapping=sequence_key_mapping,
                        predictor_name=self.program_name))

        if len(all_preds) == 0:
            logger.warning("No predictions from %s" % self.program_name)

        # Group by (peptide, offset, source) into PeptideResult
        groups = defaultdict(list)
        for pred in all_preds:
            key = (pred.peptide, pred.offset, pred.source_sequence_name)
            groups[key].append(pred)
        return [PeptideResult(preds=tuple(preds)) for preds in groups.values()]

    def _build_peptide_commands(
            self,
            peptides,
            alleles,
            allele_cli_names=None):
        self._check_peptide_inputs(peptides)
        input_filenames = create_input_peptides_files(
            peptides,
            max_peptides_per_file=self.max_peptides_per_file,
            group_by_length=self.group_peptides_by_length)
        logger.debug("Created %d input files" % len(input_filenames))
        commands = {}
        dirs = []

        for i, input_filename in enumerate(input_filenames):
            for j, allele_group in enumerate(
                    self._allele_groups(
                        n_input_files=len(input_filenames),
                        alleles=alleles)):
                if self.tempdir_flag:
                    temp_dirname = tempfile.mkdtemp(
                        prefix="tmp_%d_%d_%s" % (
                            i,
                            j,
                            self.program_name),
                        suffix="XXXXXX")
                    logger.debug(
                        "Created temporary directory %s for alleles %s",
                        temp_dirname,
                        allele_group)
                    dirs.append(temp_dirname)
                else:
                    temp_dirname = None
                output_file = tempfile.NamedTemporaryFile(
                    "w+",
                    prefix="%s_output_length_%d_%d" % (
                        self.program_name, i, j),
                    delete=False)
                commands[output_file] = self._build_command(
                    input_filename=input_filename,
                    alleles=allele_group,
                    peptide_mode=True,
                    temp_dirname=temp_dirname,
                    allele_cli_names=allele_cli_names)
        return commands, input_filenames, dirs

    def _predict_binding_predictions_for_alleles(
            self,
            peptides,
            alleles,
            allele_cli_names=None):
        commands, input_filenames, dirs = self._build_peptide_commands(
            peptides=peptides,
            alleles=alleles,
            allele_cli_names=allele_cli_names)
        results = self._run_commands_and_collect_predictions(
            commands=commands,
            input_filenames=input_filenames,
            temp_dir_list=dirs)
        self._check_results(
            results,
            peptides=peptides,
            alleles=alleles)
        return results

    def _predict_for_alleles(self, peptides, alleles, allele_cli_names=None):
        if self.parse_to_preds_fn is None:
            collection = self._predict_binding_predictions_for_alleles(
                peptides=peptides,
                alleles=alleles,
                allele_cli_names=allele_cli_names)
            return collection.to_peptide_preds(kind=self._default_pred_kind())

        commands, input_filenames, dirs = self._build_peptide_commands(
            peptides=peptides,
            alleles=alleles,
            allele_cli_names=allele_cli_names)
        return self._run_commands_and_collect_preds(
            commands=commands,
            input_filenames=input_filenames,
            temp_dir_list=dirs)

    def _check_pair_inputs(self, peptides, alleles=None):
        if alleles is None:
            pairs = list(peptides)
            peptide_list = []
            allele_list = []
            for pair in pairs:
                try:
                    peptide, allele = pair
                except (TypeError, ValueError):
                    raise ValueError(
                        "Expected (peptide, allele) pairs, got %r" % (pair,))
                peptide_list.append(peptide)
                allele_list.append(allele)
        else:
            peptide_list = list(peptides)
            if isinstance(alleles, str):
                allele_list = [alleles]
            else:
                allele_list = list(alleles)
            if len(peptide_list) != len(allele_list):
                raise ValueError(
                    "peptides length %d != alleles length %d"
                    % (len(peptide_list), len(allele_list)))

        self._check_peptide_inputs(peptide_list)
        require_iterable_of(allele_list, str, "HLA alleles")
        allele_list = [
            normalize_allele_name_or_raw(allele)
            for allele in allele_list]
        return peptide_list, allele_list

    @staticmethod
    def _result_lookup_for_allele(results, allele):
        lookup = {}
        for peptide_result in results:
            matching_preds = tuple(
                pred for pred in peptide_result.preds
                if pred.allele == allele)
            if not matching_preds:
                continue
            peptide = matching_preds[0].peptide
            existing = lookup.get(peptide)
            if existing is None:
                lookup[peptide] = PeptideResult(preds=matching_preds)
            else:
                lookup[peptide] = PeptideResult(
                    preds=existing.preds + matching_preds)
        return lookup

    def predict_pairs(self, peptides, alleles=None):
        """Predict explicit peptide/allele pairs.

        Parameters
        ----------
        peptides : list of str or iterable of (str, str)
            Peptides to score. If *alleles* is omitted, this must be an
            iterable of ``(peptide, allele)`` pairs.
        alleles : list of str, optional
            Alleles parallel to *peptides*.

        Returns
        -------
        list of PeptideResult
            One entry per input pair, in order.
        """
        peptide_list, allele_list = self._check_pair_inputs(peptides, alleles)
        if not peptide_list:
            return []

        unique_alleles = _unique_in_order(allele_list)
        allele_cli_names = self._resolve_supported_allele_cli_names(
            unique_alleles)
        indices_by_allele = defaultdict(list)
        for index, allele in enumerate(allele_list):
            indices_by_allele[allele].append(index)

        results = [None] * len(peptide_list)
        for allele in unique_alleles:
            indices = indices_by_allele[allele]
            group_peptides = _unique_in_order(
                peptide_list[index] for index in indices)
            group_results = self._predict_for_alleles(
                peptides=group_peptides,
                alleles=[allele],
                allele_cli_names=allele_cli_names)
            lookup = self._result_lookup_for_allele(group_results, allele)
            for index in indices:
                peptide = peptide_list[index]
                peptide_result = lookup.get(peptide)
                if peptide_result is None:
                    raise ValueError(
                        "Missing predictions, example peptide='%s' allele='%s'"
                        % (peptide, allele))
                results[index] = PeptideResult(preds=tuple(peptide_result.preds))
        return results

    def predict_pairs_dataframe(self, peptides, alleles=None, sample_name=""):
        """``predict_pairs()`` flattened to a DataFrame."""
        import pandas as pd
        from .pred import COLUMNS

        dfs = [
            pp.to_dataframe(sample_name)
            for pp in self.predict_pairs(peptides, alleles)]
        if not dfs:
            return pd.DataFrame(columns=COLUMNS)
        return pd.concat(dfs, ignore_index=True)

    def predict(self, peptides, n_flanks=None, c_flanks=None):
        """
        Predict for a list of peptide sequences.

        Returns list of PeptideResult. When a native parse_to_preds_fn is
        available, parses directly to Pred objects. Otherwise falls back
        to converting from BindingPrediction.
        """
        peptides, _, _ = self._check_flank_inputs(
            peptides, n_flanks, c_flanks)
        return self._predict_for_alleles(
            peptides=peptides,
            alleles=self.alleles,
            allele_cli_names=getattr(self, "_allele_cli_names", None))

    def predict_peptides(self, peptides):
        return self._predict_binding_predictions_for_alleles(
            peptides=peptides,
            alleles=self.alleles,
            allele_cli_names=getattr(self, "_allele_cli_names", None))
