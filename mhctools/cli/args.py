
# Copyright (c) 2016. Mount Sinai School of Medicine
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""
Commandline options for MHC Binding Prediction
"""

from argparse import Action
from .help import HelpArgumentParser as ArgumentParser
import inspect
import logging


from .parsing_helpers import parse_int_list


class _LazyPredictor:
    """Lazy-import factory for CLI predictors.

    NOTE: this is an instance, not a class. `isinstance(x, _LazyPredictor(...))`
    and `issubclass(cls, _LazyPredictor(...))` will NOT work against these
    factories. For type checks use `from mhctools import BigMHC` (or
    `from mhctools.bigmhc import BigMHC` / `_resolve()`). Any of those
    paths will eagerly load the underlying heavy dependency (torch/TF),
    which is unavoidable when you need the real class object — the whole
    point of `_LazyPredictor` is to avoid that load until it's needed.
    """
    def __init__(self, module_name, class_name):
        self._module_name = module_name
        self._class_name = class_name
        self._cls = None
        self.__name__ = class_name

    def _resolve(self):
        if self._cls is None:
            from importlib import import_module
            module = import_module(self._module_name, package="mhctools")
            self._cls = getattr(module, self._class_name)
        return self._cls

    def __call__(self, *args, **kwargs):
        return self._resolve()(*args, **kwargs)

    def __getattr__(self, item):
        return getattr(self._resolve(), item)

    def __eq__(self, other):
        if other is self:
            return True
        if isinstance(other, _LazyPredictor):
            return (self._module_name, self._class_name) == (
                other._module_name, other._class_name)
        # Only classes can meaningfully equal a lazy-loaded class reference.
        # Bail out without triggering the import for anything else.
        if not isinstance(other, type):
            return NotImplemented
        return self._resolve() is other

    def __hash__(self):
        return hash((self._module_name, self._class_name))


logger = logging.getLogger(__name__)

_PREDICTORS = {
    "netmhc": _LazyPredictor(".netmhc", "NetMHC"),
    "netmhc3": _LazyPredictor(".netmhc3", "NetMHC3"),
    "netmhc4": _LazyPredictor(".netmhc4", "NetMHC4"),
    "netmhcpan": _LazyPredictor(".netmhc_pan", "NetMHCpan"),
    "netmhcpan4": _LazyPredictor(".netmhc_pan4", "NetMHCpan4"),
    "netmhcpan4-ba": _LazyPredictor(".netmhc_pan4", "NetMHCpan4_BA"),
    "netmhcpan4-el": _LazyPredictor(".netmhc_pan4", "NetMHCpan4_EL"),
    "netmhcpan41": _LazyPredictor(".netmhc_pan41", "NetMHCpan41"),
    "netmhcpan41-ba": _LazyPredictor(".netmhc_pan41", "NetMHCpan41_BA"),
    "netmhcpan41-el": _LazyPredictor(".netmhc_pan41", "NetMHCpan41_EL"),
    "netmhcpan42": _LazyPredictor(".netmhc_pan42", "NetMHCpan42"),
    "netmhcpan42-ba": _LazyPredictor(".netmhc_pan42", "NetMHCpan42_BA"),
    "netmhcpan42-el": _LazyPredictor(".netmhc_pan42", "NetMHCpan42_EL"),
    "netmhcpan28": _LazyPredictor(".netmhc_pan28", "NetMHCpan28"),
    "netmhcpan3": _LazyPredictor(".netmhc_pan3", "NetMHCpan3"),
    "netmhciipan": _LazyPredictor(".netmhcii_pan", "NetMHCIIpan"),
    "netmhciipan3": _LazyPredictor(".netmhcii_pan", "NetMHCIIpan3"),
    "netmhciipan4": _LazyPredictor(".netmhcii_pan", "NetMHCIIpan4"),
    "netmhciipan4-ba": _LazyPredictor(".netmhcii_pan", "NetMHCIIpan4_BA"),
    "netmhciipan4-el": _LazyPredictor(".netmhcii_pan", "NetMHCIIpan4_EL"),
    "netmhciipan43": _LazyPredictor(".netmhcii_pan", "NetMHCIIpan43"),
    "netmhciipan43-ba": _LazyPredictor(".netmhcii_pan", "NetMHCIIpan43_BA"),
    "netmhciipan43-el": _LazyPredictor(".netmhcii_pan", "NetMHCIIpan43_EL"),
    "netmhccons": _LazyPredictor(".netmhc_cons", "NetMHCcons"),
    "netmhcstabpan": _LazyPredictor(".netmhcstabpan", "NetMHCstabpan"),
    "bigmhc": _LazyPredictor(".bigmhc", "BigMHC"),
    "bigmhc-el": _LazyPredictor(".bigmhc", "BigMHC_EL"),
    "bigmhc-im": _LazyPredictor(".bigmhc", "BigMHC_IM"),
    "caphla": _LazyPredictor(".caphla", "CapHLA"),
    "caphla-ba": _LazyPredictor(".caphla", "CapHLA_BA"),
    "caphla-el": _LazyPredictor(".caphla", "CapHLA_EL"),
    "netchop": _LazyPredictor(".netchop", "NetChop"),
    "netcleave": _LazyPredictor(".netcleave", "NetCleave"),
    "netcleave-i": _LazyPredictor(".netcleave", "NetCleave_I"),
    "netcleave-ii": _LazyPredictor(".netcleave", "NetCleave_II"),
    "pepsickle": _LazyPredictor(".pepsickle", "Pepsickle"),
    "random": _LazyPredictor(".random_predictor", "RandomBindingPredictor"),
    "smm": _LazyPredictor(".smm", "SMM"),
    "smm-pmbec": _LazyPredictor(".smm", "SMMPMBEC"),
    "netmhcpan-iedb": _LazyPredictor(".iedb", "IedbNetMHCpan"),
    "netmhccons-iedb": _LazyPredictor(".iedb", "IedbNetMHCcons"),
    "smm-iedb": _LazyPredictor(".iedb", "IedbSMM"),
    "smm-pmbec-iedb": _LazyPredictor(".iedb", "IedbSMM_PMBEC"),
    "netmhciipan-iedb": _LazyPredictor(".iedb", "IedbNetMHCIIpan"),
    "mhcflurry": _LazyPredictor(".mhcflurry", "MHCflurry"),
    "mhcflurry-affinity": _LazyPredictor(".mhcflurry", "MHCflurry_Affinity"),
    "mixmhcpred": _LazyPredictor(".mixmhcpred", "MixMHCpred"),
    "mixmhc2pred": _LazyPredictor(".mixmhc2pred", "MixMHC2pred"),
    "prime": _LazyPredictor(".prime", "PRIME"),
    "deeptap": _LazyPredictor(".deeptap", "DeepTAP"),
    "calis": _LazyPredictor(".calis", "Calis"),
    "weirdo": _LazyPredictor(".weirdo", "WEIRDO"),
    "eramer": _LazyPredictor(".eramer", "ERAMER"),
    "deepimmuno": _LazyPredictor(".deepimmuno", "DeepImmuno"),
    "tlimmuno2": _LazyPredictor(".tlimmuno2", "TLimmuno2"),
    "peptiverse": _LazyPredictor(".peptiverse", "PeptiVerse"),
    "peptiverse-cpp": _LazyPredictor(".peptiverse_cpp", "PeptiVerseCPP"),
    "plifepred2": _LazyPredictor(".plifepred2", "PlifePred2"),
    "cavaco": _LazyPredictor(".cavaco", "CavacoHalfLife"),
}


# Keep the original public dict contract: ordinary entries are actual classes
# or factories, while the existing heavyweight backends remain lazy proxies.
# Parser construction needs only names and must not initialize that registry.
_HEAVY_PREDICTORS = frozenset((
    "bigmhc", "bigmhc-el", "bigmhc-im", "caphla", "caphla-ba", "caphla-el",
    "mhcflurry", "mhcflurry-affinity",
))


def _predictor_registry():
    if "mhc_predictors" not in globals():
        globals()["mhc_predictors"] = {
            name: factory if name in _HEAVY_PREDICTORS else factory._resolve()
            for name, factory in _PREDICTORS.items()
        }
    return globals()["mhc_predictors"]


def __getattr__(name):
    if name == "mhc_predictors":
        return _predictor_registry()
    raise AttributeError("module %r has no attribute %r" % (__name__, name))


def __dir__():
    return sorted(set(globals()) | {"mhc_predictors"})


def _predictor_names():
    # Honor caller changes to the public registry after it has been requested.
    return globals().get("mhc_predictors", _PREDICTORS).keys()


def _cls_accepts(cls, param_name):
    """Check whether a predictor class/factory accepts a given __init__ parameter."""
    # Unwrap lazy predictors so inspect sees the real class.
    if isinstance(cls, _LazyPredictor):
        cls = cls._resolve()
    sig = inspect.signature(cls)
    if any(p.kind == inspect.Parameter.VAR_KEYWORD for p in sig.parameters.values()):
        return True
    return param_name in sig.parameters


def _parse_predictor_names(value):
    """Normalize a single CLI token: lowercase, split on commas, validate."""
    names = [s.strip().lower() for s in value.split(",") if s.strip()]
    for name in names:
        if name not in _predictor_names():
            from argparse import ArgumentTypeError
            raise ArgumentTypeError(
                "Unknown predictor %r. Available: %s"
                % (name, ", ".join(sorted(_predictor_names()))))
    return names


def _complete_predictor_names(prefix, **kwargs):
    head, separator, tail = prefix.rpartition(",")
    return [head + separator + name for name in sorted(_predictor_names())
            if name.startswith(tail.lower())]


class _StoreAlleleArguments(Action):
    """Accept multiple arguments while retaining the embedded CLI's string API."""

    def __call__(self, parser, namespace, values, option_string=None):
        setattr(namespace, self.dest, " ".join(values))


def add_mhc_args(arg_parser, *, predictor_required=True):
    """Add the shared predictor options to an application's parser.

    Parameters
    ----------
    arg_parser : argparse.ArgumentParser
        Parser to extend.
    predictor_required : bool, optional
        Require --mhc-predictor during parsing (default True). Set False for
        applications that can replay saved predictions. Calling
        predictors_from_args still requires an explicit predictor selection.
    """
    mhc_options_arg_group = arg_parser.add_argument_group(
        title="Prediction Options",
        description="Model selection and prediction options")

    mhc_options_arg_group.add_argument(
        "--mhc-predictor",
        nargs="+",
        type=_parse_predictor_names,
        required=predictor_required,
        help=(
            "One or more predictor names, space or comma separated. "
            "E.g. '--mhc-predictor netmhcpan42' or "
            "'--mhc-predictor netmhcpan42 bigmhc-el pepsickle' or "
            "'--mhc-predictor netmhcpan42,bigmhc-el,pepsickle'. "
            "Available: %s" % ", ".join(sorted(_predictor_names()))
        )).completer = _complete_predictor_names

    mhc_options_arg_group.add_argument(
        "--mhc-predictor-path",
        help="Path to executable program or model directory",
        default=None,
        required=False)

    mhc_options_arg_group.add_argument(
        "--mhc-predictor-models-path",
        help="Path to models to use with predictor (e.g. mhcflurry).")

    mhc_options_arg_group.add_argument(
        "--mixmhcpred-python",
        help="Isolated pandas<3 Python executable for MixMHCpred or PRIME; "
             "overrides MIXMHCPRED_PYTHON without changing host dependencies.")

    mhc_options_arg_group.add_argument(
        "--mhc-peptide-lengths",
        type=parse_int_list,
        help="Peptide lengths to predict, e.g. '9' or '8-11' or '8,9,10'. "
             "Defaults to each predictor's own default lengths, which are "
             "narrower than the lengths it supports (the NetMHC family "
             "defaults to 9 alone).")

    mhc_options_arg_group.add_argument(
        "--mhc-epitope-lengths",
        type=parse_int_list,
        help="Deprecated name for --mhc-peptide-lengths")

    mhc_options_arg_group.add_argument(
        "--mhc-alleles-file",
        help="File with one HLA allele per line")

    mhc_options_arg_group.add_argument(
        "--mhc-alleles",
        nargs="+",
        action=_StoreAlleleArguments,
        metavar="ALLELE",
        default="",
        help="One or more alleles, as separate quoted values or a quoted "
             "comma/space-separated list. Quote '*' to prevent shell expansion. "
             "E.g. --mhc-alleles 'HLA-A*02:01' 'HLA-B*07:02'. Required (here or via "
             "--mhc-alleles-file) for every allele-specific predictor; "
             "there is no default allele.")

    return mhc_options_arg_group


def make_mhc_arg_parser(**kwargs):
    parser = ArgumentParser(**kwargs)
    add_mhc_args(parser)
    return parser

def mhc_alleles_from_args(args):
    from ..allele_normalization import normalize_allele_name_or_raw

    # normalize_allele_name_or_raw keeps exotic / non-human alleles that
    # mhcgnomes can't parse (e.g. BoLA-amani.1, H-2-Qa1) instead of raising,
    # matching the commandline predictors' keep-unparseable behavior (#220).
    alleles = [
        normalize_allele_name_or_raw(allele.strip())
        for comma_group in args.mhc_alleles.split(",")
        for allele in comma_group.split(" ")
        if allele.strip()
    ]
    if args.mhc_alleles_file:
        with open(args.mhc_alleles_file, 'r') as f:
            for line in f:
                line = line.strip()
                if line:
                    alleles.append(normalize_allele_name_or_raw(line))
    if len(alleles) == 0:
        raise ValueError(
            "MHC alleles required (use --mhc-alleles or --mhc-alleles-file)")
    return alleles

def _flatten_predictor_names(args):
    """Flatten the nested list from nargs="+" + comma-splitting type function."""
    raw = args.mhc_predictor or []  # list of lists, or None when optional
    names = [name for group in raw for name in group]
    if not names:
        raise ValueError("MHC predictor required (use --mhc-predictor)")
    return names

def _build_predictor(cls, name, alleles, peptide_lengths, args):
    """Build a single predictor instance, passing only kwargs it accepts."""
    kwargs = {}
    if alleles is not None and _cls_accepts(cls, "alleles"):
        kwargs["alleles"] = alleles
    if peptide_lengths is not None and _cls_accepts(cls, "default_peptide_lengths"):
        kwargs["default_peptide_lengths"] = peptide_lengths
    if getattr(args, "mhc_predictor_models_path", None) and _cls_accepts(cls, "models_path"):
        kwargs["models_path"] = args.mhc_predictor_models_path
    if getattr(args, "mhc_predictor_path", None):
        if _cls_accepts(cls, "bigmhc_path"):
            kwargs["bigmhc_path"] = args.mhc_predictor_path
        elif _cls_accepts(cls, "caphla_path"):
            kwargs["caphla_path"] = args.mhc_predictor_path
        elif _cls_accepts(cls, "program_name"):
            kwargs["program_name"] = args.mhc_predictor_path
    if getattr(args, "mixmhcpred_python", None):
        if name == "mixmhcpred":
            kwargs["python_executable"] = args.mixmhcpred_python
        elif name == "prime":
            kwargs["mixmhcpred_python"] = args.mixmhcpred_python
    logger.info("Building predictor %s(%s)",
                getattr(cls, "__name__", name), kwargs)
    return cls(**kwargs)


def predictors_from_args(args):
    """Build predictor instances from --mhc-predictor.

    Returns a list of instantiated predictor objects.
    """
    names = _flatten_predictor_names(args)
    if getattr(args, "mixmhcpred_python", None) and not set(names) & {"mixmhcpred", "prime"}:
        raise ValueError("--mixmhcpred-python requires MixMHCpred or PRIME selection")

    mhc_predictors = _predictor_registry()

    # Fetch alleles only if at least one predictor needs them
    needs_alleles = any(_cls_accepts(mhc_predictors[n], "alleles") for n in names)
    alleles = None
    if needs_alleles:
        alleles = mhc_alleles_from_args(args)

    peptide_lengths = getattr(args, "mhc_peptide_lengths", None)
    if not peptide_lengths:
        peptide_lengths = getattr(args, "mhc_epitope_lengths", None)

    return [
        _build_predictor(mhc_predictors[name], name, alleles, peptide_lengths, args)
        for name in names
    ]


def mhc_binding_predictor_from_args(args):
    """Build a single predictor from --mhc-predictor (legacy entry point)."""
    names = _flatten_predictor_names(args)
    if len(names) != 1:
        raise ValueError(
            "mhc_binding_predictor_from_args expects exactly one predictor, "
            "got %d. Use predictors_from_args for multiple." % len(names))
    return predictors_from_args(args)[0]
