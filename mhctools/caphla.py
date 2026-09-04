# Copyright (c) 2026 Mount Sinai School of Medicine
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0

"""Thin in-process wrapper for the CapHLA EL and BA model ensembles.

CapHLA predicts peptide presentation (EL) and binding affinity (BA) for human
and mouse MHC class I and II. Its MIT-licensed model definitions and five-fold
weights are fetched as a pinned upstream snapshot; mhctools supplies an
import-safe API, deterministic in-process batching, allele normalization, and
canonical :class:`~mhctools.pred.Prediction` outputs.

Upstream: https://github.com/changyunjian/CapHLA
Cite: Chang and Wu, Briefings in Bioinformatics 2025, bbae595.
"""

import csv
from importlib import util
import os
from pathlib import Path

from mhcgnomes import Allele, Gene, Pair, parse

from .allele_normalization import (
    AlleleParseError,
    normalize_allele_name,
    parse_classi_or_classii_allele_name,
)
from .pred import COLUMNS, Kind, PeptideResult, Prediction
from .wrapper_base import NewModelPredictorMixin


CAPHLA_REVISION = "33ebdd6ce6dadbbb1c66b026ce4b5d81dbf3a831"
CAPHLA_PEPTIDE_LENGTHS = tuple(range(7, 26))
_AMINO_ACIDS = "ACDEFGHIKLMNPQRSTVWYX"
_AMINO_ACID_INDEX = {amino_acid: i for i, amino_acid in enumerate(_AMINO_ACIDS)}


def _find_caphla_home(caphla_path=None):
    """Resolve and validate a CapHLA upstream snapshot."""
    for source, path in (
            ("caphla_path argument", caphla_path),
            ("CAPHLA_HOME", os.environ.get("CAPHLA_HOME"))):
        if path:
            candidate = Path(path).expanduser().resolve()
            if not (candidate / "EL_model.py").is_file():
                raise FileNotFoundError(
                    "CapHLA directory from %s is invalid: %s" % (
                        source, candidate))
            return str(candidate)

    legacy = Path("~/CapHLA").expanduser()
    if (legacy / "EL_model.py").is_file():
        return str(legacy.resolve())

    from .artifacts import artifact_status
    managed = artifact_status("caphla")
    if managed.manager == "mhctools" and managed.status == "ready":
        return managed.path
    raise FileNotFoundError(
        "CapHLA not found. Run `mhctools fetch caphla`, set CAPHLA_HOME, "
        "or pass caphla_path= to the constructor.")


def _component_name(component, include_species):
    core = component.gene
    if component.allele_family:
        core += "*%s:%s" % (
            component.allele_family, component.allele_code)
    else:
        core += component.allele_code
    if include_species:
        species = "H2" if component.species == "H-2" else component.species
        return "%s-%s" % (species, core)
    return core


def _annotated_component_name(component, include_species):
    gene = component.gene_name
    if component.species_prefix == "HLA" and gene == "DRA":
        gene = "DRA1"
    core = gene + "*" + ":".join(component.allele_fields)
    core += "".join(component.annotations)
    if include_species:
        species = (
            "H-2" if component.species_prefix == "H2"
            else component.species_prefix)
        return "%s-%s" % (species, core)
    return core


def _mhcgnomes_fallback_names(parsed):
    """Format expression variants and gene-level aliases for CapHLA."""
    if isinstance(parsed, Allele):
        canonical = _annotated_component_name(parsed, include_species=True)
        return canonical, parsed.to_string()
    if isinstance(parsed, Pair):
        components = (parsed.alpha, parsed.beta)
        canonical = "%s-%s" % (
            _annotated_component_name(components[0], include_species=True),
            _annotated_component_name(components[1], include_species=False),
        )
        # DRA is implicit in CapHLA's DR keys.
        if components[0].gene_name == "DRA":
            key = components[1].to_string()
        else:
            key = "%s/%s" % (
                components[0].to_string(),
                _annotated_component_name(
                    components[1], include_species=False),
            )
        return canonical, key
    if isinstance(parsed, Gene):
        # mhcgnomes resolves the Qa aliases to their official gene names;
        # CapHLA's library retains the historical aliases.
        caphla_key_by_gene = {
            "H2-T23": "H2-Qa1",
            "H2-Q7": "H2-Qa2",
        }
        mhcgnomes_name = parsed.to_string()
        try:
            key = caphla_key_by_gene[mhcgnomes_name]
        except KeyError:
            raise AlleleParseError(
                "CapHLA does not support gene-level MHC name %s" %
                mhcgnomes_name)
        canonical = "H-2-%s" % parsed.gene_name
        return canonical, key
    raise AlleleParseError("Unsupported parsed MHC name: %s" % parsed)


def normalize_caphla_allele(allele):
    """Return ``(canonical_name, CapHLA_key, class)`` for an MHC allele.

    Parsing and canonical identity come from mhcgnomes. CapHLA's only naming
    differences are that DRA1 is implicit for DR, DP/DQ chains use ``/``, and
    mouse alleles use the ``H2-`` prefix.
    """
    raw = str(allele).strip()
    try:
        components = parse_classi_or_classii_allele_name(raw)
        canonical = normalize_allele_name(raw)
    except AlleleParseError:
        # The compatibility normalizer intentionally rejects expression
        # suffixes such as N/L. mhcgnomes itself parses them, and CapHLA has
        # explicit pseudosequences for two such alleles, so preserve them.
        parsed = parse(raw, infer_class2_pairing=True, raise_on_error=True)
        canonical, key = _mhcgnomes_fallback_names(parsed)
        mhc_class = "II" if str(parsed.mhc_class).startswith("II") else "I"
        return canonical, key, mhc_class
    parsed = parse(canonical, infer_class2_pairing=True, raise_on_error=True)
    mhc_class = "II" if str(parsed.mhc_class).startswith("II") else "I"

    if len(components) == 1:
        key = _component_name(components[0], include_species=True)
    elif components[0].gene == "DRA1":
        key = _component_name(components[1], include_species=True)
    else:
        key = "%s/%s" % (
            _component_name(components[0], include_species=True),
            _component_name(components[1], include_species=False),
        )
    return canonical, key, mhc_class


def _load_upstream_class(path, module_name, class_name):
    spec = util.spec_from_file_location(module_name, path)
    if spec is None or spec.loader is None:
        raise ImportError("Could not load CapHLA module at %s" % path)
    module = util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return getattr(module, class_name)


def _affinity_nm(score):
    """Invert CapHLA's ``1 - log(IC50) / log(50000)`` BA transform."""
    return 50000.0 ** (1.0 - float(score))


class CapHLA(NewModelPredictorMixin):
    """CapHLA presentation and binding-affinity predictor.

    Parameters
    ----------
    alleles : list of str
        Human HLA or mouse H-2 class-I/class-II alleles.
    mode : str
        ``"el"`` for presentation, ``"ba"`` for affinity, or ``"both"``.
    caphla_path : str, optional
        Manual CapHLA checkout. Otherwise uses ``CAPHLA_HOME``, ``~/CapHLA``,
        or the snapshot installed by :meth:`fetch`.
    device : str
        PyTorch device. Defaults to deterministic CPU inference.
    batch_size : int
        Number of explicit peptide/allele pairs per inference batch.
    """

    VALID_MODES = ("el", "ba", "both")

    @classmethod
    def fetch(cls, version=None, data_dir=None):
        """Fetch the pinned CapHLA source and EL/BA weights."""
        from .artifacts import fetch
        return fetch("caphla", version=version, data_dir=data_dir)

    def __init__(
            self,
            alleles,
            mode="both",
            caphla_path=None,
            device="cpu",
            batch_size=128):
        if isinstance(alleles, str):
            alleles = [alleles]
        if not alleles:
            raise ValueError("CapHLA requires at least one allele")
        if mode not in self.VALID_MODES:
            raise ValueError(
                "mode must be one of %s, got %r" % (self.VALID_MODES, mode))
        if batch_size <= 0:
            raise ValueError("batch_size must be positive")

        self.mode = mode
        self.caphla_home = _find_caphla_home(caphla_path)
        self.device = device
        self.batch_size = batch_size
        self._hla_library = self._load_hla_library()

        normalized = [normalize_caphla_allele(allele) for allele in alleles]
        self.alleles = [item[0] for item in normalized]
        self._allele_keys = [item[1] for item in normalized]
        self._mhc_classes = {item[2] for item in normalized}
        self._validate_allele_keys(self._allele_keys)

        self._torch = None
        self._models = {}

    def __str__(self):
        loaded = "loaded" if self._models else "not loaded"
        return "CapHLA(mode=%s, alleles=%s, %s)" % (
            self.mode, self.alleles, loaded)

    def _active_modes(self):
        return ("el", "ba") if self.mode == "both" else (self.mode,)

    def kind_support(self):
        mhc_class = (
            next(iter(self._mhc_classes))
            if len(self._mhc_classes) == 1 else "both")
        result = {}
        if "el" in self._active_modes():
            result[Kind.pMHC_presentation] = {
                "mhc_dependence": "single_allele",
                "mhc_class": mhc_class,
            }
        if "ba" in self._active_modes():
            result[Kind.pMHC_affinity] = {
                "mhc_dependence": "single_allele",
                "mhc_class": mhc_class,
            }
        return result

    def _load_hla_library(self):
        path = Path(self.caphla_home) / "HLA_library.csv"
        if not path.is_file():
            raise FileNotFoundError("CapHLA HLA_library.csv not found: %s" % path)
        with path.open(newline="") as input_file:
            rows = csv.DictReader(input_file)
            library = {
                row["Allele Name"]: row["MHC pseudo-seq"] for row in rows}
        if not library:
            raise ValueError("CapHLA HLA_library.csv is empty")
        return library

    def _validate_allele_keys(self, allele_keys):
        unsupported = [key for key in allele_keys if key not in self._hla_library]
        if unsupported:
            raise ValueError(
                "Allele(s) unsupported by CapHLA: %s" %
                ", ".join(unsupported))

    def _ensure_loaded(self):
        missing_modes = [
            mode for mode in self._active_modes() if mode not in self._models]
        if not missing_modes:
            return
        try:
            import torch
        except ImportError as error:
            raise ImportError(
                "CapHLA requires PyTorch; install mhctools[caphla]") from error

        definitions = {
            "el": ("EL_model.py", "_mhctools_caphla_el", "CapHLA_EL"),
            "ba": ("BA_model.py", "_mhctools_caphla_ba", "CapHLA_BA"),
        }
        device = torch.device(self.device)
        for mode in missing_modes:
            filename, module_name, class_name = definitions[mode]
            model_class = _load_upstream_class(
                str(Path(self.caphla_home) / filename), module_name, class_name)
            models = []
            for fold in range(5):
                params_path = Path(self.caphla_home) / "params" / (
                    "%s_fold%d.params" % (mode, fold))
                if not params_path.is_file():
                    raise FileNotFoundError(
                        "CapHLA checkpoint not found: %s" % params_path)
                model = model_class().to(device)
                try:
                    state = torch.load(
                        str(params_path), map_location=device, weights_only=True)
                except TypeError:
                    state = torch.load(str(params_path), map_location=device)
                model.load_state_dict(state)
                model.eval()
                models.append(model)
            self._models[mode] = models
        self._torch = torch

    def _encode(self, sequences, max_length):
        indices = []
        for sequence in sequences:
            padded = sequence + "X" * (max_length - len(sequence))
            indices.append([_AMINO_ACID_INDEX[amino_acid] for amino_acid in padded])
        values = self._torch.tensor(indices, dtype=self._torch.long)
        return self._torch.nn.functional.one_hot(
            values, num_classes=len(_AMINO_ACIDS))

    def _check_peptides(self, peptides):
        for peptide in peptides:
            if len(peptide) not in CAPHLA_PEPTIDE_LENGTHS:
                raise ValueError(
                    "CapHLA supports peptide lengths 7-25, got %d-mer %r" % (
                        len(peptide), peptide))
            invalid = set(peptide) - set(_AMINO_ACIDS[:-1])
            if invalid:
                raise ValueError(
                    "CapHLA peptide %r has invalid amino acid(s): %s" % (
                        peptide, "".join(sorted(invalid))))

    def _predict_raw(self, peptides, allele_keys):
        self._ensure_loaded()
        outputs = {mode: [] for mode in self._active_modes()}
        device = self._torch.device(self.device)
        with self._torch.no_grad():
            for start in range(0, len(peptides), self.batch_size):
                stop = start + self.batch_size
                peptide_tensor = self._encode(peptides[start:stop], 25).to(device)
                mhc_sequences = [
                    self._hla_library[key] for key in allele_keys[start:stop]]
                mhc_tensor = self._encode(mhc_sequences, 34).to(device)
                if "el" in outputs:
                    fold_scores = [
                        self._torch.softmax(
                            model(peptide_tensor, mhc_tensor), dim=1)[:, 1]
                        for model in self._models["el"]
                    ]
                    scores = self._torch.stack(fold_scores).mean(dim=0)
                    outputs["el"].extend(scores.cpu().tolist())
                if "ba" in outputs:
                    fold_scores = [
                        model(peptide_tensor, mhc_tensor)
                        for model in self._models["ba"]
                    ]
                    scores = self._torch.stack(fold_scores).mean(dim=0)
                    outputs["ba"].extend(scores.cpu().tolist())
        return outputs

    def predict_pairs(self, peptides, alleles=None):
        """Predict explicit peptide/allele pairs, preserving order and duplicates."""
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
            peptide_list = self._normalize_peptides(peptides)
            allele_list = [alleles] if isinstance(alleles, str) else list(alleles)
            if len(peptide_list) != len(allele_list):
                raise ValueError(
                    "peptides length %d != alleles length %d" % (
                        len(peptide_list), len(allele_list)))
        peptide_list = self._normalize_peptides(peptide_list)
        if not peptide_list:
            return []
        self._check_peptides(peptide_list)

        normalized = [normalize_caphla_allele(allele) for allele in allele_list]
        canonical_alleles = [item[0] for item in normalized]
        allele_keys = [item[1] for item in normalized]
        self._validate_allele_keys(allele_keys)
        outputs = self._predict_raw(peptide_list, allele_keys)

        results = []
        for index, (peptide, allele) in enumerate(zip(
                peptide_list, canonical_alleles)):
            predictions = []
            if "el" in outputs:
                predictions.append(Prediction(
                    kind=Kind.pMHC_presentation,
                    score=float(outputs["el"][index]),
                    peptide=peptide,
                    allele=allele,
                    predictor_name="caphla_el",
                    predictor_version=CAPHLA_REVISION,
                ))
            if "ba" in outputs:
                score = float(outputs["ba"][index])
                predictions.append(Prediction(
                    kind=Kind.pMHC_affinity,
                    score=score,
                    value=_affinity_nm(score),
                    peptide=peptide,
                    allele=allele,
                    predictor_name="caphla_ba",
                    predictor_version=CAPHLA_REVISION,
                ))
            results.append(PeptideResult(preds=tuple(predictions)))
        return results

    def predict_pairs_dataframe(self, peptides, alleles=None, sample_name=""):
        """Flatten :meth:`predict_pairs` to one canonical row per output kind."""
        import pandas as pd
        frames = [
            result.to_dataframe(sample_name)
            for result in self.predict_pairs(peptides, alleles)]
        if not frames:
            return pd.DataFrame(columns=COLUMNS)
        return pd.concat(frames, ignore_index=True)

    def predict(self, peptides):
        """Predict each peptide against every configured allele."""
        peptide_list = self._normalize_peptides(peptides)
        pairs = [
            (peptide, allele)
            for peptide in peptide_list
            for allele in self.alleles
        ]
        flat_results = self.predict_pairs(pairs)
        outputs_per_pair = len(self._active_modes())
        results = []
        index = 0
        for _ in peptide_list:
            predictions = []
            for _ in self.alleles:
                predictions.extend(flat_results[index].preds)
                index += 1
            if len(predictions) != len(self.alleles) * outputs_per_pair:
                raise RuntimeError("CapHLA returned an unexpected number of outputs")
            results.append(PeptideResult(preds=tuple(predictions)))
        return results


class CapHLA_EL(CapHLA):
    """CapHLA presentation-only ensemble."""

    def __init__(
            self, alleles, caphla_path=None, device="cpu", batch_size=128):
        CapHLA.__init__(
            self,
            alleles=alleles,
            mode="el",
            caphla_path=caphla_path,
            device=device,
            batch_size=batch_size)


class CapHLA_BA(CapHLA):
    """CapHLA binding-affinity-only ensemble."""

    def __init__(
            self, alleles, caphla_path=None, device="cpu", batch_size=128):
        CapHLA.__init__(
            self,
            alleles=alleles,
            mode="ba",
            caphla_path=caphla_path,
            device=device,
            batch_size=batch_size)
