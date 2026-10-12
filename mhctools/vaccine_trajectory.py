"""Conditional single-target SLP/mRNA trajectories through APC pMHC loading.

All numerical rates are supplied, not inferred from native predictor scores.
The linear system tracks expected copies, with translation as a birth process.
Local, node and systemic APCs and extracellular lymph/blood are distinct.
"""

from dataclasses import asdict, dataclass
import math
import json
from numbers import Integral
from pathlib import Path

from mhcgnomes import parse
import numpy as np

from .serum_contributions import EnzymeCutRate
from .serum_degradation import DegradationTarget, _input, _nonnegative


EXTRACELLULAR = ("interstitium", "afferent_lymph", "node_fluid", "blood")
REGIONS = ("local", "node", "systemic")
NODE_ANTIGEN = tuple("node_" + phase for phase in (
    "endosome", "cytosol", "er", "protein", "secretory_er", "loaded", "surface"))
NODE_ENTRY_ROUTES = ("node_fluid_uptake", "local_apc_migration", "node_apc_translation", "cross_dressing")
EXCLUSIONS = (
    "Whole-body PBPK, renal/hepatic compartments and recirculation are not resolved",
    "Sequence-derived LNP delivery, endosomal escape and secretion are not predicted",
    "Signal-peptide-derived targets and intramembrane signal-peptide processing are not modeled",
    "Active enzyme exposure, folding, glycosylation and protein binding are not predicted",
    "APC subsets, changing inflammation and cell-number dynamics are not resolved",
    "MHC saturation, competing ligands and recycling are not resolved",
    "Cross-dressing is an effective MHC-I transfer rate, not a prediction of donor compatibility or contact",
    "Membrane-anchored antigens, cell debris and extracellular vesicles are not separately resolved",
    "T-cell activation, costimulation, antibodies and clinical efficacy are not predicted",
)


def vaccine_trajectory_evidence():
    """Read the primary-source mechanism inventory; no fitted rate defaults."""
    return json.loads((Path(__file__).parent / "data" / "vaccine_trajectory_evidence.json").read_text())


@dataclass(frozen=True)
class VaccineRate:
    """An explicit effective first-order rate, including explicit exclusions.

    ``value_per_hour=None`` means unassessed. A zero is a declared assumption,
    not evidence that the mechanism is biologically absent. Sources must identify
    measurements or the caller's scenario; scores and ranks are not rates.
    """

    value_per_hour: float
    basis: str
    source: str

    def __post_init__(self):
        if self.basis not in ("measured", "assumed", "disabled", "unassessed"):
            raise ValueError("Rate basis must be measured, assumed, disabled or unassessed")
        if not isinstance(self.source, str) or not self.source.strip():
            raise ValueError("Every rate needs a source or explicit assumption")
        if self.value_per_hour is None:
            if self.basis != "unassessed":
                raise ValueError("Missing rates must be unassessed")
        else:
            object.__setattr__(self, "value_per_hour", _nonnegative(self.value_per_hour, "rate"))
            if self.basis == "unassessed" or (self.basis == "disabled" and self.value_per_hour != 0):
                raise ValueError("Rate value contradicts its evidence basis")

    @classmethod
    def from_dict(cls, value):
        return cls(value["value_per_hour"], value["basis"], value["source"])


@dataclass(frozen=True)
class VaccineDose:
    """Documented amount of the modeled starting peptide or carrier-bound RNA.

    ``ug`` requires the molecular weight of the actual starting material in
    g/mol. For RNA this includes the full transcript, not the encoded protein.
    This conversion does not estimate encapsulation, release or bioavailability.
    """

    amount: float
    unit: str
    basis: str
    source: str
    molecular_weight_g_per_mol: float = None

    def __post_init__(self):
        if self.unit not in ("copies", "pmol", "ug"):
            raise ValueError("Dose unit must be copies, pmol or ug")
        object.__setattr__(self, "amount", _nonnegative(self.amount, "dose amount"))
        if self.amount <= 0:
            raise ValueError("Dose amount must be positive")
        if self.basis not in ("measured", "assumed") or not isinstance(self.source, str) or not self.source.strip():
            raise ValueError("Dose requires measured/assumed basis and a source")
        if self.unit == "ug":
            if self.molecular_weight_g_per_mol is None or _nonnegative(self.molecular_weight_g_per_mol, "molecular weight") <= 0:
                raise ValueError("Mass dose requires a positive molecular_weight_g_per_mol for the actual material")
            object.__setattr__(self, "molecular_weight_g_per_mol", float(self.molecular_weight_g_per_mol))
        elif self.molecular_weight_g_per_mol is not None:
            raise ValueError("Molecular weight is only used with a ug dose")
        if not math.isfinite(self.copies):
            raise ValueError("Converted dose exceeds numerical range")

    @property
    def copies(self):
        """Starting copies using the exact SI Avogadro constant."""
        if self.unit == "copies":
            return float(self.amount)
        moles = self.amount * 1e-12 if self.unit == "pmol" else self.amount * 1e-6 / self.molecular_weight_g_per_mol
        return moles * 6.02214076e23


@dataclass(frozen=True)
class RouteStep:
    """One effective biological channel; rate key is its unique name."""

    name: str
    origin: str
    destination: str
    mechanism: str
    gate: str = "always"
    kind: str = "transfer"


@dataclass(frozen=True)
class VaccineCleavageRates:
    """Assessment of all modeled cut channels on one current fragment.

    Empty channels require an explicit no-cut assumption or assay source.
    ``None`` from the callback denotes missing assessment instead.
    """

    channels: tuple
    basis: str
    source: str

    def __post_init__(self):
        if self.basis not in ("measured", "assumed"):
            raise ValueError("Cleavage assessment basis must be measured or assumed")
        if not isinstance(self.source, str) or not self.source.strip():
            raise ValueError("Cleavage assessment requires its source/assumption")
        if not isinstance(self.channels, tuple) or any(not isinstance(c, EnzymeCutRate) for c in self.channels):
            raise TypeError("Cleavage channels must be a tuple of EnzymeCutRate records")


def vaccine_route_steps(delivery, mhc_class):
    """Return the included mechanism graph for local IM/SC delivery.

    Parameters
    ----------
    delivery : str
        ``slp``, ``secreted_mrna`` or ``cytosolic_mrna``. RNA routes include
        direct APC expression; secretion is specific to secreted constructs.
    mhc_class : str
        ``I`` or ``II``; each target/allele is modeled separately.
    """
    if delivery not in ("slp", "secreted_mrna", "cytosolic_mrna") or mhc_class not in ("I", "II"):
        raise ValueError("Choose slp/secreted_mrna/cytosolic_mrna and MHC class I/II")
    steps = []

    def add(name, origin, destination, mechanism, gate="always", kind="transfer"):
        steps.append(RouteStep(name, origin, destination, mechanism, gate, kind))

    add("lymph_entry", "interstitium", "afferent_lymph", "Entry into afferent lymph")
    add("vascular_absorption", "interstitium", "blood", "Competing vascular absorption")
    add("lymph_transit", "afferent_lymph", "node_fluid", "Drainage into the draining lymph node")
    add("node_outflow", "node_fluid", "blood", "Effective efferent lymph/venous return")
    for compartment, region in (("interstitium", "local"), ("node_fluid", "node"), ("blood", "systemic")):
        add(region + "_antigen_uptake", compartment, region + "_endosome", "APC internalization; not productive loading")
    for compartment in EXTRACELLULAR:
        add(compartment + "_clearance", compartment, "antigen_cleared", "Removal without a modeled target cut")

    for region in REGIONS:
        endosome, cytosol, er = (region + "_" + c for c in ("endosome", "cytosol", "er"))
        loaded, surface = region + "_loaded", region + "_surface"
        if mhc_class == "I":
            add(region + "_cross_escape", endosome, cytosol, "Endosome-to-cytosol cross-presentation")
            add(region + "_tap", cytosol, er, "TAP transport of a C-terminally complete precursor", "tap")
            add(region + "_er_loading", er, loaded, "MHC-I peptide-loading complex and groove availability", "exact")
            add(region + "_vacuolar_loading", endosome, loaded, "Optional TAP-independent vacuolar MHC-I loading", "exact")
            add(region + "_er_loss", er, "antigen_cleared", "ER precursor disposal")
        else:
            add(region + "_autophagy", cytosol, endosome, "Effective cytosolic antigen delivery to endolysosomes")
            add(region + "_ii_loading", endosome, loaded, "MHC-II loading; includes CLIP exchange, HLA-DM editing and groove availability", "class_ii")
        for compartment in (endosome, cytosol):
            add(compartment + "_loss", compartment, "antigen_cleared", "Unproductive intracellular disposal")
        add(region + "_surface_export", loaded, surface, "Loaded pMHC trafficking to APC surface")
        add(region + "_loaded_loss", loaded, "pmhc_lost", "Loaded-complex disposal before surface display")
        add(region + "_surface_loss", surface, "pmhc_lost", "Effective pMHC dissociation, turnover or APC loss")
    for compartment in ("endosome", "cytosol", "loaded", "surface") + (("er",) if mhc_class == "I" else ()):
        add("migration_" + compartment, "local_" + compartment, "node_" + compartment,
            "Antigen/pMHC carried by a migrating local APC into the draining node")

    if delivery != "slp":
        add("carrier_drainage", "rna_carrier", "node_rna_carrier", "LNP-associated RNA drainage to node")
        add("producer_rna_uptake", "rna_carrier", "producer_rna_endosome", "RNA-carrier uptake by local non-APC producer cells")
        add("local_rna_uptake", "rna_carrier", "local_rna_endosome", "RNA-carrier uptake by local APCs")
        add("node_rna_uptake", "node_rna_carrier", "node_rna_endosome", "RNA-carrier uptake by nodal APCs")
        for compartment in ("rna_carrier", "node_rna_carrier"):
            add(compartment + "_loss", compartment, "rna_lost", "Effective extracellular RNA/carrier loss")
        for cell in ("producer", "local", "node"):
            add(cell + "_rna_escape", cell + "_rna_endosome", cell + "_rna", "RNA endosomal escape into cytosol")
            add(cell + "_endosomal_rna_loss", cell + "_rna_endosome", "rna_lost", "RNA disposal before cytosolic delivery")
            add(cell + "_rna_decay", cell + "_rna", "rna_lost", "Loss of translation-competent RNA")
            product = cell + ("_protein" if delivery == "secreted_mrna" else "_cytosol")
            add(cell + "_translation", cell + "_rna", product, "Completed antigen synthesis per active transcript; RNA remains", kind="birth")
            if delivery == "cytosolic_mrna":
                continue
            add(cell + "_signal_entry", cell + "_protein", cell + "_secretory_er",
                "Effective ER translocation and removal of the supplied signal peptide", "signal")
            add(cell + "_secretion", cell + "_secretory_er", "node_fluid" if cell == "node" else "interstitium",
                "ER/Golgi trafficking and secretion of mature antigen")
            add(cell + "_secretory_loss", cell + "_secretory_er", "antigen_cleared", "Unproductive secretory-pathway loss")
            add(cell + "_failed_entry", cell + "_protein", cell + "_cytosol",
                "Failed ER entry/defective products enter cytosolic processing")
            add(cell + "_erad", cell + "_secretory_er", cell + "_cytosol", "ER-associated retrotranslocation for endogenous processing")
            if cell != "producer" and mhc_class == "II":
                add(cell + "_secretory_to_endosome", cell + "_secretory_er", cell + "_endosome", "Conditional intracellular secretory-antigen delivery to endolysosomes")
        add("producer_cytosol_loss", "producer_cytosol", "antigen_cleared", "Unproductive producer-cell antigen disposal")
        add("producer_antigen_release", "producer_cytosol", "interstitium", "Effective release of core-bearing antigen from producer cells; not active secretion")
        if mhc_class == "I":
            add("producer_tap", "producer_cytosol", "producer_er", "Producer-cell TAP transport", "tap")
            add("producer_er_loading", "producer_er", "producer_loaded", "Producer-cell MHC-I loading on the specified allele", "exact")
            add("producer_er_loss", "producer_er", "antigen_cleared", "Producer-cell ER precursor disposal")
            add("producer_surface_export", "producer_loaded", "producer_surface", "Producer-cell pMHC surface export")
            add("producer_loaded_loss", "producer_loaded", "pmhc_lost", "Producer-cell loaded-complex disposal")
            add("producer_surface_loss", "producer_surface", "pmhc_lost", "Producer-cell surface pMHC loss")
            for region in ("local", "node"):
                add("cross_dressing_" + region, "producer_surface", region + "_surface",
                    "Effective transfer of preformed producer pMHC-I to APCs; no new loading")
        migrating = ("rna_endosome", "rna") + (("protein", "secretory_er") if delivery == "secreted_mrna" else ())
        for compartment in migrating:
            add("migration_" + compartment, "local_" + compartment, "node_" + compartment,
                "Local APC migration with RNA or intracellular antigen")
    return tuple(steps)


@dataclass(frozen=True)
class VaccineTrajectoryInput:
    """Inputs for one exact target and allele, with no physiological defaults.

    Coordinates, including ``signal_end``, are zero-based, half-open on the full
    translated construct. An mRNA source amount is transcript copies, whereas
    an SLP source amount is peptide copies. Their outputs have different units.
    """

    peptide: object
    target: DegradationTarget
    allele: str
    mhc_class: str
    delivery: str
    scenario: str
    rates: dict
    initial_copies: float = None
    signal_end: int = None
    tap_max_length: int = None
    class_ii_max_length: int = None
    dose: VaccineDose = None
    administration: str = "unspecified"
    organism: str = None

    def source_copies(self):
        """Resolve explicit dose or legacy copies; otherwise normalize to one."""
        if self.dose is not None:
            if not isinstance(self.dose, VaccineDose):
                raise TypeError("dose must be a VaccineDose")
            if self.initial_copies is not None:
                raise ValueError("Specify dose or initial_copies, not both")
            return self.dose.copies
        return 1.0 if self.initial_copies is None else _nonnegative(self.initial_copies, "initial_copies")

    def validate(self):
        peptide = _input(self.peptide)
        self.target.validate(peptide.sequence)
        if not isinstance(self.scenario, str) or not self.scenario.strip():
            raise ValueError("A nonempty scenario label is required")
        steps = vaccine_route_steps(self.delivery, self.mhc_class)
        identity = parse(self.allele, raise_on_error=False)
        if identity is None or not (identity.is_class1 if self.mhc_class == "I" else identity.is_class2):
            raise ValueError("Allele must parse with MHCgnomes and match the MHC class")
        if self.source_copies() <= 0:
            raise ValueError("initial_copies must be positive")
        if self.administration not in ("IM", "SC", "unspecified"):
            raise ValueError("administration must be IM, SC or unspecified (local injection)")
        if self.organism is not None and (not isinstance(self.organism, str) or not self.organism.strip()):
            raise ValueError("organism must be a nonempty description or None")
        for name in ("signal_end", "tap_max_length", "class_ii_max_length"):
            value = getattr(self, name)
            if value is not None and (isinstance(value, bool) or not isinstance(value, Integral) or value < 1):
                raise ValueError(name + " must be a positive integer")
        if self.delivery == "secreted_mrna":
            if self.signal_end is None or self.signal_end >= len(peptide.sequence):
                raise ValueError("Secreted mRNA needs the full construct and an explicit signal_end")
            if self.target.start < self.signal_end:
                raise ValueError("Target must lie entirely in the mature antigen; signal-peptide-derived presentation is not modeled")
        elif self.signal_end is not None:
            raise ValueError("Only secreted_mrna input can declare a signal peptide")
        if self.mhc_class == "I" and (self.tap_max_length is None or self.tap_max_length < self.target.end - self.target.start):
            raise ValueError("MHC-I requires an explicit TAP precursor length bound covering the target")
        if self.mhc_class == "II" and (self.class_ii_max_length is None or self.class_ii_max_length < self.target.end - self.target.start):
            raise ValueError("MHC-II requires an explicit loadable ligand length bound covering the target")
        unknown = set(self.rates) - {step.name for step in steps}
        if unknown:
            raise ValueError("Unknown kinetic rate keys: " + ", ".join(sorted(unknown)))
        if any(not isinstance(rate, VaccineRate) for rate in self.rates.values()):
            raise TypeError("rates must contain VaccineRate records")
        return peptide

    @classmethod
    def from_dict(cls, value):
        if value.get("schema_version") != 1:
            raise ValueError("Trajectory schema_version must be 1")
        from .peptide_input import PeptideInput
        peptide = value["peptide"]
        result = cls(
            peptide=PeptideInput.from_dict(peptide) if isinstance(peptide, dict) else peptide,
            target=DegradationTarget(**value["target"]), allele=value["allele"],
            mhc_class=value["mhc_class"], delivery=value["delivery"], scenario=value["scenario"],
            rates={name: VaccineRate.from_dict(rate) for name, rate in value.get("rates", {}).items()},
            initial_copies=value.get("initial_copies"), signal_end=value.get("signal_end"),
            tap_max_length=value.get("tap_max_length"), class_ii_max_length=value.get("class_ii_max_length"),
            dose=VaccineDose(**value["dose"]) if value.get("dose") is not None else None,
            administration=value.get("administration", "unspecified"), organism=value.get("organism"))
        result.validate()
        return result


def _rna(compartment):
    return "rna" in compartment


def _cuts_apply(compartment):
    return compartment in EXTRACELLULAR or compartment in ("producer_secretory_er", "producer_cytosol", "producer_er") or any(
        compartment == region + "_" + phase
        for region in REGIONS for phase in ("endosome", "cytosol", "er", "secretory_er"))


def _destination(model, step, start, end):
    target = model.target
    if step.gate == "exact" and (start, end) != (target.start, target.end):
        return None
    if step.gate == "tap" and (end != target.end or end - start > model.tap_max_length):
        return None
    if step.gate == "class_ii" and end - start > model.class_ii_max_length:
        return None
    if step.gate == "signal":
        if start != 0:
            raise ValueError("Signal processing requires the full translated N-terminus")
        start = model.signal_end
    if step.destination in ("antigen_cleared", "rna_lost", "pmhc_lost"):
        return step.destination, 0, 0
    if _rna(step.destination):
        return step.destination, 0, 0
    if step.kind == "birth":
        start, end = 0, len(_input(model.peptide).sequence)
    return step.destination, start, end


def _advance(values, diagonal, channels, duration):
    """Exponential action by nonnegative uniformization, including RNA births.

    No fitted solver step size: the Poisson series computes exp(A*t). A birth
    increases antigen count without subtracting RNA. Cumulative counters are
    additional linear states. Sparse scatter avoids a dense fragment matrix.
    """
    outgoing = np.zeros_like(values)
    for origin, _, rate in channels:
        outgoing[origin] += rate
    if not np.all(np.isfinite(outgoing)) or not np.all(np.isfinite(diagonal)):
        raise ValueError("Combined rates exceed numerical range")
    scale = max(float(np.max(outgoing)), float(np.max(-diagonal)))
    if not scale or not duration:
        return values.copy()
    chunks = max(1, math.ceil(scale * duration / 4))
    if chunks > 100000:
        raise ValueError("Rate/time combination too large for this scenario solver")
    poisson_time = scale * duration / chunks
    origin = np.array([channel[0] for channel in channels], dtype=int)
    destination = np.array([channel[1] for channel in channels], dtype=int)
    weight = np.array([channel[2] / scale for channel in channels])
    stay = 1 + diagonal / scale
    current = values.copy()
    for _ in range(chunks):
        term = current * math.exp(-poisson_time)
        result = term.copy()
        for order in range(1, 160):
            scattered = np.bincount(destination, weights=weight * term[origin], minlength=len(values))
            term = (stay * term + scattered) * (poisson_time / order)
            result += term
            if order > 20 and np.max(np.abs(term)) <= 1e-13 * max(1.0, float(np.max(result))):
                break
        else:
            raise ValueError("Trajectory exponential series failed to converge")
        if not np.all(np.isfinite(result)):
            raise ValueError("Trajectory amounts exceed numerical range")
        current = result
    return current


def _node_entry_system(states, events, initial_copies):
    """Absorb each antigen lineage on its first entry into a nodal APC.

    RNA remains a repeatable source, including inside nodal APCs. Absorbed
    antigen cannot be recounted on secretion and reuptake. This observer is a
    separate linear system and does not remove antigen from the main trajectory.
    """
    size = len(states)
    values = np.zeros(size + len(NODE_ENTRY_ROUTES))
    values[0] = initial_copies
    diagonal = np.zeros_like(values)
    channels = []
    for origin, destination, rate, birth, name in events:
        if states[origin][0] in NODE_ANTIGEN:
            continue
        if not birth:
            diagonal[origin] -= rate
        if states[destination][0] in NODE_ANTIGEN:
            if name == "node_antigen_uptake":
                route = "node_fluid_uptake"
            elif name == "node_translation":
                route = "node_apc_translation"
            elif name == "cross_dressing_node":
                route = "cross_dressing"
            elif name and name.startswith("migration_"):
                route = "local_apc_migration"
            else:
                raise ValueError("Unclassified node APC entry: %s" % name)
            destination = size + NODE_ENTRY_ROUTES.index(route)
        channels.append((origin, destination, rate))
    return values, diagonal, channels


def simulate_vaccine_trajectory(model, times_hours, cleavage_rates=None, *, max_states=10000):
    """Follow a target to APC loading using explicitly supplied rate channels.

    Parameters
    ----------
    model : VaccineTrajectoryInput
        One exact target occurrence, allele, local delivery and rate audit.
    times_hours : iterable of float
        Finite nonnegative observation times; output preserves their order.
    cleavage_rates : callable, optional
        Called with ``(compartment, TargetFragment)`` for each reachable free
        fragment. Return VaccineCleavageRates using CURRENT-fragment bonds,
        or None when kinetics are missing. Empty channels need their own explicit
        no-cut source, not proof of stability. A serum predictor or score must
        never be silently reused in lymph.
    max_states : int
        Explicit size guard for fragment closure; no silent fragment pruning.

    Returns
    -------
    dict
        Mechanism and rate audit, fragment states and expected-copy time courses.
        Missing reachable kinetics return ``unassessed`` with no numerical curves.
        Otherwise results are conditional scenarios, not calibrated forecasts.
        RNA translation yields multiple antigen copies per input transcript;
        SLP and RNA source-normalized outputs cannot be compared as probabilities.
    """
    if not isinstance(model, VaccineTrajectoryInput):
        raise TypeError("model must be VaccineTrajectoryInput")
    peptide = model.validate()
    initial_copies = model.source_copies()
    times = tuple(_nonnegative(t, "time") for t in times_hours)
    if not times:
        raise ValueError("At least one observation time is required")
    if isinstance(max_states, bool) or not isinstance(max_states, Integral) or max_states < 1:
        raise ValueError("max_states must be a positive integer")
    from .serum_degradation import TargetFragment
    sequence, target = peptide.sequence, model.target
    steps = vaccine_route_steps(model.delivery, model.mhc_class)
    by_origin = {}
    for step in steps:
        by_origin.setdefault(step.origin, []).append(step)
    initial = ("interstitium", 0, len(sequence)) if model.delivery == "slp" else ("rna_carrier", 0, 0)
    states, index, channels, consumed, missing, cut_audit = [initial], {initial: 0}, [], [], set(), []
    events = []

    def state_index(state):
        if state not in index:
            if len(states) >= max_states:
                raise ValueError("Fragment closure exceeds max_states; no states were pruned")
            index[state] = len(states)
            states.append(state)
        return index[state]

    def channel(origin, destination, rate, *, birth=False, counter=None, name=None):
        if not rate:
            return
        destination_index = state_index(destination)
        channels.append((origin, destination_index, rate))
        events.append((origin, destination_index, rate, birth, name))
        if not birth:
            consumed.append((origin, rate))
        if counter:
            channels.append((origin, state_index((counter, 0, 0)), rate))

    cursor = 0
    while cursor < len(states):
        compartment, start, end = states[cursor]
        for step in by_origin.get(compartment, ()):
            destination = _destination(model, step, start, end)
            if destination is None:
                continue
            record = model.rates.get(step.name)
            if record is None or record.value_per_hour is None:
                missing.add("step:" + step.name)
                # Follow possible missing transport to enumerate downstream gaps.
                state_index(destination)
                continue
            if step.kind == "birth":
                counter = "translated"
            elif step.destination.endswith("_loaded") and not step.origin.endswith("_loaded"):
                counter = "producer_loaded_count" if step.destination == "producer_loaded" else "loaded"
            elif step.destination.endswith("_surface") and step.origin.endswith("_loaded"):
                counter = "producer_exported_count" if step.destination == "producer_surface" else "surface_exported"
            else:
                counter = None
            channel(cursor, destination, record.value_per_hour, birth=step.kind == "birth", counter=counter, name=step.name)
        if _cuts_apply(compartment):
            fragment = TargetFragment(sequence[start:end], start, end)
            cuts = None if cleavage_rates is None else cleavage_rates(compartment, fragment)
            if cuts is None:
                missing.add("cuts:%s:%d:%d" % (compartment, start, end))
            else:
                if not isinstance(cuts, VaccineCleavageRates):
                    raise TypeError("cleavage_rates must return VaccineCleavageRates or None")
                cut_audit.append({"compartment": compartment, "start": start, "end": end,
                                  "sequence": fragment.sequence, "channels": [asdict(cut) for cut in cuts.channels],
                                  "basis": cuts.basis, "source": cuts.source})
                for cut in cuts.channels:
                    if cut.bond >= end - start:
                        raise ValueError("Cleavage bond lies outside current fragment")
                    if compartment in tuple(region + "_er" for region in (*REGIONS, "producer")) and cut.bond != 1:
                        raise ValueError("APC ER peptide processing supports N-terminal trimming only")
                    bond = start + cut.bond
                    destination = ("target_destroyed", 0, 0) if target.start < bond < target.end else (
                        compartment, bond if bond <= target.start else start, bond if bond >= target.end else end)
                    channel(cursor, destination, cut.rate_per_hour)
        cursor += 1
    state_records = [{"compartment": c, "start": s, "end": e,
                      "sequence": sequence[s:e] if e > s else None} for c, s, e in states]
    audit = [{**asdict(step), "rate": asdict(model.rates[step.name]) if step.name in model.rates else None}
             for step in steps]
    from . import __version__

    output = {"schema_version": 1, "status": "unassessed" if missing else "conditional_scenario",
              "mhctools_version": __version__, "algorithm": "linear first-order expected copies; uniformization",
              "scenario": model.scenario, "delivery": model.delivery, "administration": model.administration,
              "organism": model.organism, "times_hours": list(times),
              "mhc_class": model.mhc_class, "allele": parse(model.allele).to_string(),
              "target": {**asdict(target), "sequence": sequence[target.start:target.end]},
              "coordinate_system": "zero-based half-open on full input; cut bonds local to current fragment",
              "source_units": "input peptide copies" if model.delivery == "slp" else "input carrier-bound RNA copies",
              "exposure_units": "copy-hours integrated from injection (time zero) to each observation",
              "node_entry_definition": "First presence of each core-bearing antigen lineage in a draining-node APC; includes endogenous production and preformed pMHC transfer",
              "initial_copies": initial_copies, "signal_end": model.signal_end,
              "dose": asdict(model.dose) if model.dose is not None else None,
              "dose_basis": "documented_amount" if model.dose is not None else (
                  "explicit_copies" if model.initial_copies is not None else "one_input_normalization"),
              "peptide_input": peptide.to_dict(),
              "tap_max_length": model.tap_max_length, "class_ii_max_length": model.class_ii_max_length,
              "steps": audit, "cuts": cut_audit, "states": state_records,
              "missing_kinetics": sorted(missing), "exclusions": list(EXCLUSIONS), "curves": [],
              "mechanism_evidence": vaccine_trajectory_evidence(),
              "assumptions": ["Dilute, well-mixed, time-invariant first-order compartments",
                              "One exact target occurrence and at most one loaded pMHC per produced antigen copy",
                              "Producer-cell pMHC uses the specified allele; cross-dressing transfers that complex without new loading",
                              "Translation is repeated antigen production, not RNA consumption",
                              "Each cut retains only the target-bearing daughter; an internal cut destroys this exact target",
                              "Loading competes with cuts; bound target leaves free-peptide cleavage and has its own turnover",
                              "Signal removal/ER entry and secretion are effective stages, not sequence predictions",
                              "APC uptake and MHC binding ranks do not establish loading rates",
                              "Transport, uptake and loading use compartment-wide rates across eligible fragments; cleavage is fragment-specific",
                              "Secretory ER/Golgi antigen cuts require a separately supplied assessment; no convertase activity is inferred",
                              "Systemic APCs pool blood-accessible uptake separately from draining-node APCs"]}
    if missing:
        return output
    # Observers integrate copies without consuming antigen. Keep their copy-hour
    # units separate from the physical/cumulative-copy states exposed to callers.
    observed_compartments = tuple(region + "_surface" for region in REGIONS) + EXTRACELLULAR
    observer_index = {name: len(states) + i for i, name in enumerate(observed_compartments)}
    solver_channels = channels + [(i, observer_index[c], 1.0)
                                  for i, (c, _, _) in enumerate(states) if c in observer_index]
    values = np.zeros(len(states) + len(observed_compartments))
    values[0] = initial_copies
    diagonal = np.zeros_like(values)
    for origin, rate in consumed:
        diagonal[origin] -= rate
    entry_values, entry_diagonal, entry_channels = _node_entry_system(states, events, initial_copies)
    previous, rows = 0.0, {}
    for time in sorted(set(times)):
        values = _advance(values, diagonal, solver_channels, time - previous)
        entry_values = _advance(entry_values, entry_diagonal, entry_channels, time - previous)
        previous = time
        amounts = {state: float(values[i]) for i, state in enumerate(states)}

        def total(predicate):
            return math.fsum(v for (c, s, e), v in amounts.items() if predicate(c, s, e))

        def counter(name):
            return amounts.get((name, 0, 0), 0.0)

        produced = (initial_copies if model.delivery == "slp" else 0.0) + counter("translated")
        retained = total(lambda c, s, e: e > s)
        destroyed, cleared, pmhc_lost = (counter(name) for name in ("target_destroyed", "antigen_cleared", "pmhc_lost"))
        surface = {region: total(lambda c, s, e: c == region + "_surface") for region in REGIONS}
        surface_auc = {region: float(values[observer_index[region + "_surface"]]) for region in REGIONS}
        extracellular_auc = {c: float(values[observer_index[c]]) for c in EXTRACELLULAR}
        node_entry = dict(zip(NODE_ENTRY_ROUTES, map(float, entry_values[len(states):])))
        arrived = math.fsum(node_entry.values())
        blood_core = total(lambda c, s, e: c == "blood" and (s, e) == (target.start, target.end))
        blood_total = total(lambda c, s, e: c == "blood")
        rows[time] = {"time_hours": time, "antigen_copies_supplied": produced,
                      "target_retained_copies": retained, "extracellular_target_copies": total(lambda c, s, e: c in EXTRACELLULAR),
                      "extracellular_parent_copies": total(lambda c, s, e: c in EXTRACELLULAR and (s, e) == (model.signal_end or 0, len(sequence))),
                      "blood_target_copies": blood_total, "blood_free_core_copies": blood_core,
                      "blood_precursor_copies": blood_total - blood_core,
                      "node_apc_target_copies": total(lambda c, s, e: c in NODE_ANTIGEN),
                      "cumulative_node_apc_entry_copies": arrived,
                      "cumulative_node_apc_entry_by_route": node_entry,
                      "node_apc_entry_fraction_of_dosed_peptide": arrived / initial_copies if model.delivery == "slp" else None,
                      "node_apc_entry_copies_per_input": arrived / initial_copies,
                      "node_apc_entry_fraction_of_supplied_antigen": arrived / produced if produced > 0 else None,
                      "producer_surface_pmhc_copies": total(lambda c, s, e: c == "producer_surface"),
                      "cumulative_producer_loaded_copies": counter("producer_loaded_count"),
                      "target_destroyed_copies": destroyed, "antigen_cleared_copies": cleared,
                      "pmhc_lost_copies": pmhc_lost, "cumulative_loaded_copies": counter("loaded"),
                      "surface_pmhc_copies": sum(surface.values()), "surface_pmhc_by_apc_location": surface,
                      "surface_pmhc_auc_copy_hours": math.fsum(surface_auc.values()),
                      "surface_pmhc_auc_by_apc_location": surface_auc,
                      "extracellular_target_auc_copy_hours": math.fsum(extracellular_auc.values()),
                      "extracellular_target_auc_by_compartment": extracellular_auc,
                      "cumulative_surface_export_copies": counter("surface_exported"),
                      "antigen_balance_error": produced - retained - destroyed - cleared - pmhc_lost,
                      "state_copies": list(map(float, values[:len(states)]))}
    output["curves"] = [rows[time] for time in times]
    return output
