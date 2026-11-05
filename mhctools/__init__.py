"""Public API, imported on first access to keep CLI startup lightweight."""

from typing import TYPE_CHECKING

# Preserve static public aliases for documentation and editor discovery.
if TYPE_CHECKING:
    from .binding_prediction import (
        BindingPrediction as BindingPrediction,
    )
    from .binding_prediction_collection import (
        BindingPredictionCollection as BindingPredictionCollection,
    )
    from .logging import (
        set_log_level as set_log_level,
    )
    from .pred import (
        FIELD_BEST_DIRECTIONS as FIELD_BEST_DIRECTIONS,
        Kind as Kind,
        MHC_CLASS_VALUES as MHC_CLASS_VALUES,
        MHC_DEPENDENCE_VALUES as MHC_DEPENDENCE_VALUES,
        MeasurementContext as MeasurementContext,
        PeptidePreds as PeptidePreds,
        PeptideResult as PeptideResult,
        Pred as Pred,
        Prediction as Prediction,
        VALUE_BEST_DIRECTIONS as VALUE_BEST_DIRECTIONS,
        VALUE_UNITS as VALUE_UNITS,
        best_direction as best_direction,
        canonical_kind as canonical_kind,
        intern_measurement_context as intern_measurement_context,
        preds_from_rows as preds_from_rows,
    )
    from .sample import (
        MultiSample as MultiSample,
    )
    from .tcr import (
        TCR as TCR,
    )
    from .allele_formats import (
        caphla_allele_key as caphla_allele_key,
        deepimmuno_allele_name as deepimmuno_allele_name,
        format_allele as format_allele,
        netmhc4_allele_name as netmhc4_allele_name,
        netmhc_allele_name as netmhc_allele_name,
        netmhciipan_allele_name as netmhciipan_allele_name,
        parse_allele as parse_allele,
        prefer_allele_spelling as prefer_allele_spelling,
    )
    from .annotate import (
        AnnotationSpec as AnnotationSpec,
        annotate_table as annotate_table,
        parse_annotation_spec as parse_annotation_spec,
    )
    from .artifacts import (
        ArtifactStatus as ArtifactStatus,
        artifact_status as artifact_status,
        data_path as data_path,
        fetch as fetch,
        list_artifacts as list_artifacts,
        managed_path as managed_path,
    )
    from .optional_backend import (
        ArtifactIdentity as ArtifactIdentity,
        BackendInventory as BackendInventory,
        BackendSpec as BackendSpec,
        backend_inventory as backend_inventory,
        inspect_artifact as inspect_artifact,
        prediction_cache_key as prediction_cache_key,
    )
    from .integrations import (
        IntegrationStatus as IntegrationStatus,
        integration_status as integration_status,
        list_integrations as list_integrations,
    )
    from .peptide_input import (
        PeptideContext as PeptideContext,
        PeptideInput as PeptideInput,
    )
    from .smm import (
        SMM as SMM,
        SMMPMBEC as SMMPMBEC,
    )
    from .iedb import (
        IedbNetMHCcons as IedbNetMHCcons,
        IedbNetMHCpan as IedbNetMHCpan,
        IedbSMM as IedbSMM,
        IedbSMM_PMBEC as IedbSMM_PMBEC,
        IedbNetMHCIIpan as IedbNetMHCIIpan,
    )
    from .mixmhcpred import (
        MixMHCpred as MixMHCpred,
        MixMHCpredAlleleInfo as MixMHCpredAlleleInfo,
        MixMHCpredArtifacts as MixMHCpredArtifacts,
        MixMHCpredResult as MixMHCpredResult,
    )
    from .mixmhc2pred import (
        MixMHC2pred as MixMHC2pred,
    )
    from .prime import (
        PRIME as PRIME,
    )
    from .deeptap import (
        DeepTAP as DeepTAP,
    )
    from .calis import (
        Calis as Calis,
    )
    from .weirdo import (
        WEIRDO as WEIRDO,
    )
    from .eramer import (
        ERAMER as ERAMER,
    )
    from .deepimmuno import (
        DeepImmuno as DeepImmuno,
    )
    from .peptiverse import (
        PeptiVerse as PeptiVerse,
    )
    from .peptiverse_cpp import (
        PeptiVerseCPP as PeptiVerseCPP,
    )
    from .plifepred2 import (
        PlifePred2 as PlifePred2,
    )
    from .cavaco import (
        CavacoHalfLife as CavacoHalfLife,
    )
    from .serum_degradation import (
        DegradationTarget as DegradationTarget,
        TargetFragment as TargetFragment,
        DegradationStep as DegradationStep,
        DegradationPath as DegradationPath,
        target_fragments as target_fragments,
        simulate_target_degradation as simulate_target_degradation,
        degradation_curve as degradation_curve,
        summarize_target_degradation as summarize_target_degradation,
    )
    from .target_cleavage import (
        annotate_target_cleavage as annotate_target_cleavage,
    )
    from .vaccine_windows import (
        assess_vaccine_windows as assess_vaccine_windows,
        vaccine_window_rows as vaccine_window_rows,
        format_vaccine_window_report as format_vaccine_window_report,
    )
    from .vaccine_exposure import (
        summarize_core_exposure as summarize_core_exposure,
    )
    from .vaccine_trajectory import (
        VaccineDose as VaccineDose,
        VaccineRate as VaccineRate,
        VaccineCleavageRates as VaccineCleavageRates,
        VaccineTrajectoryInput as VaccineTrajectoryInput,
        vaccine_route_steps as vaccine_route_steps,
        vaccine_trajectory_evidence as vaccine_trajectory_evidence,
        simulate_vaccine_trajectory as simulate_vaccine_trajectory,
    )
    from .vaccine_trajectory_report import (
        compare_vaccine_trajectories as compare_vaccine_trajectories,
        format_vaccine_comparison_report as format_vaccine_comparison_report,
        summarize_vaccine_trajectory as summarize_vaccine_trajectory,
        format_vaccine_trajectory_report as format_vaccine_trajectory_report,
    )
    from .serum_contributions import (
        EnzymeCutRate as EnzymeCutRate,
        serum_contribution_evidence as serum_contribution_evidence,
        serum_reference_kinetics as serum_reference_kinetics,
        empirical_terminal_cut_rates as empirical_terminal_cut_rates,
        simulate_enzyme_degradation as simulate_enzyme_degradation,
        enzyme_removal_effects as enzyme_removal_effects,
    )
    from .serum_calibration import (
        serum_calibration_evidence as serum_calibration_evidence,
        serum_assay_parent_reference as serum_assay_parent_reference,
    )
    from .cleavage import (
        CleavageModel as CleavageModel,
        CleavageSite as CleavageSite,
        CleavageResult as CleavageResult,
    )
    from .cleavage_batch import (
        cleavage_overlays as cleavage_overlays,
        load_cleavage_batch as load_cleavage_batch,
        normalize_cleavage_input as normalize_cleavage_input,
        predict_cleavage_batch as predict_cleavage_batch,
        write_cleavage_batch as write_cleavage_batch,
    )
    from .dpp4 import (
        DPP4qPISA as DPP4qPISA,
    )
    from .enzyme_assay import (
        AssayConditions as AssayConditions,
        estimate_assay_parent_fraction as estimate_assay_parent_fraction,
    )
    from .enzyme_predictor import (
        DPP4 as DPP4,
        FAP as FAP,
        EnzymePredictor as EnzymePredictor,
        EnzymePredictions as EnzymePredictions,
        enzyme_predictors as enzyme_predictors,
        enzymes as enzymes,
    )
    from .peptidases import (
        cleavage_models as cleavage_models,
        get_cleavage_model as get_cleavage_model,
        motif_rule as motif_rule,
        predict_cleavage as predict_cleavage,
    )
    from .eramer_cleavage import (
        ERAMERCleavage as ERAMERCleavage,
    )
    from .itcell_cleavage import (
        ITCellCleavage as ITCellCleavage,
    )
    from .phagescout import (
        PhageScout as PhageScout,
    )
    from .tlimmuno2 import (
        TLimmuno2 as TLimmuno2,
    )
    from .processing_predictor import (
        ProcessingPredictor as ProcessingPredictor,
        SCORING_MODES as SCORING_MODES,
        resolve_scoring as resolve_scoring,
        score_cterm as score_cterm,
        score_nterm_cterm as score_nterm_cterm,
        score_cterm_anti_max_internal as score_cterm_anti_max_internal,
        score_cterm_anti_mean_internal as score_cterm_anti_mean_internal,
        score_nterm_cterm_anti_max_internal as score_nterm_cterm_anti_max_internal,
        score_nterm_cterm_anti_mean_internal as score_nterm_cterm_anti_mean_internal,
    )
    from .proteasome_predictor import (
        ProteasomePredictor as ProteasomePredictor,
    )
    from .netchop import (
        NetChop as NetChop,
    )
    from .netcleave import (
        NetCleave as NetCleave,
        NetCleave_I as NetCleave_I,
        NetCleave_II as NetCleave_II,
    )
    from .pepsickle import (
        Pepsickle as Pepsickle,
        PepsickleCleavage as PepsickleCleavage,
    )
    from .vaccine_report import (
        CleavageTrack as CleavageTrack,
        EpitopeWindow as EpitopeWindow,
        MHCWindow as MHCWindow,
        RoutePolicy as RoutePolicy,
        VaccineConstruct as VaccineConstruct,
        VaccineReportInput as VaccineReportInput,
        generate_vaccine_report as generate_vaccine_report,
        compact_target_summary as compact_target_summary,
        placement_assessments as placement_assessments,
        processing_route_policy as processing_route_policy,
        select_mhc_windows as select_mhc_windows,
    )
    from .netmhc import (
        NetMHC as NetMHC,
    )
    from .netmhc3 import (
        NetMHC3 as NetMHC3,
    )
    from .netmhc4 import (
        NetMHC4 as NetMHC4,
    )
    from .netmhc_cons import (
        NetMHCcons as NetMHCcons,
    )
    from .netmhc_pan import (
        NetMHCpan as NetMHCpan,
    )
    from .netmhc_pan28 import (
        NetMHCpan28 as NetMHCpan28,
    )
    from .netmhc_pan3 import (
        NetMHCpan3 as NetMHCpan3,
    )
    from .netmhc_pan4 import (
        NetMHCpan4 as NetMHCpan4,
        NetMHCpan4_BA as NetMHCpan4_BA,
        NetMHCpan4_EL as NetMHCpan4_EL,
    )
    from .netmhc_pan41 import (
        NetMHCpan41 as NetMHCpan41,
        NetMHCpan41_BA as NetMHCpan41_BA,
        NetMHCpan41_EL as NetMHCpan41_EL,
    )
    from .netmhc_pan42 import (
        NetMHCpan42 as NetMHCpan42,
        NetMHCpan42_BA as NetMHCpan42_BA,
        NetMHCpan42_EL as NetMHCpan42_EL,
    )
    from .netmhcii_pan import (
        NetMHCIIpan as NetMHCIIpan,
        NetMHCIIpan3 as NetMHCIIpan3,
        NetMHCIIpan4 as NetMHCIIpan4,
        NetMHCIIpan4_BA as NetMHCIIpan4_BA,
        NetMHCIIpan4_EL as NetMHCIIpan4_EL,
        NetMHCIIpan43 as NetMHCIIpan43,
        NetMHCIIpan43_BA as NetMHCIIpan43_BA,
        NetMHCIIpan43_EL as NetMHCIIpan43_EL,
    )
    from .random_predictor import (
        RandomBindingPredictor as RandomBindingPredictor,
    )
    from .netmhcstabpan import (
        NetMHCstabpan as NetMHCstabpan,
    )
    from .unsupported_allele import (
        UnsupportedAllele as UnsupportedAllele,
    )
    from .cleavenet import (
        CleaveNet as CleaveNet,
        CleaveNetResult as CleaveNetResult,
        CleaveNetScore as CleaveNetScore,
    )
    from .bigmhc import (
        BigMHC as BigMHC,
        BigMHC_EL as BigMHC_EL,
        BigMHC_IM as BigMHC_IM,
    )
    from .caphla import (
        CapHLA as CapHLA,
        CapHLA_BA as CapHLA_BA,
        CapHLA_EL as CapHLA_EL,
    )
    from .mhcflurry import (
        MHCflurry as MHCflurry,
        MHCflurry_Affinity as MHCflurry_Affinity,
        mhcflurry_composite_version as mhcflurry_composite_version,
    )
    from .mixtcrpred import (
        MixTCRpred as MixTCRpred,
        MixTCRpredModel as MixTCRpredModel,
    )
    from .nettcr import (
        NetTCR as NetTCR,
    )
    from .tulip import (
        Tulip as Tulip,
    )


_LAZY_IMPORTS = {
    "BindingPrediction": (".binding_prediction", "BindingPrediction"),
    "BindingPredictionCollection": (".binding_prediction_collection", "BindingPredictionCollection"),
    "set_log_level": (".logging", "set_log_level"),
    "FIELD_BEST_DIRECTIONS": (".pred", "FIELD_BEST_DIRECTIONS"),
    "Kind": (".pred", "Kind"),
    "MHC_CLASS_VALUES": (".pred", "MHC_CLASS_VALUES"),
    "MHC_DEPENDENCE_VALUES": (".pred", "MHC_DEPENDENCE_VALUES"),
    "MeasurementContext": (".pred", "MeasurementContext"),
    "PeptidePreds": (".pred", "PeptidePreds"),
    "PeptideResult": (".pred", "PeptideResult"),
    "Pred": (".pred", "Pred"),
    "Prediction": (".pred", "Prediction"),
    "VALUE_BEST_DIRECTIONS": (".pred", "VALUE_BEST_DIRECTIONS"),
    "VALUE_UNITS": (".pred", "VALUE_UNITS"),
    "best_direction": (".pred", "best_direction"),
    "canonical_kind": (".pred", "canonical_kind"),
    "intern_measurement_context": (".pred", "intern_measurement_context"),
    "preds_from_rows": (".pred", "preds_from_rows"),
    "MultiSample": (".sample", "MultiSample"),
    "TCR": (".tcr", "TCR"),
    "caphla_allele_key": (".allele_formats", "caphla_allele_key"),
    "deepimmuno_allele_name": (".allele_formats", "deepimmuno_allele_name"),
    "format_allele": (".allele_formats", "format_allele"),
    "netmhc4_allele_name": (".allele_formats", "netmhc4_allele_name"),
    "netmhc_allele_name": (".allele_formats", "netmhc_allele_name"),
    "netmhciipan_allele_name": (".allele_formats", "netmhciipan_allele_name"),
    "parse_allele": (".allele_formats", "parse_allele"),
    "prefer_allele_spelling": (".allele_formats", "prefer_allele_spelling"),
    "AnnotationSpec": (".annotate", "AnnotationSpec"),
    "annotate_table": (".annotate", "annotate_table"),
    "parse_annotation_spec": (".annotate", "parse_annotation_spec"),
    "ArtifactStatus": (".artifacts", "ArtifactStatus"),
    "artifact_status": (".artifacts", "artifact_status"),
    "data_path": (".artifacts", "data_path"),
    "fetch": (".artifacts", "fetch"),
    "list_artifacts": (".artifacts", "list_artifacts"),
    "managed_path": (".artifacts", "managed_path"),
    "ArtifactIdentity": (".optional_backend", "ArtifactIdentity"),
    "BackendInventory": (".optional_backend", "BackendInventory"),
    "BackendSpec": (".optional_backend", "BackendSpec"),
    "backend_inventory": (".optional_backend", "backend_inventory"),
    "inspect_artifact": (".optional_backend", "inspect_artifact"),
    "prediction_cache_key": (".optional_backend", "prediction_cache_key"),
    "IntegrationStatus": (".integrations", "IntegrationStatus"),
    "integration_status": (".integrations", "integration_status"),
    "list_integrations": (".integrations", "list_integrations"),
    "PeptideContext": (".peptide_input", "PeptideContext"),
    "PeptideInput": (".peptide_input", "PeptideInput"),
    "SMM": (".smm", "SMM"),
    "SMMPMBEC": (".smm", "SMMPMBEC"),
    "IedbNetMHCcons": (".iedb", "IedbNetMHCcons"),
    "IedbNetMHCpan": (".iedb", "IedbNetMHCpan"),
    "IedbSMM": (".iedb", "IedbSMM"),
    "IedbSMM_PMBEC": (".iedb", "IedbSMM_PMBEC"),
    "IedbNetMHCIIpan": (".iedb", "IedbNetMHCIIpan"),
    "MixMHCpred": (".mixmhcpred", "MixMHCpred"),
    "MixMHCpredAlleleInfo": (".mixmhcpred", "MixMHCpredAlleleInfo"),
    "MixMHCpredArtifacts": (".mixmhcpred", "MixMHCpredArtifacts"),
    "MixMHCpredResult": (".mixmhcpred", "MixMHCpredResult"),
    "MixMHC2pred": (".mixmhc2pred", "MixMHC2pred"),
    "PRIME": (".prime", "PRIME"),
    "DeepTAP": (".deeptap", "DeepTAP"),
    "Calis": (".calis", "Calis"),
    "WEIRDO": (".weirdo", "WEIRDO"),
    "ERAMER": (".eramer", "ERAMER"),
    "DeepImmuno": (".deepimmuno", "DeepImmuno"),
    "PeptiVerse": (".peptiverse", "PeptiVerse"),
    "PeptiVerseCPP": (".peptiverse_cpp", "PeptiVerseCPP"),
    "PlifePred2": (".plifepred2", "PlifePred2"),
    "CavacoHalfLife": (".cavaco", "CavacoHalfLife"),
    "DegradationTarget": (".serum_degradation", "DegradationTarget"),
    "TargetFragment": (".serum_degradation", "TargetFragment"),
    "DegradationStep": (".serum_degradation", "DegradationStep"),
    "DegradationPath": (".serum_degradation", "DegradationPath"),
    "target_fragments": (".serum_degradation", "target_fragments"),
    "simulate_target_degradation": (".serum_degradation", "simulate_target_degradation"),
    "degradation_curve": (".serum_degradation", "degradation_curve"),
    "summarize_target_degradation": (".serum_degradation", "summarize_target_degradation"),
    "annotate_target_cleavage": (".target_cleavage", "annotate_target_cleavage"),
    "assess_vaccine_windows": (".vaccine_windows", "assess_vaccine_windows"),
    "vaccine_window_rows": (".vaccine_windows", "vaccine_window_rows"),
    "format_vaccine_window_report": (".vaccine_windows", "format_vaccine_window_report"),
    "summarize_core_exposure": (".vaccine_exposure", "summarize_core_exposure"),
    "VaccineDose": (".vaccine_trajectory", "VaccineDose"),
    "VaccineRate": (".vaccine_trajectory", "VaccineRate"),
    "VaccineCleavageRates": (".vaccine_trajectory", "VaccineCleavageRates"),
    "VaccineTrajectoryInput": (".vaccine_trajectory", "VaccineTrajectoryInput"),
    "vaccine_route_steps": (".vaccine_trajectory", "vaccine_route_steps"),
    "vaccine_trajectory_evidence": (".vaccine_trajectory", "vaccine_trajectory_evidence"),
    "simulate_vaccine_trajectory": (".vaccine_trajectory", "simulate_vaccine_trajectory"),
    "compare_vaccine_trajectories": (".vaccine_trajectory_report", "compare_vaccine_trajectories"),
    "format_vaccine_comparison_report": (".vaccine_trajectory_report", "format_vaccine_comparison_report"),
    "summarize_vaccine_trajectory": (".vaccine_trajectory_report", "summarize_vaccine_trajectory"),
    "format_vaccine_trajectory_report": (".vaccine_trajectory_report", "format_vaccine_trajectory_report"),
    "EnzymeCutRate": (".serum_contributions", "EnzymeCutRate"),
    "serum_contribution_evidence": (".serum_contributions", "serum_contribution_evidence"),
    "serum_reference_kinetics": (".serum_contributions", "serum_reference_kinetics"),
    "empirical_terminal_cut_rates": (".serum_contributions", "empirical_terminal_cut_rates"),
    "simulate_enzyme_degradation": (".serum_contributions", "simulate_enzyme_degradation"),
    "enzyme_removal_effects": (".serum_contributions", "enzyme_removal_effects"),
    "serum_calibration_evidence": (".serum_calibration", "serum_calibration_evidence"),
    "serum_assay_parent_reference": (".serum_calibration", "serum_assay_parent_reference"),
    "CleavageModel": (".cleavage", "CleavageModel"),
    "CleavageSite": (".cleavage", "CleavageSite"),
    "CleavageResult": (".cleavage", "CleavageResult"),
    "cleavage_overlays": (".cleavage_batch", "cleavage_overlays"),
    "load_cleavage_batch": (".cleavage_batch", "load_cleavage_batch"),
    "normalize_cleavage_input": (".cleavage_batch", "normalize_cleavage_input"),
    "predict_cleavage_batch": (".cleavage_batch", "predict_cleavage_batch"),
    "write_cleavage_batch": (".cleavage_batch", "write_cleavage_batch"),
    "DPP4qPISA": (".dpp4", "DPP4qPISA"),
    "AssayConditions": (".enzyme_assay", "AssayConditions"),
    "estimate_assay_parent_fraction": (".enzyme_assay", "estimate_assay_parent_fraction"),
    "DPP4": (".enzyme_predictor", "DPP4"),
    "FAP": (".enzyme_predictor", "FAP"),
    "EnzymePredictor": (".enzyme_predictor", "EnzymePredictor"),
    "EnzymePredictions": (".enzyme_predictor", "EnzymePredictions"),
    "enzyme_predictors": (".enzyme_predictor", "enzyme_predictors"),
    "enzymes": (".enzyme_predictor", "enzymes"),
    "cleavage_models": (".peptidases", "cleavage_models"),
    "get_cleavage_model": (".peptidases", "get_cleavage_model"),
    "motif_rule": (".peptidases", "motif_rule"),
    "predict_cleavage": (".peptidases", "predict_cleavage"),
    "ERAMERCleavage": (".eramer_cleavage", "ERAMERCleavage"),
    "ITCellCleavage": (".itcell_cleavage", "ITCellCleavage"),
    "PhageScout": (".phagescout", "PhageScout"),
    "TLimmuno2": (".tlimmuno2", "TLimmuno2"),
    "ProcessingPredictor": (".processing_predictor", "ProcessingPredictor"),
    "SCORING_MODES": (".processing_predictor", "SCORING_MODES"),
    "resolve_scoring": (".processing_predictor", "resolve_scoring"),
    "score_cterm": (".processing_predictor", "score_cterm"),
    "score_nterm_cterm": (".processing_predictor", "score_nterm_cterm"),
    "score_cterm_anti_max_internal": (".processing_predictor", "score_cterm_anti_max_internal"),
    "score_cterm_anti_mean_internal": (".processing_predictor", "score_cterm_anti_mean_internal"),
    "score_nterm_cterm_anti_max_internal": (".processing_predictor", "score_nterm_cterm_anti_max_internal"),
    "score_nterm_cterm_anti_mean_internal": (".processing_predictor", "score_nterm_cterm_anti_mean_internal"),
    "ProteasomePredictor": (".proteasome_predictor", "ProteasomePredictor"),
    "NetChop": (".netchop", "NetChop"),
    "NetCleave": (".netcleave", "NetCleave"),
    "NetCleave_I": (".netcleave", "NetCleave_I"),
    "NetCleave_II": (".netcleave", "NetCleave_II"),
    "Pepsickle": (".pepsickle", "Pepsickle"),
    "PepsickleCleavage": (".pepsickle", "PepsickleCleavage"),
    "CleavageTrack": (".vaccine_report", "CleavageTrack"),
    "EpitopeWindow": (".vaccine_report", "EpitopeWindow"),
    "MHCWindow": (".vaccine_report", "MHCWindow"),
    "RoutePolicy": (".vaccine_report", "RoutePolicy"),
    "VaccineConstruct": (".vaccine_report", "VaccineConstruct"),
    "VaccineReportInput": (".vaccine_report", "VaccineReportInput"),
    "generate_vaccine_report": (".vaccine_report", "generate_vaccine_report"),
    "compact_target_summary": (".vaccine_report", "compact_target_summary"),
    "placement_assessments": (".vaccine_report", "placement_assessments"),
    "processing_route_policy": (".vaccine_report", "processing_route_policy"),
    "select_mhc_windows": (".vaccine_report", "select_mhc_windows"),
    "NetMHC": (".netmhc", "NetMHC"),
    "NetMHC3": (".netmhc3", "NetMHC3"),
    "NetMHC4": (".netmhc4", "NetMHC4"),
    "NetMHCcons": (".netmhc_cons", "NetMHCcons"),
    "NetMHCpan": (".netmhc_pan", "NetMHCpan"),
    "NetMHCpan28": (".netmhc_pan28", "NetMHCpan28"),
    "NetMHCpan3": (".netmhc_pan3", "NetMHCpan3"),
    "NetMHCpan4": (".netmhc_pan4", "NetMHCpan4"),
    "NetMHCpan4_BA": (".netmhc_pan4", "NetMHCpan4_BA"),
    "NetMHCpan4_EL": (".netmhc_pan4", "NetMHCpan4_EL"),
    "NetMHCpan41": (".netmhc_pan41", "NetMHCpan41"),
    "NetMHCpan41_BA": (".netmhc_pan41", "NetMHCpan41_BA"),
    "NetMHCpan41_EL": (".netmhc_pan41", "NetMHCpan41_EL"),
    "NetMHCpan42": (".netmhc_pan42", "NetMHCpan42"),
    "NetMHCpan42_BA": (".netmhc_pan42", "NetMHCpan42_BA"),
    "NetMHCpan42_EL": (".netmhc_pan42", "NetMHCpan42_EL"),
    "NetMHCIIpan": (".netmhcii_pan", "NetMHCIIpan"),
    "NetMHCIIpan3": (".netmhcii_pan", "NetMHCIIpan3"),
    "NetMHCIIpan4": (".netmhcii_pan", "NetMHCIIpan4"),
    "NetMHCIIpan4_BA": (".netmhcii_pan", "NetMHCIIpan4_BA"),
    "NetMHCIIpan4_EL": (".netmhcii_pan", "NetMHCIIpan4_EL"),
    "NetMHCIIpan43": (".netmhcii_pan", "NetMHCIIpan43"),
    "NetMHCIIpan43_BA": (".netmhcii_pan", "NetMHCIIpan43_BA"),
    "NetMHCIIpan43_EL": (".netmhcii_pan", "NetMHCIIpan43_EL"),
    "RandomBindingPredictor": (".random_predictor", "RandomBindingPredictor"),
    "NetMHCstabpan": (".netmhcstabpan", "NetMHCstabpan"),
    "UnsupportedAllele": (".unsupported_allele", "UnsupportedAllele"),
    "CleaveNet": (".cleavenet", "CleaveNet"),
    "CleaveNetResult": (".cleavenet", "CleaveNetResult"),
    "CleaveNetScore": (".cleavenet", "CleaveNetScore"),
    "BigMHC": (".bigmhc", "BigMHC"),
    "BigMHC_EL": (".bigmhc", "BigMHC_EL"),
    "BigMHC_IM": (".bigmhc", "BigMHC_IM"),
    "CapHLA": (".caphla", "CapHLA"),
    "CapHLA_BA": (".caphla", "CapHLA_BA"),
    "CapHLA_EL": (".caphla", "CapHLA_EL"),
    "MHCflurry": (".mhcflurry", "MHCflurry"),
    "MHCflurry_Affinity": (".mhcflurry", "MHCflurry_Affinity"),
    "mhcflurry_composite_version": (".mhcflurry", "mhcflurry_composite_version"),
    "MixTCRpred": (".mixtcrpred", "MixTCRpred"),
    "MixTCRpredModel": (".mixtcrpred", "MixTCRpredModel"),
    "NetTCR": (".nettcr", "NetTCR"),
    "Tulip": (".tulip", "Tulip"),
}


def __getattr__(name):
    """PEP 562: load public objects on first access."""
    if name in _LAZY_IMPORTS:
        from importlib import import_module
        module_name, attr = _LAZY_IMPORTS[name]
        module = import_module(module_name, package=__name__)
        value = getattr(module, attr)
        globals()[name] = value
        return value
    raise AttributeError(
        "module %r has no attribute %r" % (__name__, name))

def __dir__():
    """Include public exports in interactive discovery before their first use."""
    return sorted(set(globals()) | set(_LAZY_IMPORTS))


__version__ = "4.8.3"

__all__ = [
    "CleaveNet",
    "CleaveNetResult",
    "CleaveNetScore",
    "mhcflurry_composite_version",
    "Prediction",
    "MeasurementContext",
    "PeptideContext",
    "PeptideInput",
    "Pred",  # backward compat alias
    "PeptideResult",
    "PeptidePreds",  # backward compat alias
    "Kind",
    "MHC_CLASS_VALUES",
    "MHC_DEPENDENCE_VALUES",
    "FIELD_BEST_DIRECTIONS",
    "VALUE_BEST_DIRECTIONS",
    "VALUE_UNITS",
    "best_direction",
    "canonical_kind",
    "intern_measurement_context",
    "preds_from_rows",
    "MultiSample",
    "TCR",
    "AnnotationSpec",
    "annotate_table",
    "caphla_allele_key",
    "deepimmuno_allele_name",
    "format_allele",
    "netmhc4_allele_name",
    "netmhc_allele_name",
    "netmhciipan_allele_name",
    "parse_allele",
    "prefer_allele_spelling",
    "parse_annotation_spec",
    "ArtifactStatus",
    "artifact_status",
    "data_path",
    "fetch",
    "list_artifacts",
    "managed_path",
    "ArtifactIdentity",
    "BackendInventory",
    "BackendSpec",
    "backend_inventory",
    "inspect_artifact",
    "prediction_cache_key",
    "IntegrationStatus",
    "integration_status",
    "list_integrations",
    "BindingPrediction",
    "BindingPredictionCollection",
    "set_log_level",
    "IedbNetMHCcons",
    "IedbNetMHCpan",
    "SMM",
    "SMMPMBEC",
    "IedbSMM",
    "IedbSMM_PMBEC",
    "IedbNetMHCIIpan",
    "MixMHCpred",
    "MixMHCpredAlleleInfo",
    "MixMHCpredArtifacts",
    "MixMHCpredResult",
    "MixMHC2pred",
    "PRIME",
    "DeepTAP",
    "Calis",
    "WEIRDO",
    "ERAMER",
    "DeepImmuno",
    "PeptiVerse",
    "PeptiVerseCPP",
    "PlifePred2",
    "CavacoHalfLife",
    "DegradationTarget",
    "TargetFragment",
    "DegradationStep",
    "DegradationPath",
    "target_fragments",
    "simulate_target_degradation",
    "degradation_curve",
    "summarize_target_degradation",
    "annotate_target_cleavage",
    "VaccineRate",
    "VaccineDose",
    "VaccineCleavageRates",
    "VaccineTrajectoryInput",
    "vaccine_route_steps",
    "vaccine_trajectory_evidence",
    "simulate_vaccine_trajectory",
    "summarize_core_exposure",
    "assess_vaccine_windows",
    "vaccine_window_rows",
    "format_vaccine_window_report",
    "summarize_vaccine_trajectory",
    "format_vaccine_trajectory_report",
    "compare_vaccine_trajectories",
    "format_vaccine_comparison_report",
    "EnzymeCutRate",
    "serum_contribution_evidence",
    "serum_calibration_evidence",
    "serum_assay_parent_reference",
    "compact_target_summary",
    "serum_reference_kinetics",
    "empirical_terminal_cut_rates",
    "simulate_enzyme_degradation",
    "enzyme_removal_effects",
    "CleavageModel",
    "CleavageSite",
    "CleavageResult",
    "DPP4qPISA",
    "EnzymePredictor",
    "EnzymePredictions",
    "DPP4",
    "FAP",
    "enzyme_predictors",
    "enzymes",
    "AssayConditions",
    "estimate_assay_parent_fraction",
    "cleavage_models",
    "get_cleavage_model",
    "motif_rule",
    "predict_cleavage",
    "predict_cleavage_batch",
    "normalize_cleavage_input",
    "cleavage_overlays",
    "load_cleavage_batch",
    "write_cleavage_batch",
    "ERAMERCleavage",
    "ITCellCleavage",
    "PhageScout",
    "TLimmuno2",
    "MHCflurry",
    "MHCflurry_Affinity",
    "ProcessingPredictor",
    "ProteasomePredictor",
    "SCORING_MODES",
    "resolve_scoring",
    "score_cterm",
    "score_nterm_cterm",
    "score_cterm_anti_max_internal",
    "score_cterm_anti_mean_internal",
    "score_nterm_cterm_anti_max_internal",
    "score_nterm_cterm_anti_mean_internal",
    "NetChop",
    "NetCleave",
    "NetCleave_I",
    "NetCleave_II",
    "Pepsickle",
    "NetMHC",
    "NetMHC3",
    "NetMHC4",
    "NetMHCcons",
    "NetMHCpan",
    "NetMHCpan28",
    "NetMHCpan3",
    "NetMHCpan4",
    "NetMHCpan41",
    "NetMHCpan4_BA",
    "NetMHCpan4_EL",
    "NetMHCpan41_BA",
    "NetMHCpan41_EL",
    "NetMHCpan42",
    "NetMHCpan42_BA",
    "NetMHCpan42_EL",
    "NetMHCIIpan",
    "NetMHCIIpan3",
    "NetMHCIIpan4",
    "NetMHCIIpan4_BA",
    "NetMHCIIpan4_EL",
    "NetMHCIIpan43",
    "NetMHCIIpan43_BA",
    "NetMHCIIpan43_EL",
    "NetMHCstabpan",
    "BigMHC",
    "BigMHC_EL",
    "BigMHC_IM",
    "CapHLA",
    "CapHLA_BA",
    "CapHLA_EL",
    "MixTCRpred",
    "MixTCRpredModel",
    "NetTCR",
    "Tulip",
    "RandomBindingPredictor",
    "UnsupportedAllele",
]
