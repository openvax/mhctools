import numpy as np

class PeptideBindingMeasure(object):
    """
    We could use affinity/stability measures informally but it's cleaner
    to keep track of the names, direction, and lower/upper bounds of
    each assay measure in a single place. For example, percentiles are
    *strong* if lower and should never exceed 100.0. On the other hand,
    the stability of peptide-MHC association is considered stronger for
    higher values and has no upper bound.
    """

    def __init__(
            self,
            name,
            units,
            record_field_name,
            cutoff_is_upper_bound,
            min_value = -np.inf,
            min_inclusive=False,
            max_value = np.inf,
            max_inclusive=False):
        """
        Parameters
        ----------

        name : str

        units : str

        record_field_name : str
            The output of MHC binding/stability prediction is a collection
            of per-allele records, which field should we look at in those
            records?

        cutoff_is_upper_bound : bool
            When filtering affinities by cutoff, should we pass records which
            are below the cutoff (e.g. for IC50 affinity) or above the cutoff
            (e.g. for stability)?

        min_value : float

        max_value : float
        """
        self.name = name
        self.units = units
        self.record_field_name  = record_field_name
        self.cutoff_is_upper_bound = cutoff_is_upper_bound
        self.min_value = min_value
        self.min_inclusive = min_inclusive
        self.max_value = max_value
        self.max_inclusive = max_inclusive

    def __str__(self):
        return "%s(units=%s, field=%s, cutoff_is_upper_bound=%s)" % (
            self.name, self.units, self.prediction_record_field,
            self.cutoff_is_upper_bound
        )

    def __repr__(self):
        return str(self)

    def check_binding_value(self, value):
        """
        Ensure that predicted binding value falls within the specified
        (min,max) range or raise an assertion error.
        """
        assert isinstance(value, (int, float)), \
            "Expected float for binding value, got %s : %s" % (
                value, type(value))

        if self.min_inclusive:
            assert value >= self.min_value, \
                "Given value (%s) too low (min_value=%s)" % (
                    value, self.min_value)
        else:
            assert value > self.min_value, \
                "Given value (%s) too low (min_value=%s)" % (
                    value, self.min_value)
        if self.max_inclusive:
            assert value <= self.max_value, \
                "Given value (%s) too high (max_value=%s)" % (
                    value, self.max_value)
        else:
            assert value < self.max_value, \
                "Given value (%s) too high (max_value=%s)" % (
                    value, self.max_value)

    def extract_value(self, record):
        field_name = self.record_field_name
        if hasattr(record, field_name):
           return getattr(record_field_name)
        else:
            assert isinstance(record, dict), \
                "Invalid prediction record type %s : %s" % (
                    record, type(record))
            assert field_name in record
            return record[field_name]


    def value_is_binder(self, value, cutoff):
        """
        Is the predicted binding value stronger than the given cutoff?

        Parameters
        ----------

        value : float

        cutoff : float
        """
        self.check_binding_value(value)
        if self.cutoff_is_upper_bound:
            return value <= cutoff
        else:
            return value >= cutoff

    def record_is_binder(self, record, cutoff):
        value = self.extract_value(record)
        return self.value_is_binder(value, cutoff)


# used by other modules to name prediction fields
IC50_FIELD_NAME = "MHC_IC50"
ic50_binding_measure = PeptideBindingMeasure(
    name="IC50",
    units="nM",
    record_field_name=IC50_FIELD_NAME,
    cutoff_is_upper_bound=True,
    min_value = 0.0)



PERCENTILE_RANK_FIELD_NAME = "MHC_Percentile_Rank"
percentile_binding_measure = PeptideBindingMeasure(
    name="PercentileRank",
    units="%",
    record_field_name=PERCENTILE_RANK_FIELD_NAME,
    cutoff_is_upper_bound=True,
    min_value = 0.0,
    min_inclusive=True,
    max_value=100.0,
    max_inclusive=True)

STABILITY_FIELD_NAME = "MHC_Stability"
stability_binding_measure = PeptideBindingMeasure(
    name="Stability",
    units="minutes",
    record_field_name=STABILITY_FIELD_NAME,
    cutoff_is_upper_bound=False,
    min_value=0.0,
    min_inclusive=False)
