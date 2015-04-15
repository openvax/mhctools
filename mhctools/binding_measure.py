import numpy as np

class BindingMeasure(object):
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
            bigger_is_better,
            min_value=-np.inf,
            min_inclusive=False,
            max_value=np.inf,
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

        bigger_is_better : bool
            When filtering affinities by cutoff, should we pass records which
            are below the cutoff (e.g. for IC50 affinity) or above the cutoff
            (e.g. for stability)?

        min_value : float

        max_value : float
        """
        self.name = name
        self.units = units
        self.bigger_is_better = bigger_is_better
        self.min_value = min_value
        self.min_inclusive = min_inclusive
        self.max_value = max_value
        self.max_inclusive = max_inclusive

    def __str__(self):
        return "%s(units=%s, bigger_is_better=%s)" % (
            self.name,
            self.units,
            self.bigger_is_better
        )

    def fields(self):
        return (
            self.name,
            self.units,
            self.bigger_is_better,
            self.min_value,
            self.min_inclusive,
            self.max_value,
            self.max_inclusive
        )

    def __hash__(self):
        return hash(self.fields())

    def __eq__(self, other):
        return (
            isinstance(other, BindingMeasure) and
            self.fields() == other.fields())

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

    def is_binder(self, value, binding_cutoff):
        """Is the predicted binding value stronger than the given cutoff?

        Parameters
        ----------

        value : float

        cutoff : float
        """
        self.check_binding_value(value)
        if self.bigger_is_better:
            return value <= binding_cutoff
        else:
            return value >= binding_cutoff


ic50_nM = BindingMeasure(
    name="IC50",
    units="nM",
    bigger_is_better=False,
    min_value=0.0)

stability_minutes = BindingMeasure(
    name="Stability (minutes)",
    units="minutes",
    bigger_is_better=True,
    min_value=0.0,
    min_inclusive=False)
