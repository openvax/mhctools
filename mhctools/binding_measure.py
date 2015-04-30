# Copyright (c) 2014. Mount Sinai School of Medicine
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

from __future__ import print_function, division, absolute_import

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
            default_cutoff,
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

        default_cutoff : float
            Default value separating strong from weak binders

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
        self.default_cutoff = default_cutoff
        self.min_value = min_value
        self.min_inclusive = min_inclusive
        self.max_value = max_value
        self.max_inclusive = max_inclusive

    def __str__(self):
        return "%s_%s" % (self.name, self.units)

    def fields(self):
        return (
            self.name,
            self.units,
            self.bigger_is_better,
            self.default_cutoff,
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

    def is_binder(self, value, binding_cutoff=None):
        """Is the predicted binding value stronger than (or equal to) the given
        cutoff?

        Parameters
        ----------

        value : float

        cutoff : float
        """
        if binding_cutoff is None:
            binding_cutoff = self.default_cutoff
        else:
            self.check_binding_value(value)
        if self.bigger_is_better:
            return value >= binding_cutoff
        else:
            return value <= binding_cutoff

ic50_nM = BindingMeasure(
    name="IC50",
    units="nM",
    bigger_is_better=False,
    default_cutoff=500.0,
    min_value=0.0)

stability_minutes = BindingMeasure(
    name="Stability",
    units="minutes",
    bigger_is_better=True,
    # 6 hours
    default_cutoff=360.0,
    min_value=0.0,
    min_inclusive=False)
