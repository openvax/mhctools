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
from collections import defaultdict

import pandas as pd
from sercol import Collection

from .binding_prediction import BindingPrediction

class EpitopeCollection(Collection):
    """
    Collection of BindingPrediction objects
    """
    def __init__(
            self,
            binding_predictions,
            path=None,
            distinct=True,
            sort_key=lambda x: x.value):
        Collection.__init__(
            self,
            elements=binding_predictions,
            sources=[path] if path else [],
            distinct=distinct,
            sort_key=sort_key)

    def strong_binders(self, threshold=None):
        """
        No default threshold since we're not sure if we're always going to
        be predicting IC50 affinity (as opposed to stability or some other
        criterion)
        """
        if len(self) == 0:
            return self

        def filter_fn(x):
            if threshold is None:
                return x.measure.is_binder(x.value)
            else:
                return x.measure.is_binder(x.value, threshold)
        return self.filter(filter_fn)

    def strong_binders_by_rank(self, max_rank=2.0):
        return self.filter(lambda x: x.percentile_rank <= max_rank)

    def groupby(self, key_fn):
        groups = defaultdict(list)
        for binding_prediction in self.elements:
            key = key_fn(binding_prediction)
            groups[key].append(binding_prediction)
        # want to create an EpitopeCollection for each group
        # but need to write the construction in terms of
        # self.__class__ so that this works with derived classes
        return {
            key: self.__class__(binding_predictions)
            for (key, binding_predictions)
            in groups.items()
        }

    def groupby_allele(self):
        return self.groupby(key_fn=lambda x: x.allele)

    def groupby_peptide(self):
        return self.groupby(key_fn=lambda x: x.peptide)

    def groupby_allele_and_peptide(self):
        return self.groupby(key_fn=lambda x: (x.allele, x.peptide))

    def dataframe(self):
        return pd.DataFrame(
            self.elements,
            columns=BindingPrediction._fields)
