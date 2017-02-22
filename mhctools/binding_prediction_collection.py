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

from __future__ import print_function, division, absolute_import

import pandas as pd
from sercol import Collection

from .binding_prediction import BindingPrediction

class BindingPredictionCollection(Collection):
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

    def groupby_allele(self):
        return self.groupby(key_fn=lambda x: x.allele)

    def groupby_peptide(self):
        return self.groupby(key_fn=lambda x: x.peptide)

    def groupby_allele_and_peptide(self):
        return self.groupby(key_fn=lambda x: (x.allele, x.peptide))

    def to_dataframe(self):
        """
        Converts collection of BindingPrediction objects to DataFrame
        """
        return pd.DataFrame(
            self.elements,
            columns=BindingPrediction._fields)

    def update_fields(self, **kwargs):
        """
        Changes fields corresponding to names of keyword arguments in
        each BindingPrediction.
        """
        field_dicts = [x._asdict() for x in self]
        for field_name, values in kwargs.items():
            for i, value in enumerate(values):
                field_dicts[i][field_name] = value
        new_binding_predictions = [
            BindingPrediction(**d) for d in field_dicts
        ]
        return self.clone_with_new_elements(new_binding_predictions)
