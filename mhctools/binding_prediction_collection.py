# Copyright (c) 2017. Mount Sinai School of Medicine
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
    def to_dataframe(
            self,
            columns=BindingPrediction.fields + ("length",)):
        """
        Converts collection of BindingPrediction objects to DataFrame
        """
        return pd.DataFrame.from_records(
            [tuple([getattr(x, name) for name in columns]) for x in self],
            columns=columns)
