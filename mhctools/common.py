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


def seq_to_str(obj):
    """
    Given a sequence convert it to a comma separated string.
    If, however, the argument is a single object, return its string
    representation.
    """
    if isinstance(obj, (unicode, str)):
        return obj
    elif isinstance(obj, (list, tuple)):
        return ",".join([str(x) for x in obj])
    else:
        return str(obj)

def convert_str(obj):
    """
    Given a string, convert it to an int or float if possible.
    """
    if obj is None:
        return obj
    try:
        try:
            return int(obj)
        except:
            return float(obj)
    except:
        return str(obj)
