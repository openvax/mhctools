# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

import logging
import logging.config
from importlib.resources import as_file, files


def get_logger(name):
    config_resource = files("mhctools").joinpath("logging.conf")
    with as_file(config_resource) as config_path:
        logging.config.fileConfig(str(config_path))
    return logging.getLogger(name)
