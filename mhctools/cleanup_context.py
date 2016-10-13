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
import logging
import os
import shutil


logger = logging.getLogger(__name__)


class CleanupFiles(object):
    """
    Context manager that deletes a set of files at the end of a block or
    if an exception gets raised
    """
    def __init__(
            self,
            files=[],
            filenames=[],
            dictionaries=[],
            directories=[]):
        self.files = [f for f in files]
        self.filenames = [f for f in filenames]
        self.dictionaries = [d for d in dictionaries]
        self.directories = [d for d in directories]

    def __enter__(self):
        pass

    def __exit__(self, type, value, traceback):
        files = self.files

        # extend files with all values from all dictionaries we're tracking
        for d in self.dictionaries:
            files.extend(d.values())

        for f in files:
            name = f.name if hasattr(f, 'name') else str(f)
            logger.debug("Cleaning up %s", name)
            try:
                f.close()
            except:
                pass

            try:
                os.remove(f.name)
            except:
                pass

        for name in self.filenames:
            logger.debug("Cleaning up %s", name)
            try:
                os.remove(name)
            except:
                pass

        for dirname in self.directories:
            logger.debug("Removing directory %s", dirname)
            try:
                shutil.rmtree(dirname)
            except:
                pass
