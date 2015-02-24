import logging
import os

import shutil

class CleanupFiles(object):
    """
    Context manager that deletes a set of files at the end of a block or
    if an exception gets raised
    """
    def __init__(
            self,
            files = [],
            filenames = [],
            dictionaries = [],
            directories = []):
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
            logging.info("Cleaning up %s", name)
            try:
                f.close()
            except:
                pass

            try:
                os.remove(f.name)
            except:
                pass

        for name in self.filenames:
            logging.info("Cleaning up %s", name)
            try:
                os.remove(name)
            except:
                pass

        for dirname in self.directories:
            logging.info("Removing directory %s", dirname)
            try:
                shutil.rmtree(dirname)
            except:
                pass
