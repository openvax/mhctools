from __future__ import print_function, division, absolute_import

import logging
import logging.config
import pkg_resources


def get_logger(name):
    logging.config.fileConfig(pkg_resources.resource_filename('mhctools', 'logging.conf'))
    return logging.getLogger(name)
