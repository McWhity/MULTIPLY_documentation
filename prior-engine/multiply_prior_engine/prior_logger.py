#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
    Logger for the Prior Engine for MULTIPLY.

    Copyright (C) 2018  Thomas Ramsauer
"""

import os
import yaml
import logging
import logging.config
import tempfile


__author__ = "Thomas Ramsauer"
__copyright__ = "Thomas Ramsauer"
__maintainer__ = "Thomas Ramsauer"
__email__ = "t.ramsauer@iggf.geo.uni-muenchen.de"


class PriorLogger(object):
    """ Logger for the Prior Engine for MULTIPLY.

    initializes logging.
    By default the root logger is initialized with:
      - level: WARNING;
      - handlers: console & file.

    The log file will be placed in the 'tempfile.tempdir'.

    """

    configfile = os.path.join(os.path.dirname(__file__),
                              'prior_engine_logging.yml')

    def __init__(self,
                 level: str='warning',
                 handlers: list=['console', 'file']):

        # get logging config
        with open(self.configfile, 'r') as cfg:
            config_dict = yaml.load(cfg)

        # add temporary log file name to configuration
        filename = os.path.join(tempfile.tempdir, 'prior_engine.log')
        config_dict['handlers']['file']['filename'] = filename

        # set level
        assert level.lower() in ['debug', 'info', 'warning',
                                 'error', 'critical'],\
            '{} no defined logger level.'.format(level)
        config_dict['root']['level'] = level.upper()

        # set handlers
        if not isinstance(handlers, list):
            raise TypeError('PriorLogger: Handlers must be passed as list.')
        for h in handlers:
            assert h.lower() in ['console', 'file'],\
                '{} not a valid handler.'.format(h)
        config_dict['root']['handlers'] = handlers

        # finally configure logging
        logging.config.dictConfig(config_dict)

        logging.info('------------- Logger initialized. -------------')
        logging.info('The log file can be found at {}'.format(filename))

    # def return_logger(self):
        # return logging.getLogger(__name__)
        # return self.logger
