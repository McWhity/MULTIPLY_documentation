#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
    Prior Class for MULTIPLY.

    Copyright (C) 2018  Thomas Ramsauer
"""


from abc import ABCMeta, abstractmethod
import datetime
from dateutil.parser import parse
from typing import List
import numpy as np


__author__ = ["Alexander Löw", "Thomas Ramsauer"]
__copyright__ = "Copyright 2018, Thomas Ramsauer"
__credits__ = "Alexander Löw"
__maintainer__ = "Thomas Ramsauer"
__email__ = "t.ramsauer@iggf.geo.uni-muenchen.de"


class PriorCreator(metaclass=ABCMeta):

    def __init__(self, **kwargs):
        self.ptype = kwargs.get('ptype', None)
        self.config = kwargs.get('config', None)
        self.datestr = kwargs.get('datestr', None)
        self.variable = kwargs.get('var', None)
        self._check()
        self.time_vector, self.time_vector_months = self._create_time_vector()
        self.date = self._create_datetime()
        self.date8 = int(str(self.date.date()).replace('-', ''))

    def _check(self):
        assert self.ptype is not None, 'No prior type specified.'
        # TODO make use of config optional
        assert self.config is not None, 'No config available.'
        assert self.datestr is not None, 'No datestr available.'
        assert self.variable is not None, 'No variable available.'

    def _create_time_vector(self):
        """Creates a time vector dependent on start & end time and time interval
        from config file.
        A vector containing datetime objects is written to self.time_vector.
        A vector containing months ids (1-12) for each timestep is written to
        self.time_vector_months.

        :returns: -
        :rtype: -
        """
        date_format = ('%Y-%m-%d')
        s = self.config['General']['start_time']
        e = self.config['General']['end_time']
        interval = self.config['General']['time_interval']
        if type(s) is str:
            s = datetime.datetime.strptime(s, date_format)
        if type(e) is str:
            e = datetime.datetime.strptime(e, date_format)
        t_span = (e-s).days + 1

        time_vector = [(s+(datetime.timedelta(int(x))))
                       for x in np.arange(0, t_span, interval)]

        # create list of month ids for every queried point in time:
        time_vector_months = [(s+(datetime.timedelta(int(x)))).month
                              for x in np.arange(0, t_span, interval)]
        return time_vector, time_vector_months

    def _create_datetime(self):
        # parse (dateutil) self.datestr to create datetime.datetime object
        try:
            date = parse(self.datestr)
        except TypeError as e:
            print('[WARNING] No time info was passed to Prior Class!')
            return
        # assert parsing of self.date is working.
        assert type(date) is datetime.datetime,\
            'could not parse date {}'.format(date)

        # get month id/number from self.datestr
        # self.date_month_id = self.date.month
        return date

    @abstractmethod
    def compute_prior_file(self) -> str:
        """
        Might perform some computation, then retrieves the path to a file containing the prior info
        :return:
        """

    @classmethod
    @abstractmethod
    def get_variable_names(cls) -> List[str]:
        """
        :return: A list of the variables that this prior creator is able to create priors for
        """
