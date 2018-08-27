#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
    Integrate user defined priors to MULTIPLY engine.

    Copyright (C) 2018  Thomas Ramsauer
"""

import argparse
import csv
import datetime
import logging
import os
import sys
import tempfile
import warnings
import pandas as pd
import shutil

import yaml

from .prior import PriorCreator
from .prior_engine import _get_config, PriorEngine

__author__ = "Thomas Ramsauer"
__copyright__ = "Copyright 2018 Thomas Ramsauer"
__maintainer__ = "Thomas Ramsauer"
__email__ = "t.ramsauer@iggf.geo.uni-muenchen.de"


class UserPriorCreator(PriorCreator):
    """

    """

    def __init__(self, **kwargs):
        super(UserPriorCreator, self).__init__(**kwargs)

    def RetrievePrior(self):
        """
        Initialize prior specific (climatological, ...) calculation.

        :returns:

        """

    def userprior_conversion(self):
        """Convert user defined data for compatibility?

        :returns:
        :rtype:

        """
        pass


class UserPriorInput(object):
    """
    UserPriorInput Class

    This class contains methods to 'add', 'remove' information to the prior
    configuration and an 'import' method to import prior data provided by the
    user.

    When adding a prior config to the configuration file, the script checks if
    other 'user' priors are abundant and adds suitable numbering accordingly.
    The new configuration information is written to a subsection for the
    specified variable and labled e.g. 'user1'.

    Deleting means removing the prior specific information from the
    configuration file.

    ! Currently only tabular data is supported!


    """
    def __init__(self, **kwargs):
        # config file so far only needed to verify that variables, which
        # prior information should be added for, are inferrable.
        self.configfile = None
        while self.configfile is None:
            self.configfile = kwargs.get('config', None)
            self.configfile = kwargs.get('configfile', None)
            # have a backup/default config:
            self.configfile = PriorEngine.default_config

        assert self.configfile is not None, \
            ('No configuration filename passed to function. '
             'Add \'configfile=\' to method call.')
        self.config = _get_config(self.configfile)

        self.default_variables_lower = default_variables_lower
        assert self.default_variables_lower is not None, \
            'UserPriorInput does not know about possibly to infer variables.'

    def now(self):
        return datetime.datetime.now().strftime('%Y%m%d_%H%M%S')

    def check_for_user_config(self):
        """Checks subsection of user prior config for specific information
        and extracts them (e.g. directory).

        :returns: -
        :rtype: -

        """
        self.userconf = self.config['Prior'][self.variable][self.ptype]
        for k in self.userconf.keys():
            if 'dir' in k:
                self.dir = k
            if 'other_options' in k:
                pass

    def _path_exists(self, path):
        try:
            assert os.path.exists(path), \
                ('Entered path ({}) does not exist!'.format(path))
            return path
        except AssertionError as e:
            # TODO creation of new folder?
            try:
                parser.error(e)
            except:
                raise(e)

    def _count_defined_userpriors(self, configfile):
        """Checks if a user defined prior is already defined in config, to:
        Writing user prior config:
           get count so that the next section can be written (user1, user2,...)

        :returns: count of user defined prior sections.
        :rtype: int

        """
        self.config = _get_config(configfile)
        # TODO add possibility for other names. these must be stored somewhere.
        existing_userprior = 0
        logging.debug('Existing ptypes:')
        for ptype in self.config['Prior'][self.variable].keys():
            logging.debug(ptype)
            if 'user' in ptype:
                existing_userprior += 1
        return existing_userprior

    def _generate_userconf(self, configfile: str, new_configuration: dict):
        """ generate dictionary with user prior information.

        :param configfile: filename of current configuration file
        :param: new_configuration: dictionary holding new prior information
        :returns: -
        :rtype: -

        """
        # configfile_name = kwargs.get('configfile_name', None)
        # if configfile_name is None:
        #     configfile_name = 'user_config.yml'

        # add information to config dictionary

        count = self._count_defined_userpriors(configfile)
        self.config['Prior'][self.variable].update(
            {'user{}'.format(count+1):
                new_configuration
             })

    def write_config(self, configuration, **kwargs):
        """Write configuration to a YAML file.

        :param configuration: configuration dictionary to write to file.
        :Keyword Arguments:
            * *path_to_config* (``str``) --
              path to config file. if None, a tempfile will be created.
            * *new_config_filename* (``str``) --
              Filename of new user config. Only has effect if path_to_config
              is specified.If None, a temporary filename will be used.

        :returns: config file name
        :rtype: string

        """
        path_to_config = kwargs.get('path_to_config', None)
        new_config_filename = kwargs.get('new_config_filename', None)

        if new_config_filename is not None and path_to_config is None:
            warnings.warn('Entered config file name ({}) will be omitted '
                          '--> no path specified!'
                          .format(new_config_filename), Warning)

        self.check_path_to_config_or_create(path_to_config)

        assert os.path.isfile(self.configfile)
        if self.configfile == PriorEngine.default_config:
            # create backup file
            if os.path.exists(PriorEngine.default_config):
                src = os.path.abspath(PriorEngine.default_config)
            a, b = os.path.splitext(src)
            dest = a + '_backup' + b
            logging.info('Creating {}.'.format(dest))
            print('Creating {}.'.format(dest))
            shutil.copyfile(src, dest)
            self.configfile_bakup = dest
        logging.info('User config file: {}'.format(self.configfile))

        with open(self.configfile, 'w') as cfg:
            cfg.write(yaml.dump(configuration, default_flow_style=False))
        return self.configfile

    def check_path_to_config_or_create(self, path_to_config):
        """Create self.configfile variable with path to config file

        :param path_to_config: 
        :returns: 
        :rtype: 

        """
        # check config directory
        if path_to_config is not None:
            path_to_config = self._path_exists(path_to_config)
        else:
            # create temporary files to write config to:
            temp_conf = tempfile.NamedTemporaryFile(
                prefix='PriorEngine_config_{}_'.format(self.now()),
                suffix='.yml',
                delete=False)
            path_to_config = temp_conf.name

        # if valid path entered but is dir
        if not os.path.isfile(path_to_config):
            # and entered new config file name
            if new_config_filename is not None:
                self.configfile = os.path.join(path_to_config,
                                               new_config_filename)
            # but missing new config file name (create generic based on date):
            else:
                self.configfile = os.path.join(
                    path_to_config,
                    'PriorEngine_config_{}.yml'.format(self.now()))
            try:
                # 'x': open for exclusive creation, failing if file exists
                with open(self.configfile, "x") as f:
                    pass
            except FileExistsError as e:
                self.configfile = os.path.join(
                    path_to_config,
                    "PriorEngine_config_{}.yml".format(self.now()))
                with open(self.configfile, "x") as f:
                    warnings.warn(e, Warning)
        # if path is file:
        else:
            self.configfile = path_to_config

    def show_config(self, only_prior=False):
        """Display current prior configuration. Print to stdout and return.

        :returns: (prior engine) configuration
        :rtype: dictionary

        """
        if type(only_prior) is argparse.Namespace:
            p = only_prior.only_prior
        else:
            p = only_prior
        if self.configfile is not None:
            if p:
                print('\nMULTIPLY Prior Engine Configuration \n({}):\n\n{}\n'
                      .format(self.configfile, yaml.dump(self.config['Prior'],
                              default_flow_style=False)))
                return self.config['Prior']
            else:
                print('\nMULTIPLY Configuration \n({}):\n\n{}\n'
                      .format(self.configfile, yaml.dump(self.config,
                              default_flow_style=False)))
                return self.config
        else:
            print('MULTIPLY Configuration file has not been specified yet.'
                  ' Please specify \'configfile=\' when initializing class.')
            # sys.exit()

    def remove_prior(self, variable, ptype, write=True):
        """Remove prior from configuration file.

        :param variable:
        :param ptype:
        :param write:
        :returns:
        :rtype:

        """
        try:
            # removes entry from dictionary and returns it:
            removed = self.config['Prior'][variable].pop(ptype)
            print('Removed {} prior configuration.'.format(removed))
        except KeyError as e:
            warnings.warn('{}/{} not in configuration'
                          .format(variable, ptype), Warning)
        if write:
            self.write_config(self.config)

    def add_prior(self, prior_variable, **kwargs):
        """Adds directory, which holds user prior data, to config file.
        The user defined prior data sets have to be in the common form of gdal
        compatible files (e.g geotiff, vrt, ..). The `import_prior` utility may
        therefor be utilized. The new config will be written to \
        `path_to_config` or `new_config_filename`, please see `write_config`.

        :param prior_variable: variable, which the user prior data is \
                               supporting (e.g. lai, sm)

        :Keyword Arguments:
            * *path_to_config* (``str``) --
              path to config file. if None, a tempfile will be created.
            * *new_config_filename* (``str``) --
              Filename of new user config. Only has effect if path_to_config
              is specified.If None, a temporary filename will be used.
            * *prior_directory* (``str``) --
              Directory path where user prior data is stored.

        :returns:
        :rtype:

        """
        # config file specific info (default ones used if not present):
        path_to_config = kwargs.get('path_to_config', None)
        if path_to_config is None:
            path_to_config = PriorEngine.default_config
        new_config_filename = kwargs.get('new_config_filename', None)

        # so far only directory as user defined configuration implemented
        # TODO needs more flexibility:
        prior_directory = kwargs.get('prior_directory', None)
        user_prior_file = kwargs.get('user_prior_file', None)

        # TODO a date vector for the files has to be passed on or an utility
        # needs to be created to read the correct data for date (an utility to
        # find datestrings in filenames)

        # used in _generate_userconf for location in config file
        self.variable = prior_variable
        assert self.variable.lower() in default_variables_lower

        # check prior data directory
        if prior_directory is not None:
            self.prior_directory = self._path_exists(prior_directory)
        # check user prior file
        if user_prior_file is not None:
            assert os.path.isfile(user_prior_file)

        # adding to new config
        nc = {}
        for arg in kwargs:
            if arg is not 'path_to_config' and\
               arg is not 'new_config_filename':
                nc.update({arg: kwargs[arg]})
        try:
            assert any([x is not None for x in nc.values()]), \
              "No information passed to \'add_prior\' method."
        except AssertionError as e:
            # logging.error(e)
            try:
                parser.error(e)
            except:
                raise e
        # generate new config dictionary with user info included
        self._generate_userconf(configfile=path_to_config,  # to read config
                                new_configuration=nc)  # updates self.config
        # write extended config to file
        self.write_config(path_to_config=path_to_config,
                          new_config_filename=new_config_filename,
                          configuration=self.config)
        self.show_config(only_prior=True)

    def import_prior(self, prior_variable: str, user_file: str, **kwargs):
        """Import user prior data in common MULTIPLY prior data format (gdal
        compatible file, 2 layers).
        Subroutines may be called.

        :param arg:
        :returns:
        :rtype:

        """
        # Check for files or dir of input data, and output directory; create
        # internal directory to store the converted data if not specified

        # config file specific info (default ones used if not present):
        path_to_config = kwargs.get('path_to_config', None)
        new_config_filename = kwargs.get('new_config_filename', None)

        # used in _generate_userconf for location in config file
        self.variable = prior_variable

        # check prior data directory
        if user_file is not None:
            self.user_file = self._path_exists(user_file)
        else:
            msg = 'No user file name specified!'
            try:
                parser.error(msg)
            except:
                raise(msg)
        # -------------

        # Import data with suitable method:
        filename, dtype = os.path.splitext(self.user_file)
        dtype_method = {'csv': _read_tabular,
                        'netCDF': _read_netcdf,
                        'other': _read_other}

        # *** Temporary limitation: ***
        if dtype is 'csv':
            pass
        else:
            assert False, ('Currently, only \'.csv\' files are supported as '
                           'user prior files.')

        # Load Data:
        try:
            dtype_method[dtype](data=self.user_file)
            logging.info('Imported user file {}.'.format(user_file))
        except Exception as e:
            logging.error('Could not import user file {}.'
                         ' Data type {} not (yet) supported.'
                         .format(user_file))
            raise e

        # Convert to gdal compliant file (Inference Engine requirement):
        
        # TODO write to temporary file (geotiff)?!
        # temp_user_file = tempfile.NamedTemporaryFile(
        #                      prefix='User_{}_'.format(self.now()),
        #                      suffix='.')

        # add prior to config
        try:
            self.add_prior(prior_variable=self.variable,
                           path_to_config=path_to_config,
                           new_config_filename=new_config_filename)
            # return 0
        except Exception as e:
            raise e
            # log Error

        # Methods
        def _read_tabular(data):
            with open(data, 'r') as f:
                reader = csv.reader(f)
                variable, lat, lon = tuple(reader[0])
                data = list(reader[1:])
            assert variable.lower() in default_variables_lower, \
                ('Variable {} currently not supported as user prior.'
                 .format(variable))
            self.variable = variable
            self.data = data
            self.latlon = tuple(lat, lon)
            self.data_type = 'point'
            return variable, data, lat, lon

        def _read_netcdf(data):
            # geoval? netCDF4? other? hdf5?
            assert False, 'Reading NetCDF files not implemented yet.'

    def user_prior_cli(self):
        """CLI to include configuration for user defined prior.

        :returns: configfile name (with path)
        :rtype: string

        """

        parser = argparse.ArgumentParser(
            description=('Utility to integrate User Prior data in '
                         'MULTIPLY Prior Engine'),
            prog="user_prior.py",
            # usage='%(prog)s directory [-h] [-p]'
            )
        # action = parser.add_mutually_exclusive_group(required=True)
        # action.add_argument('-A', '--add',
        #                     default=False,
        #                     action='store_true', dest='Add',
        #                     help=('Add user prior data to configuration.'))

        # action.add_argument('-I', '--import',
        #                     default=False,
        #                     action='store_true', dest='Import',
        #                     help=('Import user prior data.'))

        # action.add_argument('-R', '--remove',
        #                     default=False,
        #                     action='store_true', dest='Remove',
        #                     help=('Remove prior data from configuration.'))

        subparsers = parser.add_subparsers()

        parser_Show = subparsers.add_parser('show', aliases='S',
                                            help='Show current prior config.')
        parser_Add = subparsers.add_parser('add', aliases='A', help='Add prior'
                                           ' directory to configuration.')
        parser_Remove = subparsers.add_parser('remove', aliases='R',
                                              help='Remove prior information'
                                              ' from configuration.')
        parser_Import = subparsers.add_parser('import', aliases='I',
                                              help='Import user prior data.')

        # - Show - #

        parser_Show.set_defaults(func=self.show_config)
        parser_Show.add_argument('-p', '--only-prior', default=False,
                                 action='store_true', dest='only_prior',
                                 help=('Only show prior relevant'
                                       ' configuration.'))

        # - Add - #
        parser_Add.set_defaults(func=self.add_prior)
        parser_Add.add_argument('-c', '--path_to_config', type=str,
                                metavar='', required=False,
                                action='store', dest='path_to_config',
                                help=('Directory of new user '
                                      'config.\nIf None, a temporary file'
                                      ' location will be used.'))
        parser_Add.add_argument('-fn', '--new_config_filename', type=str,
                                metavar='', required=False,
                                action='store', dest='new_config_filename',
                                help=('Filename of new user config. Only has'
                                      ' effect if path_to_config is specified.'
                                      '\nIf None, a temporary filename will be'
                                      ' used.'))
        parser_Add.add_argument('-v', '--prior_variable', type=str,
                                metavar='',
                                action='store', dest='prior_variable',
                                required=True,
                                choices=self.default_variables_lower,
                                help=('Variable to use the prior data for.\n'
                                      'Choices are: {}'
                                      .format(self.default_variables_lower)))
        parser_Add.add_argument('-d', '--prior_directory', type=str,
                                metavar='',
                                action='store', dest='prior_directory',
                                help=('Directory which holds specific'
                                      ' user prior data.'))
        parser_Add.add_argument('-f', '--user_prior_file', type=str,
                                metavar='',
                                action='store', dest='user_prior_file',
                                help=('User defined prior file.'))

        # - Remove - #
        parser_Remove.set_defaults(func=self.remove_prior)
        parser_Remove.add_argument('-v', '--prior_variable', type=str,
                                   metavar='',
                                   action='store', dest='prior_variable',
                                   required=True,
                                   choices=self.default_variables_lower,
                                   help=('Variable for which the prior of '
                                         'prior type should be removed '
                                         'from configuration.'
                                         '\nChoices are: {}'.format(
                                             self.default_variables_lower)))
        parser_Add.add_argument('-pt', '--prior_type', type=str,
                                metavar='',
                                action='store', dest='ptype',
                                help=('Prior type. E.g. \'climatological\','
                                      '\'user1\', ... specified in config'
                                      ' file.'))

        # - Import - #
        parser_Import.set_defaults(func=self.import_prior)
        parser_Import.add_argument('-u', '--user_file', type=str,
                                   metavar='', default=None,
                                   action='store', dest='user_file',
                                   help=('User prior file.'))
        parser_Import.add_argument('-v', '--prior_variable', type=str,
                                   metavar='',
                                   action='store', dest='prior_variable',
                                   required=True,
                                   choices=self.default_variables_lower,
                                   help=('Variable to use the prior data for.'
                                         '\nChoices are: {}'.format(
                                             self.default_variables_lower)))

        parser_Import.add_argument('-c', '--path_to_config', type=str,
                                   metavar='', required=False,
                                   action='store', dest='path_to_config',
                                   help=('Directory of new user '
                                         'config.\nIf None, a temporary file'
                                         ' location will be used.'))
        parser_Import.add_argument('-fn', '--new_config_filename', type=str,
                                   metavar='', required=False,
                                   action='store', dest='new_config_filename',
                                   help=('Filename of new user config. '
                                         'Only has effect if path_to_config is'
                                         ' specified.'
                                         '\nIf None, a temporary filename will'
                                         ' be used.'))

        # - -

        args = parser.parse_args()

        try:
            args = vars(args)
            func = args.pop('func')
            func(**args)
        except Exception as e:
            print(e)
            parser.print_usage()

        # actions = [args.Add, args.Remove, args.Import]

        # args.A(Add) will be None if u(user_file) is not provided
        # required_together = ('-A', '-u')
        # if args.Add and args.user_file is None:
        #     msg = "User prior file was not specified!"
        #     logging.error(msg)
        #     parser.error(msg)

        # if all(arg is False for arg in actions):
        #     parser.error('You need to provide an action flag: -A, -I, -R')
        # if sum(actions) > 1:
        #     set_args = [name for a in actions
        #                 for name in vars(args[a])
        #                 if a is not False]
        #     parser.error('Only one action flag can be set. You set {}!'
        #                  .format(set_args))
        # if args.Add:
        #     self.add_prior(**vars(args))

        # if args.Import:
        #     self.import_prior(**vars(args))

        # if args.Remove:
        #     self.remove_prior(**vars(args))


def main():
    try:
        U = UserPriorInput(configfile="./sample_config_prior.yml")
        U.user_prior_cli()
    except ModuleNotFoundError as e:
        print(e)
        # run from outside module or install properly


if __name__ == '__main__':
    main()
