#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
    Soil Priors for Prior Engine in MULTIPLY.

    Copyright (C) 2018  Thomas Ramsauer
"""


import datetime
import glob
import os
import subprocess
import tempfile
import logging

import numpy as np
import shapely
import shapely.wkt
from netCDF4 import Dataset
from osgeo import gdal
from scipy import spatial

from .prior_creator import PriorCreator


__author__ = ["Alexander Löw", "Thomas Ramsauer"]
__copyright__ = "Copyright 2018, Thomas Ramsauer"
__credits__ = "Alexander Löw"
__maintainer__ = "Thomas Ramsauer"
__email__ = "t.ramsauer@iggf.geo.uni-muenchen.de"


class SoilMoisturePriorCreator(PriorCreator):
    """
    Soil moisture prior class.
    Calculation of climatological prior.
    """

    def __init__(self, **kwargs):
        super(SoilMoisturePriorCreator, self).__init__(**kwargs)

    @classmethod
    def get_variable_names(cls):
        return ['sm']

    def compute_prior_file(self):
        """
        Initialize prior specific (climatological, ...) calculation.

        :returns: filename of prior file
        :rtype: string
        """
        # set None as old data_dir/file may be present from loop.
        self.data_dir = None
        self.data_file = None
        self.output_directory = self.config['Prior']['output_directory']
        if not os.path.exists(self.output_directory):
            os.mkdir(self.output_directory)

        if self.ptype == 'climatology' or self.ptype == 'coarse':
            try:
                data_dir = self.config['Prior']['sm'][self.ptype]['dir']
                self.data_dir = data_dir
                assert os.path.isdir(self.data_dir), \
                    ('Directory does not exist or cannot be found: {}'
                     .format(self.data_dir))
            except KeyError as e:
                assert self.data_dir is not None, \
                  ('Cannot find directory information for '
                   '"{}" prior in config file!'.format(self.ptype))
            else:
                return self._provide_prior_file()

        elif 'user' in self.ptype:
            try:
                data_file = (self.config['Prior']['sm'][self.ptype]['file'])
                self.data_file = data_file
            except KeyError as e:
                assert self.data_file is not None, \
                  ('Cannot find file name for '
                   '"{ptype}" prior in config file (under'
                   ' \'Prior/sm/{ptype}/file:\')!'.format(ptype=self.ptype))
            else:
                return self._provide_prior_file()

        # TODO recent only place holder
        elif self.ptype == 'recent':
            self._get_recent_sm_proxy()

        else:
            msg = '{} prior for sm not implemented'.format(self.ptype)
            logging.exception(msg)
            assert False, msg
        return self._provide_prior_file()

    def _calc_climatological_prior(self):
        """
        Calculate climatological prior.
        Reads climatological file and extracts proper values for given
        timespan and -interval.
        Then converts the means and stds to state vector and covariance
        matrices.

        :returns: state vector and covariance matrix
        :rtype: tuple

        """
        self._get_climatology_file()
        self._extract_climatology()

        # TODO limit to months

        # date_format = ('%Y-%m-%d')
        s = self.config['General']['start_time']
        e = self.config['General']['end_time']
        interval = self.config['General']['time_interval']
        t_span = (e - s).days + 1
        # print(t_span)

        # create list of month ids for every queried point in time:
        idt = [(s + (datetime.timedelta(int(x)))).month
               for x in np.arange(0, t_span, interval)]
        # idt_unique = list(set(idt))

        # create nd array with correct dimensions (time, x, y):
        p = np.ndarray(shape=(len(idt), self.clim.shape[1],
                              self.clim.shape[2]),
                       dtype=float)
        std = p.copy()

        # read correspending data into mean and std arrays:
        for i in range(len(idt)):
            p[i, :, :] = self.clim[idt[i] - 1, :, :]
            std[i, :, :] = self.std[idt[i] - 1, :, :]

        # calculate uncertainty with normalization via coefficient of variation
        # TODO scale uncertainty
        sm_unc = (std / np.mean(self.clim))
        # inverse covariance matrix
        diagon = (1. / sm_unc)
        # print(diagon.shape)

        # def create_sparse_matrix(a):
        #     return sp.sparse.lil_matrix(np.eye(t_span)*a)

        # C_prior_inv = np.apply_along_axis(create_sparse_matrix, 0, diagon)
        C_prior_inv = diagon

        # DISCUSS TODO
        # rather write to self.'prior_key' to easy concatenate afterwards
        # via concat_priors.

        return p, C_prior_inv

    def _get_climatology_file(self):
        """
        Load pre-processed climatology into self.clim_data.
        Part of prior._calc_climatological_prior().

        """
        assert (self.config['Prior']['sm']['climatology']
                           ['climatology_file']) is not None,\
            'There is no climatology file specified in the config!'

        # use xarray:
        # self.clim_data = xr.open_dataset(self.config['priors']['sm_clim']
        self.clim_data = Dataset(self.config['Prior']['sm']['climatology']
                                 ['climatology_file'])

    def _provide_prior_file(self):
        """Provide variable and prior type specific prior file name to Prior Engine.

        :returns: absolute path to prior file for requested prior.
        The file is gdal-compatible to be used in inference engine - either
        GeoTiff or VRT format.
        It includes 2 bands:
         1. mean value raster
         2. uncertainty raster
        :rtype: string

        """

        if self.data_dir is not None:
            self.data_file = self._get_prior_file_from_dir(self.data_dir)
        else:
            assert self.data_file is not None

        ext = os.path.splitext(self.data_file)[-1].lower()
        if ext == 'vrt':
            logging.info('Prior file ({}, {}) is already a .vrt-file, no need '
                         'to convert.'.format(self.variable, self.ptype))
        else:
            try:
                logging.info('Trying to convert prior file to .vrt-format.')
                self.data_file = self._create_global_vrt(self.data_file)
            except:
                assert False, ("Could not create .vrt-file for {} {} prior"
                               " ({})".format(self.variable, self.ptype,
                                              self.data_file))
        return self.data_file

    def _get_prior_file_from_dir(self, directory, return_vrt=True):
        """Get filename(s) of prior file(s) from directory.
        If multiple files are found self._merge_multiple_prior_files is called.

        Currently, the following prior types are supported:
        - climatology (calculated from ESA CCI data, standard)
        - coarse (daily aggregated SMAP L4 data, standard)
        - user prior, provided through user_prior_creator

        :param directory: directory containing the files (from config)
        :returns: filename
        :rtype: string

        """
        fn = None
        if self.ptype == 'climatology':
            pattern = (r"ESA_CCI_SM_clim_{:02d}.tiff"
                       .format(self.date.month))
        elif self.ptype == 'coarse':
            pattern = (r"SMAP_daily_{:8d}.tif"
                       .format(self.date8))
        # TODO read user pattern from config file to allow defined input
        # (has to be written to the config-file in a 'config step' first)
        elif 'user' in self.ptype:
            pattern = (r"user_{}.tiff$")
        elif self.ptype == 'recent':
            pattern = (r"recent_prior_{}.tiff$"
                       .format(self.date8))
        else:
            pattern = (r"*")

        fn_list = sorted(glob.glob('{}'.format(os.path.join(
            os.path.abspath(directory), pattern), recursive=True)))

        # AssertionError is caught by the prior engine:
        assert fn_list is not None and len(fn_list) > 0, \
            ('Did not find {} {} '
                'prior files in {} (pattern: \'{}\')!'
                .format(self.variable, self.ptype,
                        os.path.abspath(directory), pattern))
        if len(fn_list) > 1:
            fn = self._merge_multiple_prior_files(fn_list)
        else:
            fn = fn_list[0]

        self._check_gdal_compliance(fn)
        return '{}'.format(fn)

    def _merge_multiple_prior_files(self, fn_list):
        """Merge files if more than one is available for current time step.
        should be obsolete.

        :param fn_list: file list to process
        :returns: file name of merged file
        :rtype: string

        """
        # create list of alphabet for gdal funciton call
        abc = [chr(i) for i in range(ord('A'), ord('Z')+1)]
        mean_instr, unc_instr, calc_instr = '', '', ''
        # create input strings for gdal calculate call
        for i, f in enumerate(fn_list):
            self._check_gdal_compliance(f)
            mean_instr += ('-{abc} {fn} --{abc}=1'
                           .format(abc=abc[i], fn=f))
            unc_instr += ('-{abc} {fn} --{abc}=2'
                          .format(abc=abc[i], fn=f))
        calc_instr = '+'.join(map(str, abc[:i+1]))
        # create temporary files to write mean mean&unc to
        mean_tf = tempfile.NamedTemporaryFile(suffix='_mean.vrt')
        unc_tf = tempfile.NamedTemporaryFile(suffix='_unc.vrt')
        # create means of input file mean and uncertainty files
        # TODO replace subprocess call.. rasterio?
        subprocess.run('gdal_calc.py {} --outfile={} --overwrite '
                       '--calc="({})/{}"'
                       .format(mean_instr, mean_tf.name, calc_instr,
                               str(len(fn_list))),
                       shell=True, check=True)
        subprocess.run('gdal_calc.py {} --outfile={} --overwrite '
                       '--calc="({})/{}"'
                       .format(unc_instr, unc_tf.name, calc_instr,
                               str(len(fn_list))),
                       shell=True, check=True)
        # write combined/averaged mean&uncertainty info to generic file
        directory = os.path.abspath(os.path.dirname(f))
        out_fn = os.path.join(directory, self.ptype+self.date8)
        out_vrt = gdal.BuildVRT(out_fn+'.vrt', [mean_tf.name, unc_tf.name],
                                separate=True)
        out_ds = gdal.Translate(out_fn+'.tif', out_vrt)
        logging.info('Created {} from following files: {}.'
                     .format(out_fn+'.vrt', fn_list))
        # close/delete temporary files
        mean_tf.close()
        unc_tf.close()
        return out_ds

    def _create_global_vrt(self, fn, local=True):
        """Create VRT file for file.

        By default, the .vrt-file will be written to a local temporary
        directory. If `local` is set to False, the file is written to the
        directory the input file (fn) currently lives in.

        :param fn: file name
        :param local: create temporary local vrt.
        :returns: file name of created vrt, or initial file name if no success.
        :rtype: string

        """
        # TODO should it be an option in config if vrt is created and where?
        logging.info('Creating vrt file from {}.'.format(fn))
        self._check_gdal_compliance(fn)
        try:
            temp_fn = ('{}_prior_{}_{}.vrt'
                       .format(self.variable,
                               self.ptype,
                               self.date8))
            out_fn = os.path.join(self.output_directory, temp_fn)
            vrt_options = gdal.BuildVRTOptions(outputBounds=(-180, -90, 180, 90))
            gdal.BuildVRT(out_fn, fn, options=vrt_options)
            res = '{}'.format(out_fn)
            assert os.path.isfile(res), "{} is not a file.".format(res)
            self._check_gdal_compliance(res)
            return res
        except Exception as e:
            logging.warning('Cannot create .vrt file'
                           ' {} - returning {}.'.format(res, fn))
            return '{}'.format(fn)

    def _check_gdal_compliance(self, fn):
        try:
            ds = gdal.Open(fn)
            assert ds is not None, \
                ('GDAL: Check: Cannot open file ({})'.format(fn))
            ds = None
        except AssertionError as e:
            logging.error(e)
            raise

    def _extract_climatology(self):
        """
        Extract climatology values for ROI.
        Part of _clac_climatological_prior().

        """
        clim = self.clim_data.variables['sm'][:]
        std = self.clim_data.variables['sm_stdev'][:]
        lats = self.clim_data.variables['lat'][:]
        lons = self.clim_data.variables['lon'][:]

        ROI_wkt = shapely.wkt.loads(self.config['General']['roi'])

        # minx, miny, maxx, maxy:
        minx, miny, maxx, maxy = ROI_wkt.bounds
        ROI = [(minx, miny), (maxx, maxy)]

        # TODO insert check for ROI bounds:
        # if POI[0] < np.amin(lats) or POI[0] > np.amax(lats) or\
        #    POI[1] < np.amin(lons) or POI[1] > np.amax(lons):
        #     raise ValueError("POI's latitude and longitude out of bounds.")

        # stack all raveled lats and lons from climatology
        # combined_LAT_LON results in e.g.
        # [[ 34.625 -10.875]
        #  [ 34.625 -10.625]
        #  [ 34.625 -10.375]...]
        combined_LAT_LON = np.dstack([lats.ravel(), lons.ravel()])[0]
        mytree = spatial.cKDTree(combined_LAT_LON)
        idx, idy = [], []
        # for all coordinate pairs in ROI search indexes in combined_LAT_LON:
        for i in range(len(ROI)):
            dist, indexes = mytree.query(ROI[i])
            x, y = tuple(combined_LAT_LON[indexes])
            idx.append(indexes % clim.shape[2])
            idy.append(int(np.ceil(indexes / clim.shape[2])))

        # TODO check assignment
        # print(idx, idy)

        # extract sm data
        sm_area_ = clim[:, min(idy):max(idy) + 1, :]
        sm_area = sm_area_[:, :, min(idx):max(idx) + 1]

        # extract sm stddef data
        sm_area_std_ = std[:, min(idy):max(idy) + 1, :]
        sm_area_std = sm_area_std_[:, :, min(idx):max(idx) + 1]
        # sm_area_std = np.std(sm_area, axis=(1, 2))
        # sm_area_mean = np.mean(sm_area, axis=(1, 2))

        # TODO Respect spatial resolution in config file.
        #      Adjust result accordingly.

        # print(sm_area)
        self.clim = sm_area
        self.std = sm_area_std

    def _get_recent_sm_proxy(self):
        assert False, "recent sm proxy not implemented"


class MapPriorCreator(PriorCreator):
    """
    Prior which is based on a LC map and a LUT
    """

    def __init__(self, **kwargs):
        """
        Parameters
        ----------
        lut_file : str
            filename of LUT file
        lc_file : str
            filename of landcover file
        """
        super(MapPriorCreator, self).__init__(**kwargs)
        self.lut_file = kwargs.get('lut_file', None)
        assert self.lut_file is not None, 'LUT needs to be provided'

        self.lc_file = kwargs.get('lc_file', None)
        assert self.lc_file is not None, 'LC file needs to be provided'

        # check that files exist
        assert os.path.exists(self.lc_file)
        assert os.path.exists(self.lut_file)


class RoughnessPriorCreator(MapPriorCreator):

    def __init__(self, **kwargs):
        super(RoughnessPriorCreator, self).__init__(**kwargs)

    def calc(self):
        if self.ptype == 'climatology':
            self._read_lut()
            self._read_lc()
            self._map_lut()
            self.file = self.save()
        else:
            assert False

    def _read_lut(self):
        self.lut = 'abc'

    def _read_lc(self):
        self.lc = 'efg'

    def _map_lut(self):
        """
        should do the mapping of s, l, ACL type
        """
        self.roughness = self.lut + self.lc

    def save(self):
        """
        save mapped roughness data to file
        """
        return tempfile.mktemp(suffix='.nc')

    @classmethod
    def get_variable_names(cls):
        return ['roughness']

    def compute_prior_file(self):
        assert False, 'roughness prior not implemented'
