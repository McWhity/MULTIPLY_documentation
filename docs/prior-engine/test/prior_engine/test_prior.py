# import sys
import os
import pytest
import sys

myPath = os.path.dirname(os.path.abspath(__file__))
sys.path.append(myPath + '/../../multiply_prior_engine/')
from multiply_prior_engine.prior_creator import PriorCreator
from multiply_prior_engine.prior_engine import PriorEngine
from multiply_prior_engine.soilmoisture_prior_creator import SoilMoisturePriorCreator


def test_priorengine_init():
    P = PriorEngine(config='./test/prior_engine/test_config_prior.yml',
                    datestr='2017-01-01',
                    variables=['sm'])
    assert P.configfile is not None
    assert type(P.configfile) is str


def test_priorengine_get_priors_sm():
    P = PriorEngine(config='./test/prior_engine/test_config_prior.yml',
                    datestr='2017-03-01',
                    variables=['sm'])
    assert type(P.get_priors()) is dict


def test_priorengine_get_priors_lai():
    P = PriorEngine(config='./test/prior_engine/test_config_prior.yml',
                    datestr='2017-01-01',
                    variables=['lai'])
    assert type(P.get_priors()) is dict


def test_priorengine_get_priors_n():
    P = PriorEngine(config='./test/prior_engine/test_config_prior.yml',
                    datestr='2017-01-01',
                    variables=['n'])
    assert type(P.get_priors()) is dict


def test_sm_prior_init():
    with pytest.raises(AssertionError,
                       message=("Expecting AssertionError \
                                --> no config specified")):
        SoilMoisturePriorCreator()


def test_sm_prior_no_ptype():
    with pytest.raises(AssertionError,
                       message=("Expecting AssertionError \
                                --> no config specified")):
        SoilMoisturePriorCreator()


def test_sm_prior_invalid_ptype():
    with pytest.raises(AssertionError,
                       message=("Expecting AssertionError \
                                --> no config specified")):
        SoilMoisturePriorCreator(ptype='climatologi')


def test_calc_variable():
    with pytest.raises(AssertionError,
                       message=("Expecting AssertionError \
                                --> no variable specified")):
        P = PriorEngine(config='./test/prior_engine/test_config_prior.yml',
                        datestr='2017-03-01')
        P.get_priors()


"""

Test prior types:

"""


def test_climatology_prior():
    P = PriorEngine(config='./test/prior_engine/test_config_prior.yml',
                    datestr='2017-03-01',
                    variables=['sm'])
    SoilMoisturePriorCreator(config=P.config,
                             datestr='2017-01-01',
                             ptype='climatology',
                             var=["sm"])


def test_recent_prior():
    P = PriorEngine(config='./test/prior_engine/test_config_recent_prior.yml',
                    datestr='2017-03-01',
                    variables=['sm'])
    with pytest.raises(AssertionError,
                       message=("Expecting AssertionError"),
                       match=r'.*recent.*'):
        P.get_priors()


def test_recent_prior2():
    P = PriorEngine(config='./test/prior_engine/test_config_recent_prior.yml',
                    datestr='2017-03-01',
                    variables=['sm'])
    with pytest.raises(AssertionError,
                       message=("Expecting AssertionError \
                                --> recent prior not implemented")):
        S = SoilMoisturePriorCreator(config=P.config,
                                     datestr='2017-03-01',
                                     var="sm",
                                     ptype='recent')
        S.compute_prior_file()


def test_user_prior_initialization():
    P = PriorEngine(config='./test/prior_engine/test_config_user_prior.yml',
                    datestr='2017-03-01',
                    variables=['sm'])
    S = SoilMoisturePriorCreator(config=P.config,
                                 datestr='2017-03-01',
                                 var="sm",
                                 ptype='user1')
    with pytest.raises(AssertionError,
                       message=("Expecting AssertionError \
                                --> wrong path in config file."),
                       match=r'.*path/to/file*'):
        S.compute_prior_file()


def test_calc_config():
    P = PriorEngine(config='./test/prior_engine/test_config_prior.yml',
                    datestr='2017-03-01',
                    variables=['sm'])
    S = SoilMoisturePriorCreator(config=P.config,
                                 datestr='2017-03-01',
                                 ptype='climatology',
                                 var="sm")
    assert type(S.config) is dict


def test_sm_prior_missing_datestr():
    P = PriorEngine(config='./test/prior_engine/test_config_prior.yml',
                    datestr='2017-03-01',
                    variables=['sm'])
    with pytest.raises(AssertionError,
                       message=("Expecting AssertionError \
                                --> no datestr & variable specified")):
        SoilMoisturePriorCreator(config=P.config,
                                 ptype='climatology')

# def test_roughness():
#     lut_file = tempfile.mktemp(suffix='.lut')
#    lc_file = tempfile.mktemp(suffix='.nc')
#     gen_file(lut_file)
#     gen_file(lc_file)
#     P = RoughnessPrior(ptype='climatology',
#                        lut_file=lut_file,
#                        lc_file=lc_file)
