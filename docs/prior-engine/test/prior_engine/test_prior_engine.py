import shutil
from multiply_prior_engine import PriorEngine

CONFIG_FILE = './test/prior_engine/prior_engine_test_config.yml'


def test_prior_engine():
    prior_engine = PriorEngine(datestr="2017-03-01",
                               variables=['sm', 'lai', 'cab'],
                               config=CONFIG_FILE)
    priors = prior_engine.get_priors()
    assert 3, len(priors.keys())
    assert 'cab' in priors.keys()
    assert 1 == len(priors['cab'])
    assert 'database' in priors['cab']
    assert ('./test/output_dir/Priors_cab_060_global.vrt' == priors['cab']['database'])
    assert 'lai' in priors.keys()
    assert 1 == len(priors['lai'])
    assert 'database' in priors['lai']
    assert ('./test/output_dir/Priors_lai_060_global.vrt' == priors['lai']['database'])
    assert 'sm' in priors.keys()
    assert 'climatology' in priors['sm']
    assert 'coarse' in priors['sm']
    assert ('./test/output_dir/sm_prior_climatology_20170301.vrt' == priors['sm']['climatology'])

