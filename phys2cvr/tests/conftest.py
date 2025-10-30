"""
This configuration test module was taken from phys2bids.
Credit to the original author(s) and to the phys2bids community.
"""

import os
import ssl
from urllib.request import urlretrieve

import nibabel as nib
import numpy as np
import pytest


def fetch_file(osf_id, path, filename):
    """
    Fetch file located on OSF and downloads to `path`/`filename`.

    Parameters
    ----------
    osf_id : str
        Unique OSF ID for file to be downloaded. Will be inserted into relevant
        location in URL: https://osf.io/{osf_id}/download
    path : str
        Path to which `filename` should be downloaded. Ideally a temporary
        directory
    filename : str
        Name of file to be downloaded (does not necessarily have to match name
        of file on OSF)

    Returns
    -------
    full_path : str
        Full path to downloaded `filename`
    """
    # This restores the same behavior as before.
    # this three lines make tests downloads work in windows
    if os.name == 'nt':
        orig_sslsocket_init = ssl.SSLSocket.__init__
        ssl.SSLSocket.__init__ = (
            lambda *args, cert_reqs=ssl.CERT_NONE, **kwargs: orig_sslsocket_init(
                *args, cert_reqs=ssl.CERT_NONE, **kwargs
            )
        )
        ssl._create_default_https_context = ssl._create_unverified_context
    url = f'https://osf.io/{osf_id}/download'
    full_path = os.path.join(path, filename)
    if not os.path.isfile(full_path):
        urlretrieve(url, full_path)
    return full_path


@pytest.fixture(scope='session')
def testdir(tmp_path_factory):
    """Test path that will be used to download all files."""
    return tmp_path_factory.getbasetemp()


@pytest.fixture
def atlas(testdir):
    return fetch_file('h6nj7', testdir, 'atlas.nii.gz')


@pytest.fixture
def atlastime(testdir):
    return fetch_file('ts6a8', testdir, 'ats.nii.gz')


@pytest.fixture
def mean_fc(testdir):
    return fetch_file('jrg8d', testdir, 'mean_fc_matlab.tsv')


@pytest.fixture
def sdi(testdir):
    return fetch_file('rs4dn', testdir, 'SDI_matlab.tsv')


@pytest.fixture
def sc_mtx(testdir):
    return fetch_file('vwh75', testdir, 'sc.mat')


@pytest.fixture
def timeseries(testdir):
    return fetch_file('ay8df', testdir, 'func.mat')


@pytest.fixture
def nifti_data(testdir):
    # Create a test nifti file with random data
    data = np.random.rand(10, 10, 10)
    fname = os.path.join(testdir, 'test.nii.gz')
    img = nib.Nifti1Image(data, np.eye(4))
    nib.save(img, fname)
    return fname, data


@pytest.fixture
def nifti_mask(testdir):
    # Create a test nifti mask file with random binary data
    data = np.random.randint(0, 2, size=(10, 10, 10))
    fname = os.path.join(testdir, 'test_mask.nii.gz')
    img = nib.Nifti1Image(data, np.eye(4))
    nib.save(img, fname)
    return fname, data
