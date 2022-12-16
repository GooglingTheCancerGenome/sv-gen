"""Test suite for validator/__init__.py."""
import pytest

from validator import load_configfile
from yatiml import RecognitionError


def test_simulation__svtype_zerocount_exception(shared_datadir):
    with pytest.raises(RecognitionError):
        load_configfile(shared_datadir / 'analysis_svtype_zerocount_error.yaml')


def test_simulation__svtype_paramcount_exception(shared_datadir):
    with pytest.raises(RecognitionError):
        load_configfile(shared_datadir / 'analysis_svtype_paramcount_error.yaml')


def test_simulation__svtype_paramtype_exception(shared_datadir):
    with pytest.raises(RecognitionError):
        load_configfile(shared_datadir / 'analysis_svtype_paramtype_error.yaml')

