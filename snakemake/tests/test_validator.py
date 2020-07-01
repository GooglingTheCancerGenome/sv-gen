"""Test suite for validator/__init__.py."""
import pytest

from validator import load_configfile
from yatiml import RecognitionError


def test_seqids_empty(shared_datadir):
    conf = load_configfile(shared_datadir / 'analysis_seqids_empty.yaml')
    assert conf._seqids() == 2


def test_input__seqids_unknown_exception(shared_datadir):
    with pytest.raises(ValueError):
        load_configfile(shared_datadir / 'analysis_seqids_unknown_error.yaml')


def test_input__fastafilext_exception(shared_datadir):
    with pytest.raises(ValueError):
        load_configfile(shared_datadir / 'analysis_fastafilext_error.yaml')


def test_simulation__svtype_zerocount_exception(shared_datadir):
    with pytest.raises(ValueError):
        load_configfile(shared_datadir / 'analysis_svtype_zerocount_error.yaml')


def test_simulation__svtype_paramcount_exception(shared_datadir):
    with pytest.raises(RecognitionError):
        load_configfile(shared_datadir / 'analysis_svtype_paramcount_error.yaml')


def test_simulation__svtype_paramtype_exception(shared_datadir):
    with pytest.raises(RecognitionError):
        load_configfile(shared_datadir / 'analysis_svtype_paramtype_error.yaml')


def test_simulation__svtype_tra_exception(shared_datadir):
    with pytest.raises(ValueError):
        load_configfile(shared_datadir / 'analysis_svtype_tra_error.yaml')
