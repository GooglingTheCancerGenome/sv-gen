"""Test suite for validator/__init__.py."""
import pytest
import pyfaidx

from validator import load_configfile


def test_input__seqidcount_exception(shared_datadir):
    with pytest.raises(ValueError):
        load_configfile(shared_datadir / 'analysis_seqidcount_error.yaml')


def test_input__seqidnotfound_exception(shared_datadir):
    with pytest.raises(ValueError):
        load_configfile(shared_datadir / 'analysis_seqidnotfound_error.yaml')


def test_input__fastafilext_exception(shared_datadir):
    with pytest.raises(ValueError):
        load_configfile(shared_datadir / 'analysis_fastafilext_error.yaml')


def test_simulation__svtype_exception(shared_datadir):
    with pytest.raises(ValueError):
        load_configfile(shared_datadir / 'analysis_svtype_error.yaml')


def test_simulation__translocation_exception(shared_datadir):
    with pytest.raises(ValueError):
        load_configfile(shared_datadir / 'analysis_trans_error.yaml')
