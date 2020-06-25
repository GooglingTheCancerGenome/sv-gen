"""Test suite for validator/__init__.py."""
import pytest
import pyfaidx

from validator import load_configfile, SeqIdCountError, SeqIdNotFoundError, \
                      SimulateSvTypeError, SimulateTransError, FastaFilextError


def test_input__countseqid_exception():
    with pytest.raises(SeqIdCountError):
        load_configfile('tests/analysis_seqidcount_error.yaml')


def test_input__seqidnotfound_exception():
    with pytest.raises(SeqIdNotFoundError):
        load_configfile('tests/analysis_seqidnotfound_error.yaml')


def test_input__fastafilext_exception():
    with pytest.raises(FastaFilextError):
        load_configfile('tests/analysis_fastafilext_error.yaml')


def test_simulation__svtype_exception():
    with pytest.raises(SimulateSvTypeError):
        load_configfile('tests/analysis_svtype_error.yaml')


def test_simulation__translocation_exception():
    with pytest.raises(SimulateTransError):
        load_configfile('tests/analysis_trans_error.yaml')
