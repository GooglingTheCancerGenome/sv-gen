"""Test suite for validator/__init__.py."""
import pytest
import pyfaidx

from validator import load_configfile, SeqIdCountError, SeqIdNotFoundError, \
                      SimulateSvTypeError, SimulateTransError, FastaFilextError


def test_input__countseqid_exception(shared_datadir):
    with pytest.raises(SeqIdCountError):
        load_configfile(shared_datadir / 'analysis_seqidcount_error.yaml')


def test_input__seqidnotfound_exception(shared_datadir):
    with pytest.raises(SeqIdNotFoundError):
        load_configfile(shared_datadir / 'analysis_seqidnotfound_error.yaml')


def test_input__fastafilext_exception(shared_datadir):
    with pytest.raises(FastaFilextError):
        load_configfile(shared_datadir / 'analysis_fastafilext_error.yaml')


def test_simulation__svtype_exception(shared_datadir):
    with pytest.raises(SimulateSvTypeError):
        load_configfile(shared_datadir / 'analysis_svtype_error.yaml')


def test_simulation__translocation_exception(shared_datadir):
    with pytest.raises(SimulateTransError):
        load_configfile(shared_datadir / 'analysis_trans_error.yaml')
