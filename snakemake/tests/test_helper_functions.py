"""Test suite for helper_functions.py."""
import pytest
import psutil as ps
import helper_functions as hf


def test_get_reference():
    result = hf.get_reference()
    expected = 'data/test.fasta'
    assert result == expected


def test_get_reference__filenotfound_exception():
    with pytest.raises(Exception):
        hf.config.input.fasta = 'data/nowhere.fasta'
        hf.get_reference()


def test_get_reference__filextnotfound_exception():
    with pytest.raises(Exception):
        hf.config.filext.fasta = '.fa'
        hf.get_reference()


def test_get_genotype():
    result = hf.get_genotype()
    expected = ['hmz', 'hmz-sv', 'htz-sv']
    assert set(result) == set(expected)


@pytest.fixture
def set_svtype():
    conf = hf.config.simulation.svtype
    conf.indel.count = 10  # default
    conf.tra.count = 10
    return conf


def test_get_svtype(set_svtype):
    hf.config.simulation.svtype = set_svtype
    result = hf.get_svtype()
    expected = 'indel_tra'
    assert result == expected


def test_get_nthreads():
    result = hf.get_nthreads()
    expected = ps.cpu_count(logical=True)
    assert result == expected


def test_get_mem():
    tolerance = 0.01
    result = hf.get_mem()
    expected = int(ps.virtual_memory().free / hf.get_nthreads() / 2**20)
    assert result == pytest.approx(expected, tolerance)
