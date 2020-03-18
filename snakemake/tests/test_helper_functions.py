import os
import pytest
import psutil as ps
import helper_functions as hf


def test_get_reference():
    result = hf.get_reference()
    expected = 'data/chr22_44-45Mb.GRCh37.fasta'
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


def test_get_svtype():
    result = hf.get_svtype()
    expected = 'indel'
    assert result == expected


def test_get_nthreads():
    result = hf.get_nthreads()
    expected = ps.cpu_count(logical=True)
    assert result == expected


def test_get_mem():
    result = hf.get_mem()
    expected = int(ps.virtual_memory().free / hf.get_nthreads() / 2**20)
    assert result == expected
