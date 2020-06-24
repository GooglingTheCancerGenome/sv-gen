"""Helper functions."""
import os
import psutil as ps

from validator import load_configfile

config = load_configfile('analysis.yaml')


def get_reference():
    """
    Get reference genome in FASTA format.

    :returns: filepath
    """
    return config.input.fasta


def get_genotype():
    """
    Get a list of genotypes.

    :returns (list) genotypes
    """
    return [str(g.value) for g in list(config.output.genotype)]


def get_svtype():
    """
    Get one or more SV types with non-zero counts.

    :returns: (str) SV type(s)
    """
    types = []
    for sv, params in config.simulation.svtype.__dict__.items():
        if params.count > 0:
            types.append(sv)
    types.sort()
    return '_'.join(types)


def get_nthreads(logical=True):
    """
    Get the number of threads used by `samtools` and `bwa`.

    :returns: (int) threads (default: -1 = number of logical cores)
    """
    n = int(config.threads)
    if n > 0:
        return n
    return ps.cpu_count(logical)


def get_mem():
    """
    Get free memory per core used by `samtools sort`.

    :returns: (int) in MB
    """
    if config.memory < 0:
        return int(ps.virtual_memory().free / get_nthreads() / 2**20)
    return config.memory


def get_tmpspace():
    """
    Get the amount of temporary space used by `samtools sort`.

    :returns (int) in MB
    """
    return config.tmpspace
