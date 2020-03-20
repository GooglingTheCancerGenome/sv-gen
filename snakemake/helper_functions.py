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
    fname = config.input.fasta
    fext = config.filext.fasta
    if not os.path.exists(fname):
        raise FileNotFoundError("FASTA file '{}' not found.".format(fname))
    if not fname.endswith(fext):
        raise ValueError("FASTA file extension '{}' not registered.".format(
            os.path.splitext(fname)[-1]))
    if os.path.getsize(fname) == 0:
        raise OSError("FASTA file '{}' is empty.".format(fname))
    return fname


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
    mem = ps.virtual_memory().free / get_nthreads() / 2**20
    return int(mem)
