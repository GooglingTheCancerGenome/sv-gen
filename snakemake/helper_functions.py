import os
import sys
import psutil as ps

from snakemake import load_configfile
from csv import DictReader
from validator import validate_configfile


config_file = validate_configfile('analysis.yaml')
config = load_configfile(config_file)


def get_filext(fmt):
    """Get file extension(s) given file type/format:
        ['fasta', 'fasta_idx', 'bam', 'bam_idx', 'fastq', 'vcf', 'bed']
    :param fmt: (str) input file format
    :returns: (str) file extension
    """
    if fmt not in config['filext'].keys():
         raise ValueError("Input file format '{}' not supported."
            .format(fmt.lower()))
    return config['filext'][fmt]


def get_reference():
    """Get reference genome in FASTA format.
    :returns: filepath
    """
    fname = config['input']['fasta']
    sfx = get_filext('fasta')
    if not os.path.exists(fname):
        raise FileNotFoundError("FASTA file '{}' not found.".format(fname))
    if not fname.endswith(sfx):
        raise ValueError("FASTA file extension '{}' not registered."
            .format(os.path.splitext(fname)[-1]))
    if os.path.getsize(fname) == 0:
        raise OSError("FASTA file '{}' is empty.".format(fname))
    return fname


def get_outdir():
    """Get output directory.
    :returns: (str) output path
    """
    return config['output']['basedir']


def get_svtype():
    """Get one or more SV types set to non-zero count.
    :returns: (str) SV type(s)
    """
    types = []
    for sv, arr in config['simulation']['sv_type'].items():
        count = arr[0]
        if count > 0:
            types.append(sv)
    types.sort()
    return '_'.join(types)


def get_nthreads(logical=True):
    """Get the number of threads used by `samtools` and `bwa`.
    :returns: (int) threads (default: -1 = number of logical cores)
    """
    n = int(config['threads'])
    if n > 0:
        return n
    else:
        return ps.cpu_count(logical)


def get_mem():
    """Get free memory per core used by `samtools sort`.
    :returns: (int) in MB
    """
    mem = ps.virtual_memory().free / get_nthreads() / 2**20
    return int(mem)
