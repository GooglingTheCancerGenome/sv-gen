import os
import sys

from snakemake import load_configfile
from csv import DictReader

config = load_configfile('analysis.yaml')


def get_filext(fmt):
    """Get file extension(s) given file type/format:
        ['fasta', 'fasta_idx', 'bam', 'bam_idx', 'fastq']
    :param fmt: (str) input file format
    :returns: (str) file extension
    """
    if fmt not in config['file_exts'].keys():
         raise ValueError("Input file format '{}' not supported.".format(fmt.lower()))
    return config['file_exts'][fmt]


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
    :returns: (str) outdir
    """
    if 'basedir' not in config['output']:
        raise OSError('Output basedir is not set.')
    return config['output']['basedir']


def get_svtype():
    types = []
    for sv, arr in config['sim_genomes']['sv_type'].items():
        count = arr[0]
        if count > 0:
            types.append(sv)
    types.sort()
    return '_'.join(types)


def get_nthreads():
    """Get the max. number of threads used by the tools (i.e. samtools and bwa).
    :returns: (int) number of threads (default: # cores available)
    """
    if 'threads' not in config:
        raise KeyError("Missing key 'threads' from config.")

    n = config['threads']
    try:
        int(n)
    except ValueError:
        sys.exit('Use INT to set the number of threads.')
    if int(n) > 0:
        return config['threads']
    else:
        return os.cpu_count()
