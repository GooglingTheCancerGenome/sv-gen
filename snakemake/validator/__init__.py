"""YAML validator for analysis.yaml file."""
import os
import enum
import logging

from typing import Dict, List, Union
from ruamel import yaml
from pyfaidx import Fasta

import yatiml


# Create document classes
class Input:
    """
    Input class to hold the 'input' node.
    """
    def __init__(self, fasta: str, seqids: List[Union[int, str]]) -> None:
        self.fasta = fasta
        self.seqids = seqids


class Genotype(enum.Enum):
    """
    Genotype class to hold the 'genotype' node.
    """
    HOMOZYG_NOSV = 'hmz'
    HOMOZYG_SV = 'hmz-sv'
    HETEROZYG_SV = 'htz-sv'

    @classmethod
    def _yatiml_savorize(cls, node: yatiml.Node) -> None:
        yaml_to_py = {v._value_: v._name_ for v in cls.__members__.values()}
        if node.is_scalar(str):
            node.set_value(yaml_to_py.get(node.get_value()))


class Output:
    """
    Output class to hold the 'output' node.
    """
    def __init__(self, basedir: str, genotype: List[Genotype]) -> None:
        self.basedir = basedir
        self.genotype = genotype


class FileExtension:
    """
    FileExtension class to hold the 'filext' node.
    """
    def __init__(self, fasta: str, fasta_idx: List[str], fastq: str, bam: str,
                 bam_idx: str, bed: str, vcf: str) -> None:
        self.fasta = fasta
        self.fasta_idx = fasta_idx
        self.bam = bam
        self.bam_idx = bam_idx
        self.fastq = fastq
        self.bed = bed
        self.vcf = vcf


class Edit:
    """
    Edit class to hold the attributes of each SV type in the 'svtype' node.
    """
    def __init__(self, count: int, min_len: int, max_len: int) -> None:
        self.count = count
        self.min_len = min_len
        self.max_len = max_len

    @classmethod
    def _yatiml_recognize(cls, node: yatiml.UnknownNode) -> None:
        node.require_sequence()
        err_msg = ('Edit descriptions must be arrays of INTs: count, \
                   min_len and max_len')
        if len(node.yaml_node.value) != 3:
            raise yatiml.RecognitionError(err_msg)
        for item in node.yaml_node.value:
            if (not isinstance(item, yaml.ScalarNode) or
               item.tag != 'tag:yaml.org,2002:int'):
                raise yatiml.RecognitionError(err_msg)

    @classmethod
    def _yatiml_savorize(cls, node: yatiml.Node) -> None:
        if node.is_sequence():
            items = node.seq_items()
            node.make_mapping()
            node.set_attribute('count', items[0].get_value())
            node.set_attribute('min_len', items[1].get_value())
            node.set_attribute('max_len', items[2].get_value())


class SvType:
    """
    SvType class to hold the 'svtype' node.
    """
    def __init__(self, dup: Edit, inv: Edit, tra: Edit, indel: Edit,
                 invdel: Edit, invdup: Edit) -> None:
        self.dup = dup
        self.inv = inv
        self.tra = tra
        self.indel = indel
        self.invdel = invdel
        self.invdup = invdup


class Read:
    """
    Read class to hold DNA read attribute(s).
    """
    def __init__(self, length: List[int]) -> None:
        self.length = length


class Insert:
    """
    Insert class to hold DNA insert attribute(s)
    """
    def __init__(self, stdev: List[int], length: List[int]) -> None:
        self.stdev = stdev
        self.length = length


class Simulation:
    """
    Simulation class to hold the 'simulation' node.
    """
    def __init__(self, config: str, svtype: SvType, seed: int, profile: str,
                 coverage: List[int], read: Read, insert: Insert) -> None:
        self.config = config
        self.svtype = svtype
        self.seed = seed
        self.profile = profile
        self.coverage = coverage
        self.read = read
        self.insert = insert


class Analysis:
    """
    Analysis class to hold all nodes.
    """
    def __init__(self, threads: int, memory: int, tmpspace: int, input: Input,
                 output: Output, filext: FileExtension, simulation: Simulation) \
                 -> None:
        self.threads = threads
        self.memory = memory
        self.tmpspace = tmpspace
        self.input = input
        self.output = output
        self.filext = filext
        self.simulation = simulation
        self._validate()

    def _seqids(self):
        fname = self.input.fasta
        with Fasta(fname) as fasta:
            headers = fasta.keys()
            n = len(headers)
            n_sel = len(self.input.seqids)
            if n_sel == 0:
                return n
            for s in self.input.seqids:
                if str(s) not in headers:
                    raise ValueError("SeqID '{}' is not in the FASTA file '{}'."
                        .format(s, fname))
            return n_sel

    def _filext(self):
        fname = self.input.fasta
        fext = self.filext.fasta
        if not fname.endswith(fext):
            raise ValueError("FASTA file extension '{}' is not registered."
                .format(os.path.splitext(fname)[-1]))

    def _simulation(self):
        count = sum([sv.count for sv in vars(self.simulation.svtype).values()])
        if count < 1:
            raise ValueError("Select at least one SV type by setting its count to non-zero value.")
        if self.simulation.svtype.tra.count > 0 and self._seqids() == 1:
            raise ValueError("Two or more chromosomes are required to simulate translocations.")

    def _validate(self):
        self._seqids()
        self._filext()
        self._simulation()


# Create loader
class MyLoader(yatiml.Loader):
    """
    MyLoader class.
    """
    pass


def load_configfile(yaml_file: str) -> Dict:
    with open(yaml_file, 'r') as conf:
        return yaml.load(conf, MyLoader)


yatiml.logger.setLevel(logging.DEBUG)
yatiml.add_to_loader(MyLoader, [
    Input, Genotype, Output, FileExtension, Edit, SvType, Read, Insert,
    Simulation, Analysis
])
yatiml.set_document_type(MyLoader, Analysis)
