import enum
import logging

from typing import Dict, List, Union
from ruamel import yaml

import yatiml


# Create document classes
class Input:
    def __init__(self, fasta: str, seqids: List[Union[int, str]]) -> None:
        self.fasta = fasta
        self.seqids = seqids


class Genotype(enum.Enum):
    HOMOZYG_NOSV = 'hmz'
    HOMOZYG_SV = 'hmz-sv'
    HETEROZYG_SV = 'htz-sv'

    @classmethod
    def _yatiml_savorize(cls, node: yatiml.Node) -> None:
        yaml_to_py = {v._value_: v._name_ for v in cls.__members__.values()}
        if node.is_scalar(str):
            node.set_value(yaml_to_py.get(node.get_value()))


class Output:
    def __init__(self, basedir: str, genotype: List[Genotype]) -> None:
        self.basedir = basedir
        self.genotype = genotype


class FileExtension:
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
            if (not isinstance(item, yaml.ScalarNode) or \
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
    def __init__(self, dup: Edit, inv: Edit, tra: Edit, indel: Edit,
                 invdel: Edit, invdup: Edit) -> None:
        self.dup = dup
        self.inv = inv
        self.tra = tra
        self.indel = indel
        self.invdel = invdel
        self.invdup = invdup


class Read:
    def __init__(self, length: List[int]) -> None:
        self.length = length


class Insert:
    def __init__(self, stdev: int, length: List[int]) -> None:
        self.stdev = stdev
        self.length = length


class Simulation:
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
    def __init__(self, threads: int, input: Input, output: Output,
                 filext: FileExtension, simulation: Simulation) -> None:
        self.threads = threads
        self.input = input
        self.output = output
        self.filext = filext
        self.simulation = simulation


# Create loader
class MyLoader(yatiml.Loader):
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
