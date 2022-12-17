"""YAML validator for analysis.yaml file."""
import enum
import logging
from typing import Dict, List, Optional, Union

import yatiml
from ruamel import yaml


# Create document classes
class Input:
    """
    Input class to hold the 'input' node.
    """

    def __init__(self, fasta: str, seqids: Optional[List[Union[int, str]]] = None) -> None:
        self.fasta = fasta
        self.seqids = seqids

    @classmethod
    def _yatiml_recognize(cls, node: yatiml.UnknownNode) -> None:
        node.require_attribute('fasta', str)


class Genotype(enum.Enum):
    """
    Genotype class to hold the 'genotype' node.
    """
    HOMOZYG_NOSV = 'hmz'
    HOMOZYG_SV = 'hmz-sv'
    HETEROZYG_SV = 'htz-sv'

    @classmethod
    def _yatiml_recognize(cls, node: yatiml.UnknownNode) -> None:
        node.require_scalar(str)

    @classmethod
    def _yatiml_savorize(cls, node: yatiml.Node) -> None:
        yaml_to_py = {v.value: v.name for v in cls.__members__.values()}
        if node.is_scalar(str):
            val = node.yaml_node.value
            if val not in yaml_to_py:
                raise yatiml.SeasoningError(
                    "Invalid genotype: '{}'".format(val))
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
    Edit class holds the attributes for each 'svtype' node.
    """

    def __init__(self, count: int, min_len: int, max_len: int) -> None:
        self.count = count
        self.min_len = min_len
        self.max_len = max_len

    @classmethod
    def _yatiml_recognize(cls, node: yatiml.UnknownNode) -> None:
        node.require_sequence()
        err_msg = (
            'svtype attributes must be an array of INTs: [count, min_len, max_len]')
        tag = 'tag:yaml.org,2002:int'
        if len(node.yaml_node.value) != 3:
            raise yatiml.RecognitionError(err_msg)
        for item in node.yaml_node.value:
            if (not isinstance(item, yaml.ScalarNode) or item.tag != tag):
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

    @classmethod
    def _yatiml_recognize(cls, node: yatiml.UnknownNode) -> None:
        node.require_mapping()
        svtype_count = 0
        for item in node.yaml_node.value:
            svtype_count += int(item[1].value[0].value)
        if svtype_count == 0:
            raise yatiml.RecognitionError(
                "At least one SV type must have non-zero count.")


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


# Create loader
load = yatiml.load_function(Analysis, Simulation, Insert, Read, SvType, Edit,
                            FileExtension, Output, Genotype, Input)
yatiml.logger.setLevel(logging.DEBUG)


def load_configfile(yaml_file: str) -> Dict:
    """
    Redefined Snakemake function.
    """
    with open(yaml_file, 'r') as conf:
        return load(conf)
