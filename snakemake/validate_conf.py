from ruamel import yaml
from typing import Dict, List, Union

import logging
import yatiml
import enum


# Create document classes
class Input:
    def __init__(
            self,
            fasta: str,
            seqids: List[Union[int, str]]) -> None:
        self.fasta = fasta
        self.seqids = seqids


class Genotype(enum.Enum):
    HOMOZYG_NOSV = 'hmz'
    HOMOZYG_SV = 'hmz-sv'
    HETEROZYG_SV = 'htz-sv'

    @classmethod
    def yatiml_savorize(cls, node: yatiml.Node) -> None:
        yaml_to_python = {
                v._value_: v._name_ for v in cls.__members__.values()}
        if node.is_scalar(str):
            node.set_value(yaml_to_python.get(node.get_value()))


class Output:
    def __init__(
            self,
            basedir: str,
            genotype: List[Genotype]) -> None:
        self.basedir = basedir
        self.genotype = genotype


class FileExtension:
    def __init__(
            self,
            fasta: str,
            fasta_idx: List[str],
            fastq: str,
            bam: str,
            bam_idx: str,
            bed: str,
            vcf: str) -> None:
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
    def yatiml_recognize(cls, node: yatiml.UnknownNode) -> None:
        node.require_sequence()
        err_msg = ('Edit descriptions must be arrays of INTs: count, \
                   min_len and max_len')
        if len(node.yaml_node.value) != 3:
            raise yatiml.RecognitionError(err_msg);
        for item in node.yaml_node.value:
            if (not isinstance(item, yaml.ScalarNode) or
               item.tag != 'tag:yaml.org,2002:int'):
                raise yatiml.RecognitionError(err_msg)

    @classmethod
    def yatiml_savorize(cls, node: yatiml.Node) -> None:
        if node.is_sequence():
            items = node.seq_items()
            node.make_mapping()
            node.set_attribute('count', items[0].get_value())
            node.set_attribute('min_len', items[1].get_value())
            node.set_attribute('max_len', items[2].get_value())


class SvType:
    def __init__(
            self,
            dup: Edit,
            inv: Edit,
            tra: Edit,
            indel: Edit,
            invdel: Edit,
            invdup: Edit) -> None:
        self.dup = dup
        self.inv = inv
        self.tra = tra
        self.indel = indel
        self.invdel = invdel
        self.invdup = invdup


class Simulation:
    def __init__(
            self,
            config: str,
            sv_type: SvType,
            seed: int,
            profile: str,
            coverage: List[int],
            read: Dict[str, List[int]],
            insert: Dict[str, Union[int, List[int]]]) -> None:
        self.config = config
        self.sv_type = sv_type
        self.seed = seed
        self.profile = profile
        self.coverage = coverage
        self.read = read
        self.insert = insert


class Analysis:
    def __init__(
            self,
            threads: int,
            input: Input,
            output: Output,
            filext: FileExtension,
            simulation: Simulation) -> None:
        self.threads = threads
        self.input = input
        self.output = output
        self.filext = filext
        self.simulation = simulation


# Create loader
class MyLoader(yatiml.Loader):
    pass


yatiml.logger.setLevel(logging.DEBUG)
yatiml.add_to_loader(MyLoader, [Input, Genotype, Output, FileExtension, Edit,
                     SvType, Simulation, Analysis])
yatiml.set_document_type(MyLoader, Analysis)

with open('analysis.yaml', 'r') as conf:
    doc = yaml.load(conf, Loader=MyLoader)
    print(doc)
    print(doc.threads)
    print(doc.input.fasta)
    print(doc.input.seqids)
    print(doc.output.basedir)
    print([str(g.value) for g in doc.output.genotype])
    print(doc.filext.fasta)
    print(doc.filext.fasta_idx)
    print(doc.filext.fastq)
    print(doc.filext.bam)
    print(doc.filext.bam_idx)
    print(doc.filext.bed)
    print(doc.filext.vcf)
    print(doc.simulation.config)
    print(doc.simulation.sv_type.dup.count,
          doc.simulation.sv_type.dup.min_len,
          doc.simulation.sv_type.dup.max_len)
    print(doc.simulation.sv_type.inv.count,
          doc.simulation.sv_type.inv.min_len,
          doc.simulation.sv_type.inv.max_len)
    print(doc.simulation.sv_type.tra.count,
          doc.simulation.sv_type.tra.min_len,
          doc.simulation.sv_type.tra.max_len)
    print(doc.simulation.sv_type.indel.count,
          doc.simulation.sv_type.indel.min_len,
          doc.simulation.sv_type.indel.max_len)
    print(doc.simulation.sv_type.invdel.count,
          doc.simulation.sv_type.invdel.min_len,
          doc.simulation.sv_type.invdel.max_len)
    print(doc.simulation.sv_type.invdup.count,
          doc.simulation.sv_type.invdup.min_len,
          doc.simulation.sv_type.invdup.max_len)
    print(doc.simulation.seed)
    print(doc.simulation.profile)
    print(doc.simulation.coverage)
    print(doc.simulation.read)
    print(doc.simulation.insert)
