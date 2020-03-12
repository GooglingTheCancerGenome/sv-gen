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



class SVType(enum.Enum):
    DUP = 'DUP'
    INV = 'INV'
    TRA = 'TRA'
    INDEL = 'INDEL'
    INV_DEL = 'INV-DEL'
    INV_DUP = 'INV-DUP'

    @classmethod
    def yatiml_savorize(cls, node: yatiml.Node) -> None:
        yaml_to_python = {
                v._value_: v._name_ for v in cls.__members__.values()}
        if node.is_scalar(str):
            node.set_value(yaml_to_python.get(node.get_value()))


class SimGenome:
    def __init__(
            self,
            config: str,
            sv_type: Dict[str, List[int]]) -> None:  # keys: restrict to SVTypes
                                                     # values: restrict to Tuple[int, int, int]
        self.config = config
        self.sv_type = sv_type


class SimRead:
    def __init__(
            self,
            seed: int,
            profile: str,
            coverage: List[int],
            read: Dict[str, List[int]],
            insert: Dict[str, Union[int, List[int]]]) -> None:
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
            file_exts: FileExtension,
            sim_genomes: SimGenome,
            sim_reads: SimRead) -> None:
        self.threads = threads
        self.input = input
        self.output = output
        self.file_exts = file_exts
        self.sim_genomes = sim_genomes
        self.sim_reads = sim_reads


# Create loader
class MyLoader(yatiml.Loader):
    pass


yatiml.logger.setLevel(logging.DEBUG)
yatiml.add_to_loader(MyLoader, [Input, Genotype, Output, FileExtension,
                     SVType, SimGenome, SimRead, Analysis])
yatiml.set_document_type(MyLoader, Analysis)

# Load YAML
yaml_text = """
threads: -1
input:
  fasta: data/chr22_44-45Mb.GRCh37.fasta
  seqids: [22]
output:
  basedir: data/out
  genotype:
    - hmz
    - hmz-sv
    - htz-sv
file_exts:
  fasta: .fasta
  fasta_idx:
    - .fasta.ann
    - .fasta.amb
    - .fasta.bwt
    - .fasta.pac
    - .fasta.sa
  fastq: .fq
  bam: .bam
  bam_idx: .bam.bai
  bed: .bed
  vcf: .vcf
sim_genomes:
  config: survivor.cfg
  sv_type:
    DUP: [0, 100, 10000]
    INV: [0, 600, 800]
    TRA: [0, 1000, 3000]
    INDEL: [10, 20, 500]
    INV-DEL: [0, 600, 800]
    INV-DUP: [0, 600, 800]
sim_reads:
  seed: 1000
  profile: HSXt
  coverage: [10, 30]
  read:
    len: [150]
  insert:
    stdev: 10
    len: [500]
"""
#with open('analysis.yaml', 'r') as conf:
doc = yaml.load(yaml_text, Loader=MyLoader)
print(doc)
print(doc.threads)
print(doc.input.fasta)
print(doc.input.seqids)
print(doc.output.basedir)
print([str(g.value) for g in doc.output.genotype])
print(doc.file_exts.fasta)
print(doc.file_exts.fasta_idx)
print(doc.file_exts.fastq)
print(doc.file_exts.bam)
print(doc.file_exts.bam_idx)
print(doc.file_exts.bed)
print(doc.file_exts.vcf)
print(doc.sim_genomes.config)
print(doc.sim_genomes.sv_type)
print(doc.sim_reads.seed)
print(doc.sim_reads.profile)
print(doc.sim_reads.coverage)
print(doc.sim_reads.read)
print(doc.sim_reads.insert)
