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


class Analysis:
    def __init__(
            self,
            threads: int,
            input: Input,
            output: Output,
            file_exts: FileExtension) -> None:
        self.threads = threads
        self.input = input
        self.output = output
        self.file_exts = file_exts


# Create loader
class MyLoader(yatiml.Loader):
    pass


yatiml.logger.setLevel(logging.DEBUG)
yatiml.add_to_loader(MyLoader, [Input, Genotype, Output, FileExtension, Analysis])
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
"""
doc = yaml.load(yaml_text, Loader=MyLoader)
print(doc)
print(doc.threads)
print(doc.input.fasta)
print(doc.input.seqids)
print(doc.output.basedir)
print("\n".join([str(g.value) for g in doc.output.genotype]))
print(doc.file_exts.fasta)
print(doc.file_exts.fasta_idx)
print(doc.file_exts.fastq)
print(doc.file_exts.bam)
print(doc.file_exts.bam_idx)
print(doc.file_exts.bed)
print(doc.file_exts.vcf)
