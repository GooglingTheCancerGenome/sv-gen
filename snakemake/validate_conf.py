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
    hmz = 'homozygous_without_sv'
    hmz_sv = 'homozygous_with_sv'
    htz_sv = 'heterozygous_with_sv'


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
            filext: FileExtension) -> None:
        self.threads = threads
        self.input = input
        self.output = output
        self.filext = filext  # file_exts->filext


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
    - hmz_sv
    - htz_sv
filext:
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
print("\n".join([str(g) for g in doc.output.genotype]))
print(doc.filext.fasta)
print(doc.filext.fasta_idx)
print(doc.filext.fastq)
print(doc.filext.bam)
print(doc.filext.bam_idx)
print(doc.filext.bed)
print(doc.filext.vcf)
