from ruamel import yaml
from typing import Dict, List, Union
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
    HOMO_NOSV = 'hmz'
    HOMO_SV = 'hmz-sv'
    HETERO_SV = 'htz-sv'

    @classmethod
    def yatiml_savorize(cls, node: yatiml.Node) -> None:
        yaml_to_python = {
                'hmz': 'HOMO_NOSV',
                'hmz-sv': 'HOMO_SV',
                'htz-sv': 'HETERO_SV'}

        if node.is_scalar(str):
            yaml_str = node.get_value()
            if yaml_str in yaml_to_python:
                node.set_value(yaml_to_python[yaml_str])
            else:
                raise yatiml.SeasoningError('Invalid Genotype value')

    @classmethod
    def yatiml_sweeten(cls, node: yatiml.Node) -> None:
        python_to_yaml = {
                'HOMO_NOSV': 'hmz',
                'HOMO_SV': 'hmz-sv',
                'HETERO_SV': 'htz-sv'}
        python_str = node.get_value()
        if python_str in python_to_yaml:
            node.set_value(python_to_yaml[python_str])
        else:
            raise yatiml.SeasoningError('Invalid Genotype value')


class Output:
    def __init__(
            self,
            basedir: str,
            genotype: Genotype) -> None:
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
print(type(doc))
print(doc)
print(doc.threads)
print(doc.input.fasta)
print(doc.input.seqids)
print(doc.output.basedir)
print(doc.output.genotype)
print(doc.filext.fasta)
print(doc.filext.fasta_idx)
print(doc.filext.fastq)
print(doc.filext.bam)
print(doc.filext.bam_idx)
print(doc.filext.bed)
print(doc.filext.vcf)
