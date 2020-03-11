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


class Output:
    def __init__(
            self,
            basedir: str,
            genotype: List[str]) -> None:
        self.basedir = basedir
        self.genotype = genotype


class Analysis:
    def __init__(self, input: Input, output: Output) -> None:
        self.input = input
        self.output = output

# Create loader
class MyLoader(yatiml.Loader):
    pass

yatiml.add_to_loader(MyLoader, [Input, Genotype, Output, Analysis])
yatiml.set_document_type(MyLoader, Analysis)

# Load YAML
yaml_text = """
input:
    fasta: data/chr22_44-45Mb.GRCh37.fasta
    seqids: [22]
output:
  basedir: data/out
  genotype:
    - hmz
    - hmz-sv
    - htz-sv
"""
doc = yaml.load(yaml_text, Loader=MyLoader)
print(type(doc))
print(doc)
print(doc.input.fasta)
print(doc.input.seqids)
print(doc.output.basedir)
print(doc.output.genotype)
