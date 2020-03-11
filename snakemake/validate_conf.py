from ruamel import yaml
from typing import Dict, List, Union
import yatiml


# Create document class
class Input:
    def __init__(
            self,
            input: Dict[str, str],
            fasta: str,
            seqids: List[Union[int, str]]) -> None:
        self.input = input
        self.fasta = fasta
        self.seqids = seqids

# Class Output:
#   def __init__(self, basedir: data/out  # relative or absolute path
#   genotype:   # diploid genomes
#     - hmz     # homozygous
#     - hmz-sv  # homozygous with SVs
#     - htz-sv  # heterozygous with SVs

# Create loader
class MyLoader(yatiml.Loader):
    pass

yatiml.add_to_loader(MyLoader, Input)
yatiml.set_document_type(MyLoader, Input)

# Load YAML
yaml_text = """
input:
    fasta: data/chr22_44-45Mb.GRCh37.fasta
    seqids: [22]
"""
doc = yaml.load(yaml_text, Loader=MyLoader)
print(type(doc))
print(doc)
print(doc.fasta)
print(doc.seqids)
