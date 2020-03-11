from ruamel import yaml
from typing import Dict, List, Union
import yatiml


# Create document class
class Input:
    def __init__(
            self,
            fasta: str,
            seqids: List[Union[int, str]]) -> None:
        self.fasta = fasta
        self.seqids = seqids

class Analysis:
    def __init__(self, input: Input) -> None:
        self.input = input

# Create loader
class MyLoader(yatiml.Loader):
    pass

yatiml.add_to_loader(MyLoader, [Input, Analysis])
yatiml.set_document_type(MyLoader, Analysis)

# Load YAML
yaml_text = """
input:
    fasta: data/chr22_44-45Mb.GRCh37.fasta
    seqids: [22]
"""
doc = yaml.load(yaml_text, Loader=MyLoader)
print(type(doc))
print(doc)
print(doc.input.fasta)
print(doc.input.seqids)
