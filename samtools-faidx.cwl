#!/usr/bin/env cwl-runner
#
# To use it as stand alone tool. The working directory should not have input .fa file
#    example: "./samtools-faidx.cwl --input=./test-files/mm10.fa"
# As part of a workflow should be no problem at all

cwlVersion: v1.0
class: CommandLineTool

requirements:
- class: InlineJavascriptRequirement
- class: InitialWorkDirRequirement
  listing:
  - entry: $(inputs.input)
    entryname: $(inputs.input.path.split('/').slice(-1)[0])
inputs:
  input:
    type: File
    doc: <file.fa|file.fa.gz>
  region:
    type: string?
    inputBinding:
      position: 2

outputs:
  index:
    type: File
    outputBinding:
      glob: $(inputs.input.path.split('/').slice(-1)[0]) #+'.fai')
    secondaryFiles:
    - .fai
    - .gzi

baseCommand:
- samtools
- faidx

arguments:
- valueFrom: $(inputs.input.path.split('/').slice(-1)[0])
  position: 1



