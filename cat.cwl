cwlVersion: "cwl:v1.0"

class: CommandLineTool
baseCommand: cat
stdout: merge.fasta
inputs:
  file_a:
    type: File
    inputBinding:
      position: 1
  file_b:
    type: File
    inputBinding:
      position: 2



outputs:
  - id: example_out
    type: File
    outputBinding:
      glob: merge.fasta
