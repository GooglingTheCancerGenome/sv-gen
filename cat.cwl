cwlVersion: "cwl:v1.0"

requirements:
- class: InlineJavascriptRequirement


class: CommandLineTool
baseCommand: cat
inputs:
  file_a:
    type: File
    inputBinding:
      position: 2
  file_b:
    type: File
    inputBinding:
      position: 3
outputs:
  - id: con_fasta
    type: File
    outputBinding:
      glob:  $('merged_'+inputs.file_b.basename)

stdout: $('merged_'+inputs.file_b.basename)