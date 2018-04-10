cwlVersion: "cwl:v1.0"

class: CommandLineTool

requirements:
- class: InlineJavascriptRequirement


inputs:
  Reference_file:
    type: File
    inputBinding:
      position: 1

  Parameter_file:
    type: File
    inputBinding:
      position: 2

  snp_mutation_freq:
    type: int
    inputBinding:
      position: 3

  read_type:
    type: int
    inputBinding:
      position: 4

  output_prefix:
    type: string
    inputBinding:
      position: 5


outputs:
   - id: fasta_out
     type: File
     outputBinding:
       glob: $(inputs.output_prefix).fasta
   - id: bed_out
     type: File
     outputBinding:
       glob: $(inputs.output_prefix).bed
   - id: vcf_out
     type: File
     outputBinding:
       glob: $(inputs.output_prefix).vcf



baseCommand: [SURVIVOR, simSV]
