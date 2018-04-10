cwlVersion: "cwl:v1.0"
class: Workflow
inputs:
  RUNMODE:
    type: string

  SURVPARAM:
    type: string

  Reference_file:
    type: File

  Parameter_file:
    type: File

  snp_mutation_freq:
    type: int

  read_type:
    type: int

  output_prefix:
    type: string

outputs:
  fasta:
    type: File
    outputSource: runsv/fasta_out

steps:
  configsv:
    run: configSURVIVOR.cwl
    in:
       RUNMODE: RUNMODE
       SURVPARAM: SURVPARAM
    out:
       [example_out]
  runsv:
    run: SURVIVOR_simSV.cwl
    in:
       Reference_file: Reference_file
       Parameter_file: configsv/example_out
       snp_mutation_freq: snp_mutation_freq
       read_type: read_type
       output_prefix: output_prefix
    out: [fasta_out]

