cwlVersion: "cwl:v1.0"
class: Workflow
inputs:
  file_a:
     type: File
  file_b:
     type: File
  RUNMODE:
    type: string

  SURVPARAM:
    type: string

  snp_mutation_freq:
    type: int

  read_type:
    type: int

  output_prefix:
    type: string

outputs:
  fasta1:
    type: File
    outputSource: catfasta/con_fasta
  fasta2:
    type: File
    outputSource: catfasta2/con_fasta
  faidx_out:
    type: File
    outputSource: faidx/index
  faidx2_out:
    type: File
    outputSource: faidx2/index


steps:
  catfasta:
    run: cat.cwl
    in:
       file_a: file_a
       file_b: file_b
    out:
       [con_fasta]
  faidx:
    run: samtools-faidx.cwl
    in:
       input: catfasta/con_fasta
    out:
       [index]
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
       Reference_file: catfasta/con_fasta
       Parameter_file: configsv/example_out
       snp_mutation_freq: snp_mutation_freq
       read_type: read_type
       output_prefix: output_prefix
    out: [fasta_out]

  catfasta2:
    run: cat.cwl
    in:
       file_a: file_a
       file_b: runsv/fasta_out
    out:
       [con_fasta]
  faidx2:
    run: samtools-faidx.cwl
    in:
       input: catfasta2/con_fasta
    out:
       [index]



