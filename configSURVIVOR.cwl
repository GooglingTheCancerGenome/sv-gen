#!/usr/bin/env cwl-runner
# This tool description was generated automatically by argparse2tool ver. 0.4.3-2
# To generate again: $ configSURVIVOR.py --generate_cwl_tool
# Help: $ configSURVIVOR.py --help_arg2cwl

cwlVersion: "cwl:v1.0"

class: CommandLineTool
baseCommand: ['configSURVIVOR.py']

doc: |
  Provide SURVIVOR runtime parameters.

inputs:
  RUNMODE:
    type: string
    default: inv
    inputBinding:
      prefix: --RUNMODE

  SURVPARAM:
    type: string
    default: 5000/500/5000
    inputBinding:
      prefix: --SURVPARAM
outputs:
  - id: example_out
    type: File
    outputBinding:
      glob: sv_template.txt

