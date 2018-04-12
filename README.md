# sv-gen-workflow
CWL workflows and helper script to generate artificial structural variant (SV) data used to train the CNN model.

* chr17 of the hg19 Human genome assembly is used as the test chromosome
* SURVIVOR is used to introduce SVs only on one of the two homologous chromosomes
* The ART read simulator is used to generate paired-end reads
* bwa-mem is used to map the reads on the hg19 assembly

Output of the workflow are two categories of Tumor/Normal pairs of BAM alignment files: Somatic and Germline.

* The Somatic category contains a Tumor BAM file with 50% tumor purity and a Normal BAM file which does not contain SVs. The model is trained to recognize SVs that are present in the Tumor but not in the Normal BAM files as 'somatic'.
* The Germline category contains two BAM files with the same set of SVs, so the model is trained to recognize SVs that are present in both the Tumor and Normal BAM files as 'germline'.


# Requirements SURVIVOR pipeline
Ensure that the following tools are in your $PATH:
- SURVIVOR
- Samtools

Create a python virtualenv and install the following packages:
- argparse2tool
- cwltool 

Example:

```
virtualenv -p python2 venv_cwl
source /venv_cwl/bin/activate
pip install cwltool argparse argparse2tool

````
Run Survivor pipeline

````
PATH=$PATH:/path/to/sv-gen-workflow/tools cwltool --preserve-environment PATH run_SURVIVOR.cwl json/run_SURVIVOR.json

```

