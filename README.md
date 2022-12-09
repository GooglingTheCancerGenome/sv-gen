# sv-gen

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3725664.svg)](https://doi.org/10.5281/zenodo.3725664)
[![CI](https://github.com/GooglingTheCancerGenome/sv-gen/actions/workflows/ci.yaml/badge.svg?branch=dev)](https://github.com/GooglingTheCancerGenome/sv-gen/actions/workflows/ci.yaml)
[![Codacy Badge](https://app.codacy.com/project/badge/Grade/7d9a698a93fa44ec8ad79b96842d48ee)](https://www.codacy.com/gh/GooglingTheCancerGenome/sv-gen/dashboard?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=GooglingTheCancerGenome/sv-gen&amp;utm_campaign=Badge_Grade)
[![Codacy Badge](https://app.codacy.com/project/badge/Coverage/7d9a698a93fa44ec8ad79b96842d48ee)](https://www.codacy.com/gh/GooglingTheCancerGenome/sv-gen/dashboard?utm_source=github.com&utm_medium=referral&utm_content=GooglingTheCancerGenome/sv-gen&utm_campaign=Badge_Coverage)

Structural variants (SVs) are an important class of genetic variation implicated in a wide array of genetic diseases. _sv-gen_ is a Snakemake-based workflow to generate artificial short-read alignments based on a reference genome with(out) SVs. The workflow is easy to use and deploy on any Linux-based machine. In particular, the workflow supports automated software deployment, easy configuration and addition of new analysis tools as well as enables to scale from a single computer to different HPC clusters with minimal effort.

## Dependencies

-   [Python 3](https://www.python.org/)
-   [Conda](https://conda.io/) - package/environment management system
-   [Snakemake](https://snakemake.readthedocs.io/) - workflow management system
-   [Xenon CLI](https://github.com/NLeSC/xenon-cli) - command-line interface to compute and storage resources
-   [jq](https://stedolan.github.io/jq/) - command-line JSON processor (optional)
-   [YAtiML](https://github.com/yatiml/yatiml) - library for YAML type inference and schema validation

The workflow ([DAG](/doc/sv-gen.svg)) includes the following tools:

-   [SURVIVOR](https://github.com/fritzsedlazeck/SURVIVOR)
-   [ART](https://www.niehs.nih.gov/research/resources/software/biostatistics/art/)
-   [BWA](https://github.com/lh3/bwa)
-   [Samtools](https://github.com/samtools/samtools)

The software dependencies and versions can be found in the conda `environment.yaml` files ([1](/environment.yaml), [2](/workflow/environment.yaml)).

**1. Clone this repo.**

```bash
git clone https://github.com/GooglingTheCancerGenome/sv-gen.git
cd sv-gen
```

**2. Install dependencies.**

```bash
# download Miniconda3 installer
wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh
# install Conda (respond by 'yes')
bash miniconda.sh
# update Conda
conda update -y conda
# install Mamba
conda install -n base -c conda-forge -y mamba
# create a new environment with dependencies & activate it
mamba env create -n wf -f environment.yaml
conda activate wf
```

**3. Configure the workflow.**

-   **config files**:
    -   [`analysis.yaml`](/config/analysis.yaml) - analysis-specific settings
    -   [`environment.yaml`](/workflow/environment.yaml) - software dependencies and versions

**4. Execute the workflow.**

```bash
cd workflow
# 'dry' run only checks I/O files
snakemake -np

# run the workflow locally
snakemake --use-conda --cores
```

_Submit jobs to Slurm/GridEngine-based cluster_

```bash
SCH=slurm   # or gridengine
snakemake --use-conda --latency-wait 30 --jobs \
--cluster "xenon scheduler $SCH --location local:// submit --name smk.{rule} --inherit-env --max-run-time 5 --working-directory . --stderr stderr-%j.log --stdout stdout-%j.log" &>smk.log&
```

_Query job accounting information_

```bash
SCH=slurm   # or gridengine
xenon --json scheduler $SCH --location local:// list --identifier [jobID] | jq ...
``` 
