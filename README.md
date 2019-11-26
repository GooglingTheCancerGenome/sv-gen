# sv-gen

[![Build Status](https://travis-ci.org/GooglingTheCancerGenome/sv-gen.svg?branch=master)](https://travis-ci.org/GooglingTheCancerGenome/sv-gen)

Snakemake-based workflow to generate artificial structural variant (SV) data.

## Dependencies

-   python (>=3.6)
-   [conda](https://conda.io/) (>=4.5)
-   [snakemake](https://snakemake.readthedocs.io/) (>=4.8)
-   [xenon-cli](https://github.com/NLeSC/xenon-cli) (3.0.4)

The workflow includes the following tools:

-   [SURVIVOR](https://github.com/fritzsedlazeck/SURVIVOR) (1.0.6)
-   [ART](https://www.niehs.nih.gov/research/resources/software/biostatistics/art/) (2016-06-05)
-   [BWA](https://github.com/lh3/bwa) (0.7.17)
-   [Samtools](https://github.com/samtools/samtools) (1.9)

**1. Clone this repo.**

```bash
git clone https://github.com/GooglingTheCancerGenome/sv-gen.git
cd sv-gen
```

**2. Install dependencies.**

```bash
wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh  # download Miniconda installer (with Python 3)
bash miniconda.sh  # install Conda (accept defaults)
# add the following line to ~/.bashrc & source env
# export PATH="$HOME/miniconda3/bin:$PATH"
source ~/.bashrc
conda update -y conda  # update Conda
conda env create -n wf -f environment.yaml
conda activate wf
cd snakemake
```

**3. Configure the workflow.**

-   **config files**:
    -   `analysis.yaml` - analysis-specific settings
    -   `environment.yaml` - software dependencies and versions

**4. Execute the workflow.**

```bash
snakemake -np  # 'dry' run only checks I/O files
snakemake --use-conda  # run simulations locally
```

_Submit jobs to Grid Engine-based cluster_

```bash
snakemake --use-conda --latency-wait 30 --jobs 17 \
--cluster 'xenon scheduler gridengine --location local:// submit --name smk.{rule} --inherit-env --cores-per-task 1 --max-run-time 5 --working-directory . --stderr stderr-%j.log --stdout stdout-%j.log' &>smk.log&
```

_Submit jobs to Slurm-based cluster_

```bash
snakemake --use-conda --latency-wait 30 --jobs 17 \
--cluster 'xenon scheduler slurm --location local:// submit --name smk.{rule} --inherit-env --cores-per-task 1 --max-run-time 5 --working-directory . --stderr stderr-%j.log --stdout stdout-%j.log' &>smk.log&
```
