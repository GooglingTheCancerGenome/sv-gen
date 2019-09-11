# sv-gen-workflow

Snakemake-based workflow to generate artificial structural variant (SV) data.

## Dependencies

-   python (>=3.6)
-   [conda](https://conda.io/) (>=4.5)
-   [snakemake](https://snakemake.readthedocs.io/) (>=4.8)

The workflow installs the following tools:

-   [SURVIVOR](https://github.com/fritzsedlazeck/SURVIVOR) (1.0.6)
-   [ART](https://www.niehs.nih.gov/research/resources/software/biostatistics/art/) (2016-06-05)
-   [BWA](https://github.com/lh3/bwa) (0.7.17)
-   [Sambamba](https://github.com/biod/sambamba) (0.7.0) or
-   [Samtools](https://github.com/samtools/samtools) (1.9)

**1. Clone this repo.**

```bash
git clone https://github.com/GooglingTheCancerGenome/sv-gen-workflow.git
cd sv-gen-workflow/snakemake
```

**2. Install dependencies.**

```bash
wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh  # download Miniconda installer (with Python 3)
bash miniconda.sh  # install Conda (accept defaults)
# add the following line to ~/.bashrc & source env
# export PATH="$HOME/miniconda3/bin:$PATH"
source ~/.bashrc
conda update -y conda  # update Conda
conda create -y -n wf2 && source activate wf2  # create & activate new env
conda install -y -c bioconda snakemake
```

**3. Configure the workflow.**

-   **config files**:
    -   `analysis.yaml` - analysis-specific settings
    -   `environment.yaml` - software dependencies and versions

**4. Execute the workflow.**

```bash
snakemake -np  # 'dry' run only checks I/O files
snakemake --use-conda
```
