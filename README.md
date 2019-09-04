# sv-gen-workflow

Snakemake-based workflow to generate artificial structural variant (SV) data.

### Dependencies

-   python (>=3.6)
-   [conda](https://conda.io/) (>=4.5)
-   [snakemake](https://snakemake.readthedocs.io/) (>=4.8)

The workflow installs the following tools:

-   [SURVIVOR](https://github.com/fritzsedlazeck/SURVIVOR) (1.0.5)
-   [ART](https://www.niehs.nih.gov/research/resources/software/biostatistics/art/) (2016-06-05)
-   [BWA](https://github.com/lh3/bwa) (0.7.17)
-   [Sambamba](https://github.com/biod/sambamba) (0.7.0) or
-   [Samtools](https://github.com/samtools/samtools) (1.9)
