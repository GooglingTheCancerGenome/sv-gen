rule art_illumina:
    input:
        fasta = os.path.join(config.output.basedir, '{svtype}',
                             '{genotype}' + config.filext.fasta)
    output:
        fastq1 = os.path.join(config.output.basedir, '{svtype}',
                              'rlen_{read_len}', 'ilen_{insert_len}', 'isd_{insert_sd}',
                              '{genotype}_1' + config.filext.fastq),
        fastq2 = os.path.join(config.output.basedir, '{svtype}',
                              'rlen_{read_len}', 'ilen_{insert_len}', 'isd_{insert_sd}',
                              '{genotype}_2' + config.filext.fastq)
    params:
        seed = config.simulation.seed,
        profile = config.simulation.profile,
        cov = max(config.simulation.coverage),
        prefix = os.path.join(config.output.basedir, '{svtype}', 'rlen_{read_len}',
                              'ilen_{insert_len}', 'isd_{insert_sd}', '{genotype}_')
    conda:
        "../environment.yaml"
    shell:
        """
        set -xe

        art_illumina \
            -M \
            -na \
            -p \
            -f {params.cov} \
            -l {wildcards.read_len} \
            -m {wildcards.insert_len} \
            -s {wildcards.insert_sd} \
            -ss {params.profile} \
            -rs {params.seed} \
            -i "{input.fasta}" \
            -o "{params.prefix}"
        """
