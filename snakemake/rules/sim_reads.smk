rule art_illumina:
    input:
        fasta = os.path.join(config.output.basedir, '{svtype}',
                             '{genotype}' + config.filext.fasta)
    output:
        fastq1 = os.path.join(config.output.basedir, '{svtype}',
                              'r{read_len}_i{insert_len}',
                              '{genotype}_1' + config.filext.fastq),
        fastq2 = os.path.join(config.output.basedir, '{svtype}',
                              'r{read_len}_i{insert_len}',
                              '{genotype}_2' + config.filext.fastq)
    params:
        seed = config.simulation.seed,
        profile = config.simulation.profile,
        insert_stdev = config.simulation.insert.stdev,
        coverage = max(config.simulation.coverage),
        prefix = os.path.join(config.output.basedir, '{svtype}',
                              'r{read_len}_i{insert_len}', '{genotype}_')
    conda:
        "../environment.yaml"
    shell:
        """
        set -xe

        art_illumina \
            -M \
            -na \
            -p \
            -f {params.coverage} \
            -l {wildcards.read_len} \
            -m {wildcards.insert_len} \
            -s {params.insert_stdev} \
            -ss {params.profile} \
            -rs {params.seed} \
            -i "{input.fasta}" \
            -o "{params.prefix}"
        """
