rule art_illumina:
    input:
        fasta = os.path.join(get_outdir(), '{svtype}',
                             '{genotype}' + get_filext('fasta'))
    output:
        fastq1 = os.path.join(get_outdir(), '{svtype}',
                              'r{read_len}_i{insert_len}',
                              '{genotype}_1' + get_filext('fastq')),
        fastq2 = os.path.join(get_outdir(), '{svtype}',
                              'r{read_len}_i{insert_len}',
                              '{genotype}_2' + get_filext('fastq'))
    params:
        seed = config['sim_reads']['seed'],
        profile = config['sim_reads']['profile'],
        insert_stdev = config['sim_reads']['insert']['stdev'],
        coverage = max(config['sim_reads']['coverage']),
        prefix = os.path.join(get_outdir(), '{svtype}',
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
