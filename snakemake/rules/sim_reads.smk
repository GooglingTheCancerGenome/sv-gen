rule art_illumina:
    input:
        os.path.join("{basedir}", "{genotype}.fasta")
    output:
        fastq1 = os.path.join("{basedir}", "{genotype}_1.fq"),
        fastq2 = os.path.join("{basedir}", "{genotype}_2.fq")
    params:
        seed = config['sim_reads']['seed'],
        profile = config['sim_reads']['profile'],
        read_len = config['sim_reads']['read_len'],
        insert_len = config['sim_reads']['insert_len'],
        coverage = max(config['sim_reads']['coverage']),
        prefix = os.path.join("{basedir}", "{genotype}_")
    conda:
        "../environment.yaml"
    shell:
        """
        set -xe

        art_illumina \
            -ss {params.profile} \
            -M \
            -i "{input}" \
            -p \
            -l {params.read_len[0]} \
            -f {params.coverage} \
            -s {params.insert_len[0]} \
            -m {params.insert_len[1]} \
            -na \
            -rs {params.seed} \
            -o "{params.prefix}"
        """
