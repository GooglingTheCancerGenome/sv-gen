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
        stdev = config['sim_reads']['stdev'],
        coverage = config['sim_reads']['coverage'],
        insert_size = config['sim_reads']['insert_size'],
        prefix = os.path.join("{basedir}", "{genotype}_")
    conda:
        "../environment.yaml"
    shell:
        """
        set -xe

        art_illumina \
            -ss {params.profile} \
            -M \
            -i {input} \
            -p \
            -l {params.read_len} \
            -f {params.coverage} \
            -m {params.insert_size} \
            -s {params.stdev} \
            -na \
            -rs {params.seed} \
            -o {params.prefix}
        """
