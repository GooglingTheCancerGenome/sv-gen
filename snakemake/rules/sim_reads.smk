rule art_illumina:
    input:
        fasta = "{prefix}.fasta"
    params:
        seed = config['sim_reads']['seed'],
        profile = config['sim_reads']['profile'],
        read_len = config['sim_reads']['read_len'],
        stdev = config['sim_reads']['stdev'],
        coverage = config['sim_reads']['coverage'],
        insert_size = config['sim_reads']['insert_size']
    output:
        fastq1 = '{prefix}_1.fq',
        fastq2 = '{prefix}_2.fq'
    conda:
        "../environment.yaml"
    shell:
        """
        set -xe
        PREFIX="$(basename "{input.fasta}" .fasta)"
        OUTDIR="$(dirname "{input.fasta}")"

        art_illumina \
            -ss {params.profile} \
            -M \
            -i {input.fasta} \
            -p \
            -l {params.read_len} \
            -f {params.coverage} \
            -m {params.insert_size} \
            -s {params.stdev} \
            -na \
            -rs {params.seed} \
            -o ${{OUTDIR}}/${{PREFIX}}_
        """
