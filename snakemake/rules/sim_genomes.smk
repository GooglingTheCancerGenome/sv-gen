rule survivor_config:
    output:
        config['input']['config']
    params:
        n_trans = 0
    conda:
        "../environment.yaml"
    shell:
        """
        set -xe
        SURVIVOR simSV "{output}"
        sed -E -i.org "s/^(TRANSLOCATION_number:)\s+[0-9]+/\\1 {params.n_trans}/" "{output}"
        cat "{output}"
        """

rule survivor_simsv:
    input:
        config = config['input']['config'],
        fasta = config['input']['fasta']
    output:
        "{basedir}/{genotype}/{prefix}.fasta"
    params:
        sfx = '.org'
    conda:
        "../environment.yaml"
    shell:
        """
        set -xe
        OUTDIR="$(dirname "{output}")"
        PREFIX="${{OUTDIR}}/$(basename "{output}" .fasta)"

        # write diploid genomes into FASTA
        # sequence IDs must be unique -> [SEQID].[N]
        if [ "{wildcards.genotype}" == "hmzsv" ]; then
            SURVIVOR simSV "{input.fasta}" "{input.config}" 0 0 "${{PREFIX}}"
            sed -E -i{params.sfx} "s:^(>.*):\\1\.1:" "{output}"
            sed -E "s:^(>.*):\\1\.2:" "{output}{params.sfx}" >> "{output}"
            rm -f "{output}{params.sfx}"
        elif [ "{wildcards.genotype}" == "htzsv" ]; then
            SURVIVOR simSV "{input.fasta}" "{input.config}" 0 0 "${{PREFIX}}"
            sed -E -i{params.sfx} "s:^(>.*):\\1\.1:" "{output}"
            sed -E "s:^(>.*):\\1\.2:" "{input.fasta}" >> "{output}"
            rm -f "{output}{params.sfx}"
        elif [ "{wildcards.genotype}" == "hmz" ]; then
            sed -E "s:^(>.*):\\1\.1:" "{input.fasta}" > "{output}"
            sed -E "s:^(>.*):\\1\.2:" "{input.fasta}" >> "{output}"
        fi
        """
