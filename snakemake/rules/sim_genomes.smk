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

rule samtools_faidx:
    input:
        fasta = config['input']['fasta']
    output:
        region = config['input']['region'],
        fasta = os.path.splitext(config['input']['region'])[0] + ".fasta"
    params:
        chr = "\n".join([str(c) for c in config['input']['chr']])
    conda:
        "../environment.yaml"
    shell:
        """
        set -xe
        echo "{params.chr}" > "{output.region}"
        if [ "{params.chr}" == "" ]; then
            ln -sr "{input.fasta}" "{output.fasta}"
        else
            samtools faidx "{input.fasta}" -r "{output.region}" -o "{output.fasta}"
        fi
        """

rule survivor_simsv:
    input:
        config = config['input']['config'],
        fasta = os.path.splitext(config['input']['region'])[0] + ".fasta"
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
