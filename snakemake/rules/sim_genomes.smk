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

        SURVIVOR simSV "{output}" &&
        sed -E -i.org "s/^(TRANSLOCATION_number:)\s+[0-9]+/\\1 {params.n_trans}/" "{output}" &&
        cat "{output}"
        """

rule samtools_faidx:
    input:
        config['input']['fasta']
    output:
        [os.path.join("{basedir}", "seqids.") + sfx for sfx in ("txt", "fasta")]
    params:
        seqids = "\n".join([str(c) for c in config['input']['seqids']])
    conda:
        "../environment.yaml"
    shell:
        """
        set -xe
        echo "{params.seqids}" > "{output[0]}"
        if [ "{params.seqids}" == "" ]; then
            ln -sr "{input}" "{output[1]}"
        else
            samtools faidx "{input}" -r "{output[0]}" -o "{output[1]}"
        fi
        """

rule survivor_simsv:
    input:
        config = config['input']['config'],
        fasta = os.path.join("{basedir}", "seqids.fasta")
    output:
        fasta = os.path.join("{basedir}", "{genotype}.fasta")
    params:
        sfx = '.org',
        prefix = os.path.join("{basedir}", "{genotype}")
    conda:
        "../environment.yaml"
    shell:
        """
        set -xe

        # write FASTA files with diploid genomes
        # N.B.: sequence IDs must be unique -> [SEQID].[N]
        if [ "{wildcards.genotype}" == "hmz" ]; then
            sed -E "s:^(>.*):\\1\.1:" "{input.fasta}" > "{output.fasta}"
            sed -E "s:^(>.*):\\1\.2:" "{input.fasta}" >> "{output.fasta}"
        elif [ "{wildcards.genotype}" == "hmz-sv" ]; then
            SURVIVOR simSV "{input.fasta}" "{input.config}" 0 0 "{params.prefix}"
            sed -E -i{params.sfx} "s:^(>.*):\\1\.1:" "{output.fasta}"
            sed -E "s:^(>.*):\\1\.2:" "{output.fasta}{params.sfx}" >> "{output.fasta}"
            rm -f "{output.fasta}{params.sfx}"
        elif [ "{wildcards.genotype}" == "htz-sv" ]; then
            SURVIVOR simSV "{input.fasta}" "{input.config}" 0 0 "{params.prefix}"
            sed -E -i{params.sfx} "s:^(>.*):\\1\.1:" "{output.fasta}"
            sed -E "s:^(>.*):\\1\.2:" "{input.fasta}" >> "{output.fasta}"
            rm -f "{output.fasta}{params.sfx}"
        fi
        """
