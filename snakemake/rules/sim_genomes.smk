rule survivor_config:
    output:
        config['input']['config']
    params:
        matrix = config['sim_genomes']['sv_type']
    conda:
        "../environment.yaml"
    shell:
        """
        set -xe

        SURVIVOR simSV "{output}" &&
        sed -i.org \
            -E "s/^(DUPLICATION_number:)\s+[0-9]+/\\1 {params.matrix[DUP][0]}/;\
                s/^(DUPLICATION_minimum_length:)\s+[0-9]+/\\1 {params.matrix[DUP][1]}/;\
                s/^(DUPLICATION_maximum_length:)\s+[0-9]+/\\1 {params.matrix[DUP][2]}/;\
                s/^(INDEL_number:)\s+[0-9]+/\\1 {params.matrix[INDEL][0]}/;\
                s/^(INDEL_minimum_length:)\s+[0-9]+/\\1 {params.matrix[INDEL][1]}/;\
                s/^(INDEL_maximum_length:)\s+[0-9]+/\\1 {params.matrix[INDEL][2]}/;\
                s/^(INVERSION_number:)\s+[0-9]+/\\1 {params.matrix[INV][0]}/;\
                s/^(INVERSION_minimum_length:)\s+[0-9]+/\\1 {params.matrix[INV][1]}/;\
                s/^(INVERSION_maximum_length:)\s+[0-9]+/\\1 {params.matrix[INV][2]}/;\
                s/^(TRANSLOCATION_number:)\s+[0-9]+/\\1 {params.matrix[TRA][0]}/;\
                s/^(TRANSLOCATION_minimum_length:)\s+[0-9]+/\\1 {params.matrix[TRA][1]}/;\
                s/^(TRANSLOCATION_maximum_length:)\s+[0-9]+/\\1 {params.matrix[TRA][2]}/"\
            "{output}"
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
