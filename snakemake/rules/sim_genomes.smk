rule survivor_config:
    output:
        config = os.path.join(get_outdir(), '{svtype}',
                              config['simulation']['config'])
    params:
        matrix = config['simulation']['sv_type']
    conda:
        "../environment.yaml"
    shell:
        """
        set -xe

        SURVIVOR simSV "{output.config}"
        sed -i.org \
            -E "s/^(DUPLICATION_number:)\s+[0-9]+/\\1 {params.matrix[dup][0]}/;\
                s/^(DUPLICATION_minimum_length:)\s+[0-9]+/\\1 {params.matrix[dup][1]}/;\
                s/^(DUPLICATION_maximum_length:)\s+[0-9]+/\\1 {params.matrix[dup][2]}/;\
                s/^(INVERSION_number:)\s+[0-9]+/\\1 {params.matrix[inv][0]}/;\
                s/^(INVERSION_minimum_length:)\s+[0-9]+/\\1 {params.matrix[inv][1]}/;\
                s/^(INVERSION_maximum_length:)\s+[0-9]+/\\1 {params.matrix[inv][2]}/;\
                s/^(INDEL_number:)\s+[0-9]+/\\1 {params.matrix[indel][0]}/;\
                s/^(INDEL_minimum_length:)\s+[0-9]+/\\1 {params.matrix[indel][1]}/;\
                s/^(INDEL_maximum_length:)\s+[0-9]+/\\1 {params.matrix[indel][2]}/;\
                s/^(TRANSLOCATION_number:)\s+[0-9]+/\\1 {params.matrix[tra][0]}/;\
                s/^(TRANSLOCATION_minimum_length:)\s+[0-9]+/\\1 {params.matrix[tra][1]}/;\
                s/^(TRANSLOCATION_maximum_length:)\s+[0-9]+/\\1 {params.matrix[tra][2]}/;\
                s/^(INV_del_number:)\s+[0-9]+/\\1 {params.matrix[invdel][0]}/;\
                s/^(INV_del_minimum_length:)\s+[0-9]+/\\1 {params.matrix[invdel][1]}/;\
                s/^(INV_del_maximum_length:)\s+[0-9]+/\\1 {params.matrix[invdel][2]}/;\
                s/^(INV_dup_number:)\s+[0-9]+/\\1 {params.matrix[invdup][0]}/;\
                s/^(INV_dup_minimum_length:)\s+[0-9]+/\\1 {params.matrix[invdup][1]}/;\
                s/^(INV_dup_maximum_length:)\s+[0-9]+/\\1 {params.matrix[invdup][2]}/"\
            "{output}"
        cat "{output}"
        """

rule samtools_faidx:
    input:
        fasta = get_reference()
    output:
        seqids = os.path.join(get_outdir(), 'seqids.txt'),
        fasta = os.path.join(get_outdir(), 'seqids' + get_filext('fasta'))
    params:
        seqids = "\n".join([str(c) for c in config['input']['seqids']])
    conda:
        "../environment.yaml"
    shell:
        """
        set -xe
        echo "{params.seqids}" > "{output.seqids}"
        if [ "{params.seqids}" == "" ]; then
            ln -sr "{input.fasta}" "{output.fasta}"
        else
            samtools faidx "{input.fasta}" -r "{output.seqids}" -o "{output.fasta}"
        fi
        """

rule survivor_simsv:
    input:
        config = os.path.join(get_outdir(), '{svtype}',
                              config['simulation']['config']),
        fasta = os.path.join(get_outdir(), 'seqids' + get_filext('fasta'))
    output:
        fasta = os.path.join(get_outdir(), '{svtype}',
                             '{genotype}' + get_filext('fasta')),
        vcf = os.path.join(get_outdir(), '{svtype}',
                             '{genotype}' + get_filext('vcf')),
        bed = os.path.join(get_outdir(), '{svtype}',
                             '{genotype}' + get_filext('bed'))
    params:
        sfx = '.org',
        prefix = os.path.join(get_outdir(), '{svtype}', '{genotype}')
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
            # write dummy VCF and BED output files (no SVs)
            touch "{output.vcf}" "{output.bed}"
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
