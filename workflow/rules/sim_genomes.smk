rule survivor_config:
    output:
        config = os.path.join(config.output.basedir, '{svtype}',
                              config.simulation.config)
    params:
        m = config.simulation.svtype
    conda:
        "../environment.yaml"
    shell:
        """
        set -xe

        SURVIVOR simSV "{output.config}"
        sed -i.org \
            -E "s/^(DUPLICATION_number:)\s+[0-9]+/\\1 {params.m.dup.count}/;\
                s/^(DUPLICATION_minimum_length:)\s+[0-9]+/\\1 {params.m.dup.min_len}/;\
                s/^(DUPLICATION_maximum_length:)\s+[0-9]+/\\1 {params.m.dup.max_len}/;\
                s/^(INVERSION_number:)\s+[0-9]+/\\1 {params.m.inv.count}/;\
                s/^(INVERSION_minimum_length:)\s+[0-9]+/\\1 {params.m.inv.min_len}/;\
                s/^(INVERSION_maximum_length:)\s+[0-9]+/\\1 {params.m.inv.max_len}/;\
                s/^(INDEL_number:)\s+[0-9]+/\\1 {params.m.indel.count}/;\
                s/^(INDEL_minimum_length:)\s+[0-9]+/\\1 {params.m.indel.min_len}/;\
                s/^(INDEL_maximum_length:)\s+[0-9]+/\\1 {params.m.indel.max_len}/;\
                s/^(TRANSLOCATION_number:)\s+[0-9]+/\\1 {params.m.tra.count}/;\
                s/^(TRANSLOCATION_minimum_length:)\s+[0-9]+/\\1 {params.m.tra.min_len}/;\
                s/^(TRANSLOCATION_maximum_length:)\s+[0-9]+/\\1 {params.m.tra.max_len}/;\
                s/^(INV_del_number:)\s+[0-9]+/\\1 {params.m.invdel.count}/;\
                s/^(INV_del_minimum_length:)\s+[0-9]+/\\1 {params.m.invdel.min_len}/;\
                s/^(INV_del_maximum_length:)\s+[0-9]+/\\1 {params.m.invdel.max_len}/;\
                s/^(INV_dup_number:)\s+[0-9]+/\\1 {params.m.invdup.count}/;\
                s/^(INV_dup_minimum_length:)\s+[0-9]+/\\1 {params.m.invdup.min_len}/;\
                s/^(INV_dup_maximum_length:)\s+[0-9]+/\\1 {params.m.invdup.max_len}/"\
            "{output}"
        cat "{output}"
        """

rule samtools_faidx:
    input:
        fasta = get_reference()
    output:
        seqids = os.path.join(config.output.basedir, 'seqids.txt'),
        fasta = os.path.join(config.output.basedir, 'seqids' + config.filext.fasta)
    params:
        seqids = "\n".join([str(s) for s in get_seqids()])
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
        config = os.path.join(config.output.basedir, '{svtype}',
                              config.simulation.config),
        fasta = os.path.join(config.output.basedir, 'seqids' + config.filext.fasta)
    output:
        fasta = os.path.join(config.output.basedir, '{svtype}',
                             '{genotype}' + config.filext.fasta),
        vcf = os.path.join(config.output.basedir, '{svtype}',
                             '{genotype}' + config.filext.vcf),
        bed = os.path.join(config.output.basedir, '{svtype}',
                             '{genotype}' + config.filext.bed)
    params:
        sfx = '.org',
        prefix = os.path.join(config.output.basedir, '{svtype}', '{genotype}')
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
