rule bwa_index:
    input:
        fasta = os.path.join(get_outdir(), 'seqids' + get_filext('fasta'))
    output:
        fastai = [os.path.join(get_outdir(), 'seqids') +
                  e for e in get_filext('fasta_idx')]
    conda:
        "../environment.yaml"
    shell:
        """
        set -xe

        bwa index "{input.fasta}"
        ls {output.fastai}
        """

rule bwa_mem:
    input:
        fasta = os.path.join(get_outdir(), 'seqids' + get_filext('fasta')),
        fastai = [os.path.join(get_outdir(), 'seqids') +
                  e for e in get_filext('fasta_idx')],
        fastq1 = os.path.join(get_outdir(), '{svtype}',
                              '{read_len}_{insert_len}',
                              '{genotype}_1' + get_filext('fastq')),
        fastq2 = os.path.join(get_outdir(), '{svtype}',
                              '{read_len}_{insert_len}',
                              '{genotype}_2' + get_filext('fastq'))
    output:
        bam = os.path.join(get_outdir(), '{svtype}',
                           '{read_len}_{insert_len}',
                           'cov' + str(max(config['sim_reads']['coverage'])),
                           '{genotype}' + get_filext('bam'))
    params:
        read_group = "@RG\\tID:{0}\\tLB:{0}\\tSM:{0}".format('{genotype}')
    conda:
        "../environment.yaml"
    threads:
        get_nthreads()
    shell:
        """
        set -xe

        bwa mem \
            -t {threads} \
            -R "{params.read_group}" \
            "{input.fasta}" "{input.fastq1}" "{input.fastq2}" | \
        samtools sort -@ {threads} -o "{output.bam}"
        rm -f "{input.fastq1}" "{input.fastq2}"
        """

rule samtools_view:
    input:
        bam = os.path.join(get_outdir(), '{svtype}',
                           '{read_len}_{insert_len}',
                           'cov' + str(max(config['sim_reads']['coverage'])),
                           '{genotype}' + get_filext('bam'))
    output:
        bam = os.path.join(get_outdir(), '{svtype}',
                           '{read_len}_{insert_len}', 'cov{cov}',
                           '{genotype}' + get_filext('bam'))
    conda:
        "../environment.yaml"
    threads:
        get_nthreads()
    shell:
        """
        set -xe

        FRAC=$(LC_ALL=C printf "%.2f" $(bc <<< "scale=2; {wildcards.cov} / 100"))
        samtools view -@ {threads} -s ${{FRAC}} "{input.bam}" -o "{output.bam}"
        """

rule samtools_index:
    input:
        bam = os.path.join(get_outdir(), '{svtype}',
                           '{read_len}_{insert_len}', 'cov{cov}',
                           '{genotype}' + get_filext('bam'))
    output:
        bai = os.path.join(get_outdir(), '{svtype}',
                           '{read_len}_{insert_len}', 'cov{cov}',
                           '{genotype}' + get_filext('bam_idx'))
    conda:
        "../environment.yaml"
    threads:
        get_nthreads()
    shell:
        """
        set -xe

        samtools index -@ {threads} -b "{input.bam}" "{output.bai}"
        """
