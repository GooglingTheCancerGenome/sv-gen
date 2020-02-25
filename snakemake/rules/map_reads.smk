rule bwa_index:
    input:
        fasta = os.path.join("{basedir}", "seqids.fasta")
    output:
        fastai = [os.path.join("{basedir}", "seqids") +
                  e for e in config['file_exts']['fasta_idx']]
    conda:
        "../environment.yaml"
    shell:
        """
        set -xe

        bwa index "{input.fasta}" &&
        ls {output.fastai}
        """

rule bwa_mem:
    input:
        fasta = os.path.join("{basedir}", "seqids.fasta"),
        fastai = [os.path.join("{basedir}", "seqids") +
                  e for e in config['file_exts']['fasta_idx']],
        fastq1 = os.path.join("{basedir}", "r{read_len}_i{insert_len}", "{genotype}_1.fq"),
        fastq2 = os.path.join("{basedir}", "r{read_len}_i{insert_len}", "{genotype}_2.fq")
    output:
        bam = os.path.join("{basedir}", "r{read_len}_i{insert_len}", "cov" + str(max(config['sim_reads']['coverage'])), "{genotype}.bam")
    params:
        read_group = "@RG\\tID:{0}\\tLB:{0}\\tSM:{0}".format("{genotype}")
    conda:
        "../environment.yaml"
    shell:
        """
        set -xe

        bwa mem \
            -R "{params.read_group}" \
            "{input.fasta}" "{input.fastq1}" "{input.fastq2}" | \
        samtools sort -o "{output.bam}" && \
        rm -f "{input.fastq1}" "{input.fastq2}"
        """

rule samtools_view:
    input:
        bam = os.path.join("{basedir}", "r{read_len}_i{insert_len}", "cov" + str(max(config['sim_reads']['coverage'])), "{genotype}.bam")
    output:
        bam = os.path.join("{basedir}", "r{read_len}_i{insert_len}", "cov{cov}", "{genotype}.bam")
    conda:
        "../environment.yaml"
    shell:
        """
        set -xe

        FRAC=$(LC_ALL=C printf "%.2f" $(bc <<< "scale=2; {wildcards.cov} / 100"))
        samtools view -s ${{FRAC}} "{input.bam}" -o "{output.bam}"
        """

rule samtools_index:
    input:
        bam = os.path.join("{basedir}", "r{read_len}_i{insert_len}", "cov{cov}", "{genotype}.bam")
    output:
        bai = os.path.join("{basedir}", "r{read_len}_i{insert_len}", "cov{cov}", "{genotype}.bam.bai")
    conda:
        "../environment.yaml"
    shell:
        """
        set -xe

        samtools index -b "{input.bam}" "{output.bai}"
        """
