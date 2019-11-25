rule bwa_index:
    input:
        fasta = os.path.join("{basedir}", "{genotype}.fasta")
    output:
        fastai = [os.path.join("{basedir}", "{genotype}") +
                  e for e in config['file_exts']['fastai'][1:]]
    conda:
        "../environment.yaml"
    shell:
        """
        set -xe

        bwa index "{input.fasta}"
        """
    
rule bwa_mem:
    input:
        fastai = [os.path.join("{basedir}", "{genotype}") +
                  e for e in config['file_exts']['fastai']],
        fastq1 = os.path.join("{basedir}", "{genotype}_1.fq"),
        fastq2 = os.path.join("{basedir}", "{genotype}_2.fq")
    output:
        bam = os.path.join("{basedir}", "{genotype}.bam")
    params:
        read_group = "@RG\\tID:{0}\\tLB:{0}\\tSM:{0}".format("{genotype}")
    conda:
        "../environment.yaml"
    shell:
        """
        set -xe

        bwa mem \
            -R "{params.read_group}" \
            "{input.fastai[0]}" "{input.fastq1}" "{input.fastq2}" | \
        samtools sort -o "{output.bam}"
        """

rule samtools_view:
    input:
        bam = os.path.join("{basedir}", "{genotype}.bam")
    output:
        bam = os.path.join("{basedir}", "{genotype}.cov{cov}.bam")
    conda:
        "../environment.yaml"
    shell:
        """
        set -xe

        FRAC=$(LC_ALL=C printf "%.2f" $(bc <<< "scale=2; {wildcards.cov} / 100"))
        samtools view -s ${{FRAC}} "{input.bam}" -o "{output.bam}"
        """
