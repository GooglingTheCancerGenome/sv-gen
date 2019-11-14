rule bwa_index:
    input:
        os.path.join("{basedir}", "{genotype}.fasta")
    output:
        [os.path.join("{basedir}", "{genotype}") +
         e for e in config['file_exts']['fastai'][1:]]
    conda:
        "../environment.yaml"
    shell:
        """
        set -xe

        bwa index "{input}"
        """
    
rule bwa_mem:
    input:
        fastai = [os.path.join("{basedir}", "{genotype}") + e
                  for e in config['file_exts']['fastai']],
        fastq1 = os.path.join("{basedir}", "{genotype}_1.fq"),
        fastq2 = os.path.join("{basedir}", "{genotype}_2.fq")
    output:
        os.path.join("{basedir}", "{genotype}.sam")
    params:
        read_group = "@RG\\tID:{0}\\tLB:{0}\\tSM:{0}".format("{genotype}")
    conda:
        "../environment.yaml"
    shell:
        """
        set -xe

        bwa mem \
            -R "{params.read_group}" \
            "{input.fastai[0]}" "{input.fastq1}" "{input.fastq2}" > "{output}"
        """
