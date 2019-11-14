rule bwa_index:
    input:
        os.path.join("{basedir}", "{genotype}.fasta")
    output:
        [os.path.join("{basedir}", "{genotype}.fasta.") + ext for ext in ['bwt', 'amb', 'ann', 'pac', 'sa']]
    conda:
        "../environment.yaml"
    shell:
        """
        set -xe

        bwa index "{input}" 
        """
    
rule bwa_mem:
    input:
        fasta = os.path.join("{basedir}", "{genotype}.fasta"),
        fastq1 = os.path.join("{basedir}", "{genotype}_1.fq"),
        fastq2 = os.path.join("{basedir}", "{genotype}_2.fq")
    params:
        read_group = "@RG\\tID:{0}\\tLB:{0}\\tSM:{0}".format("{genotype}")
    output:
        os.path.join("{basedir}", "{genotype}.sam")
    conda:
        "../environment.yaml"
    shell:
        """
        set -xe

        bwa mem \
            -R "{params.read_group}" \
            "{input.fasta}" "{input.fastq1}" "{input.fastq2}" > "{output}"
        """
