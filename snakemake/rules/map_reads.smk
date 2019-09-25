rule bwa_index:
    input:
        fasta = config['input']['fasta']
    output:
        fasta_idx = [config['input']['fasta'] + "." + ext for ext in ['bwt', 'amb', 'ann', 'pac', 'sa']]
    conda:
        "../environment.yaml"
    shell:
        """
        set -xe

        bwa index "{input.fasta}" 
        """
    
rule bwa_mem:
    input:
        fasta = config['input']['fasta'],
        fasta_idx = [config['input']['fasta'] + "." + ext for ext in ['bwt', 'amb', 'ann', 'pac', 'sa']],
        fastq1 = "{path}/{prefix}_1.fq",
        fastq2 = "{path}/{prefix}_2.fq"
    params:
        read_group = "@RG\\tID:{0}\\tLB:{0}\\tSM:{0}".format("{prefix}")
    output:
        sam = "{path}/{prefix}.sam"
    conda:
        "../environment.yaml"
    shell:
        """
        set -xe

        bwa mem \
            -R "{params.read_group}" \
            "{input.fasta}" "{input.fastq1}" "{input.fastq2}" > "{output.sam}"
        """
