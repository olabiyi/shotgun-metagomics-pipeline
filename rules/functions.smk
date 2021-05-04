from glob import glob

def get_fq1(wildcards):
    # code that returns a list of fastq files for read 1 based on *wilcards.sample* e.g
    return sorted(glob(wilcards.sample + '*_R1_001.fastq.gz'))

def get_fq2(wildcards):
    # code that returns a list of fastq files for read 2 based on *wilcards.sample* e.g
    return sorted(glob(wilcards.sample + '*_R2_001.fastq.gz'))

GENOME='genome.fa'

# usage

rule bwa_map:
    input:
        fq1=get_fq1,
        fq2=get_fq2
    output:
        bam='{sample}.bam'
    shell:
        "bwa mem {GENEOME} {input} | samtools sort -o {output}"
