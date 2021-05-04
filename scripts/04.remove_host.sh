#!/usr/bin/env bash
#$ -S /bin/bash
#$ -q bioinfo.q
#$ -V
#$ -cwd
#$ -N remove-contaminants
#$ -pe shared 72

set -e

# Your project directory
PROJECT_DIR='/gpfs0/bioinfo/projects/Amit_Gross/03.Metagenomics_stormwater_biofiltration'
HOST_FASTA="${PROJECT_DIR}/genomes/PHIX_and_host_sequences.fna"
OUT_DIR="${PROJECT_DIR}/03.remove_host/"
FASTQ_SUFFIX='_host_removed.fastq.qz'

# Build bowtie database
/gpfs0/bioinfo/apps/bowtie2/bowtie2-2.3.5-linux-x86_64/bowtie2-build ${HOST_FASTA} host_DB

# The directory contain the reads that were quality trimmed using trimmomatic
TRIM_DIR="${PROJECT_DIR}/02.trim_reads/data/trimmo/Trim_imported_reads"

FILES=($(ls -1 ${TRIM_DIR}/*.fq))

# This assumes that the reads quality trimming was managed using NeatSeq_Flow
SAMPLES=($(basename -a -s ".Single.import.fastq.trimmo.fq" ${FILES[*]}))

source activate /gpfs0/bioinfo/apps/Miniconda2/Miniconda_v4.3.21/envs/NeatSeq_Flow
# To avoid conflict between differnt perl versions that may be installed
export PERL5LIB="/gpfs0/bioinfo/apps/Miniconda2/Miniconda_v4.3.21/envs/NeatSeq_Flow/lib/site_perl/5.26.2/x86_64-linux-thread-multi"


[ -d fastq_files/ ] || mkdir fastq_files/

# Remove Phix and Carex DNA
parallel \
	--jobs 0 \
	--link  \
	"bowtie2 -p 8 -x host_DB -U {1} --un-gz fastq_files/{2}${FASTQ_SUFFIX}  > {2}_mapped_and_unmapped.sam" \
	::: ${FILES[*]} ::: ${SAMPLES[*]}


# Remove sam file to free memory
rm -rf ./*.sam
