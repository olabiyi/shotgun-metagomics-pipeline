#!/usr/bin/env bash
#$ -S /bin/bash
#$ -q bioinfo.q
#$ -V
#$ -cwd
#$ -N trim-reads
#$ -pe shared 72
#$ -M obadbotanist@yahoo.com


# AUTHOR: Olabiyi Aderemi Obayomi
# EMAIL: obadbotanist@yahoo.com


set -eo pipefail
source activate /gpfs0/bioinfo/users/obayomi/miniconda3/envs/bioinfo
export PERL5LIB='/gpfs0/bioinfo/users/obayomi/miniconda3/envs/bioinfo/lib/site_perl/5.26.2/x86_64-linux-thread-multi' 

PROJECT_DIR='/gpfs0/bioinfo/projects/Amit_Gross/03.Metagenomics_stormwater_biofiltration/'
PAIRED='True'
RAW_DATA_DIR="${PROJECT_DIR}/03.remove_host/fastq_files/"
TRIM_DIR="${PROJECT_DIR}/02.trim_reads"
FORWARD_SUFFIX="_R1.fastq.gz"
REVERSE_SUFFIX="_R2.fastq.gz"
ADAPTORS='/fastspace/bioinfo_apps/Trimmomatic-0.32/adapters/TruSeq3-PE.fa'
METADATA="${PROJECT_DIR}/metadata.tsv"
SAMPLES=($(cat ${METADATA} | awk 'NR>1{print $1}'))
MIN_LENGTH=50

# Trim multiple fastq files

if [ ${PAIRED} == 'True' ]; then

	parallel --jobs 5 	"trimmomatic-0.32.jar  PE -phred33 \
					${RAW_DATA_DIR}/{}${FORWARD_SUFFIX} ${RAW_DATA_DIR}/{}${REVERSE_SUFFIX} \
					${TRIM_DIR}/{}.F.trim.fastq ${TRIM_DIR}/{}.F.un.trim.fastq \
					${TRIM_DIR}/{}.R.trim.fastq ${TRIM_DIR}/{}.R.un.trim.fastq \
					ILLUMINACLIP:${ADAPTORS}:2:30:10 \
					LEADING:3 \
					TRAILING:3 \
					SLIDINGWINDOW:4:20 \
					MINLEN:${MIN_LENGTH}" ::: ${SAMPLES[*]}

else

	parallel --jobs 5	"trimmomatic-0.32.jar  SE -phred33 \
					${RAW_DATA_DIR}/{}${FORWARD_SUFFIX} \
					${TRIM_DIR}/{}.trim.fastq ${TRIM_DIR}/{}.un.trim.fastq \
					ILLUMINACLIP:${ADAPTORS}:2:30:10 \
					LEADING:3 \
					TRAILING:3 \
					SLIDINGWINDOW:4:20 \
					MINLEN:${MIN_LENGTH}" ::: ${SAMPLES[*]}


fi
