#!/usr/bin/bash
#$ -q bioinfo.q
#$ -V
#$ -cwd
#$ -N concat-reads
#$ -pe shared 72
#$ -M obadbotanist@yahoo.com

# AUTHOR: Olabiyi Aderemi Obayomi
# EMAIL: obadbotanist@yahoo.com


set -e

source /gpfs0/bioinfo/apps/Miniconda2/Miniconda_v4.3.21/bin/activate /gpfs0/bioinfo/apps/Miniconda2/Miniconda_v4.3.21/envs/Metagenomics
export PERL5LIB='/gpfs0/bioinfo/apps/Miniconda2/Miniconda_v4.3.21/envs/Metagenomics/lib/site_perl/5.26.2/x86_64-linux-thread-multi'

PROJECT_DIR='/gpfs0/bioinfo/projects/Amit_Gross/03.Metagenomics_stormwater_biofiltration/'
METADATA="${PROJECT_DIR}/metadata.tsv"

# This assumes that the first column contains your sample names
SAMPLES=($(cat ${METADATA} | awk 'NR>1{print $1}'))
RAW_DATA_DIR="${PROJECT_DIR}/03.remove_host/fastq_files/"
OUT_DIR="${PROJECT_DIR}/04.concated_reads/"


# Assumes xyou fastqz files begin with the sample name
function cat_fastq(){

	local SAMPLE=$1
        local READ_DIR=$2
        local OUT_DIR=$3

	FILES=($(find ${READ_DIR} -name "${SAMPLE}*" -type f))

        cat ${FILES[*]} > ${OUT_DIR}/${SAMPLE}.concated.fastq

}


[-f ${OUT_DIR} ] || mkdir ${OUT_DIR}
export -f cat_fastq
parallel --jobs 5 "cat_fastq {} ${RAW_DATA_DIR} ${OUT_DIR}" ::: ${SAMPLES[*]}
