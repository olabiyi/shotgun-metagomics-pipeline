#!/usr/bin/env bash
#$ -S /bin/bash
#$ -q bioinfo.q
#$ -V
#$ -cwd
#$ -N find_protein
#$ -pe shared 72

# A script to find a set of proteins of interest
# AUTHOR: Olabiyi Aderemi Obayomi
# Email: obadbotanist@yahoo.com 


set -e

# Set the varaziable appropriately
source /gpfs0/bioinfo/apps/Miniconda2/Miniconda_v4.3.21/bin/activate /gpfs0/bioinfo/apps/Miniconda2/Miniconda_v4.3.21/envs/Metagenomics
export PERL5LIB='/gpfs0/bioinfo/apps/Miniconda2/Miniconda_v4.3.21/envs/Metagenomics/lib/site_perl/5.26.2/x86_64-linux-thread-multi'

PROJECT_DIR='/gpfs0/bioinfo/projects/Amit_Gross/03.Metagenomics_stormwater_biofiltration/'
METADATA="${PROJECT_DIR}/metadata.tsv"
# Assumes your sample names correspon to the first column of your metadata file
SAMPLES=($(cat ${METADATA} | awk 'NR>1{print $1}'))
RAW_DATA_DIR="${PROJECT_DIR}/03.remove_host/fastq_files/"
OUT_DIR="${PROJECT_DIR}/04.read_based_analysis/10.find_protein/"
DB_DIR="${PROJECT_DIR}/proteins_of_interest/"
PROTEINS="${DB_DIR}/proteins_of_interest.ffa"
DB="${DB_DIR}/proteins_of_interest.dmnd"
FASTQ_SUFFIX='_host_removed.fastq.gz'


# make diamond database
[ -f ${DB} ] || diamond makedb --in ${PROTEINS} --db ${DB%%.dmnd*}


# function to query a diamond database (protein or nucleotide) for sequences of interest.
# It returns a tabular summary of counts of sequences per hit at a chosen similarity threshold.
# By default hits with e-value <= 1e-5 and a minimum query coverage of 60% percent at a specified 
# sequence similarity threshod are returned.
function find_protein(){
	
	local READ_DIR=$1
	local SAMPLE=$2
	local OUT_DIR=$3
	local DB=$4
	local SEARCH=$5
	local SUFFIX=$6
	local PERCENT_ID=$7
	
	# Create a Directory for export
	[ -d ${OUT_DIR}/Export_${PERCENT_ID} ] || mkdir ${OUT_DIR}/Export_${PERCENT_ID}


	# Diamond blastx against proteins_of_interest
	diamond ${SEARCH} \
		--db ${DB} \
		--query ${READ_DIR}/${SAMPLE}${SUFFIX} \
		--out ${OUT_DIR}/Export_${PERCENT_ID}/${SAMPLE}_matches_${PERCENT_ID}.tsv \
		--outfmt 6 qseqid sseqid stitle pident length mismatch gapopen qstart qend sstart send evalue bitscore \
		--max-target-seqs 1 \
		--evalue 1e-5 \
		--id ${PERCENT_ID} \
		--threads 10 \
		--query-cover 60 \
		--more-sensitive

	awk -F "\t" '{print $3}' ${OUT_DIR}/Export_${PERCENT_ID}/${SAMPLE}_matches_${PERCENT_ID}.tsv |sort | uniq -c |sort -nr \
	> ${OUT_DIR}/${SAMPLE}/${SAMPLE}_count_${PERCENT_ID}.txt

}


export -f  find_protein

parallel --jobs 0 "[ -d ${OUT_DIR}/{} ] || mkdir ${OUT_DIR}/{}" ::: ${SAMPLES[*]} 

# 95% identity
parallel --jobs 10  "find_protein ${RAW_DATA_DIR} {} ${OUT_DIR} ${DB%%.dmnd*} blastx ${FASTQ_SUFFIX} 95"  ::: ${SAMPLES[*]}
# 90% identity
parallel --jobs 10  "find_protein ${RAW_DATA_DIR} {} ${OUT_DIR} ${DB%%.dmnd*} blastx ${FASTQ_SUFFIX} 90"  ::: ${SAMPLES[*]}
# 85% indentity
parallel --jobs 10  "find_protein ${RAW_DATA_DIR} {} ${OUT_DIR} ${DB%%.dmnd*} blastx ${FASTQ_SUFFIX} 85"  ::: ${SAMPLES[*]}

 
