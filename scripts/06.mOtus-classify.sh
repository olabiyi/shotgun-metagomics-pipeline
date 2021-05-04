#!/usr/bin/env bash
#$ -S /bin/bash
#$ -q bioinfo.q
#$ -V
#$ -cwd
#$ -N classify_mOtus
#$ -pe shared 72
#$ -M obadbotanist@yahoo.com


# AUTHOR: Olabiyi Aderemi Obayomi
# Email: obadbotanist@yahoo.com
set -e

source /gpfs0/bioinfo/users/obayomi/miniconda3/bin/activate base
export PERL5LIB="/gpfs0/bioinfo/users/obayomi/miniconda3/lib/5.32.0"

PROJECT_DIR='/gpfs0/bioinfo/projects/Amit_Gross/03.Metagenomics_stormwater_biofiltration/'
RAW_DATA_DIR="${PROJECT_DIR}/03.remove_host/fastq_files/"
SAMPLES=($(cat ${PROJECT_DIR}/metadata.tsv | awk 'NR>1{print $1}'))
MOTUS_CLASSIFY="${PROJECT_DIR}/04.read_based_analysis/07.mOtus-classify/"
MOTUS_SNV="${PROJECT_DIR}/04.read_based_analysis/08.mOtus-snv/"
FASTQ_SUFFIX='_host_removed.fastq.gz'

TAXON_LEVELS=(phylum class order family genus mOTU)
MOTUS_BIN='/gpfs0/bioinfo/users/obayomi/miniconda3/bin/motus'


function motus_profile() {
	
	local SAMPLE=$1
	local MOTUS_CLASSIFY=$2
	local RAW_DATA_DIR=$3
	local TAXON_LEVELS=(phylum class order family genus mOTU)

	for TAXON_LEVEL in ${TAXON_LEVELS[*]}; do
			${MOTUS_BIN} profile \
				-k ${TAXON_LEVEL} \
                        	-c \
                        	-t 10 \
                        	-s ${RAW_DATA_DIR}/${SAMPLE}${FASTQ_SUFFIX} \
                        	-o ${MOTUS_CLASSIFY}/${SAMPLE}/${SAMPLE}.${TAXON_LEVEL}.motus \
                        	-n ${SAMPLE}
	done

}

export -f motus_profile

# Make directories
parallel --jobs 5 "[ -d ${MOTUS_CLASSIFY}/{} ] || mkdir ${MOTUS_CLASSIFY}/{}" ::: ${SAMPLES[*]}

# Profile taxonomy - out results as counts not relative abundance using the -c argument
parallel --jobs 5 "motus_profile {} ${MOTUS_CLASSIFY} ${RAW_DATA_DIR}" ::: ${SAMPLES[*]} 

# Merge profiles
for TAXON_LEVEL in ${TAXON_LEVELS[*]}; do

	FILES=($(find ${MOTUS_CLASSIFY}/ -type f -name "*.${TAXON_LEVEL}.motus"))
	FILES=$(echo ${FILES[*]} |sed 's/ /,/g')
	# Merge profiles
	${MOTUS_BIN} merge -i ${FILES} -o ${MOTUS_CLASSIFY}/merged.${TAXON_LEVEL}.motus

done

# Generating single nucleotide variant (SNV) profiles using MGs
# Map reads
#parallel --jobs 5 "${MOTUS_BIN} map_snv -s ${RAW_DATA_DIR}/{}${FASTQ_SUFFIX}  > ${MOTUS_SNV}/{}.bam" ::: ${SAMPLES[*]}
# Call variants
#${MOTUS_BIN} snv_call -d ${MOTUS_SNV} -o ${MOTUS_SNV}
