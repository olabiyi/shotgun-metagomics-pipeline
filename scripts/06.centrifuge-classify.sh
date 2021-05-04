#!/usr/bin/env bash
#$ -S /bin/bash
#$ -q bioinfo.q
#$ -V
#$ -cwd
#$ -N classify_centrifuge
#$ -pe shared 72
#$ -M obadbotanist@yahoo.com

set -e

source /gpfs0/bioinfo/apps/Miniconda2/Miniconda_v4.3.21/bin/activate /gpfs0/bioinfo/apps/Miniconda2/Miniconda_v4.3.21/envs/Metagenomics
export PERL5LIB='/gpfs0/bioinfo/apps/Miniconda2/Miniconda_v4.3.21/envs/Metagenomics/lib/site_perl/5.26.2/x86_64-linux-thread-multi'
PROJECT_DIR='/gpfs0/bioinfo/projects/Amit_Gross/03.Metagenomics_stormwater_biofiltration/'
CENTRIFUGE_DB='/gpfs0/bioinfo/databases/Centrifuge/db_arch_bact_vir_2018/'
CENTRIFUGE_CLASSIFY="${PROJECT_DIR}/04.read_based_analysis/09.centrifuge-classify/"
METADATA="${PROJECT_DIR}/metadata.tsv"
FASTQ_SUFFIX='_host_removed.fastq.gz'
# This assumes that the first column contains your sample names
SAMPLES=($(cat ${METADATA} | awk 'NR>1{print $1}'))
RAW_DATA_DIR="${PROJECT_DIR}/03.remove_host/fastq_files/"

KRONA_TAXONOMY='/gpfs0/bioinfo/apps/Miniconda2/Miniconda_v4.3.21/envs/Metagenomics/databases/krona/taxonomy'

function centrifuge_classify(){
	
	local SAMPLE=$1
	local OUT_DIR=$2
	local FORWARD=$3
	local REVERSE=$4
	local DATABASE=$5
	local THREADS=$6


	# Classify single-end or paired-end reads using centrifuge
	# Single-end
	if [ ${REVERSE} == 'NULL' ]; then
		centrifuge \
			-x ${DATABASE} \
			--seed 4 \
			-q \
			--threads ${THREADS} \
			-S ${OUT_DIR}/${SAMPLE}/${SAMPLE}.centrifuge.out \
			-U ${FORWARD} 
	else
	# Paired-end
                centrifuge \
                        -x ${DATABASE} \
                        --seed 4 \
                        -q \
                        --threads ${THREADS} \
                        -S ${OUT_DIR}/${SAMPLE}/${SAMPLE}.centrifuge.out \
                        -1 ${FORWARD} \
                        -2 ${REVERSE}
 
	fi
	
	# Generate report
	centrifuge-kreport \
    		-x ${DATABASE} \
     		${OUT_DIR}/${SAMPLE}/${SAMPLE}.centrifuge.out  \
    		> ${OUT_DIR}/${SAMPLE}/${SAMPLE}.centrifuge.out.report

	# Create text file for Krona
	cut -f 1,3 ${OUT_DIR}/${SAMPLE}/${SAMPLE}.centrifuge.out \
        > ${OUT_DIR}/${SAMPLE}/${SAMPLE}.centrifuge.out.4krona

}

export -f centrifuge_classify

# Create directory for each sample
parallel --jobs 0 "[ -d ${CENTRIFUGE_CLASSIFY}/{} ] ||  mkdir ${CENTRIFUGE_CLASSIFY}/{}" ::: ${SAMPLES[*]}

# Run centrifuge
# For single-end reads
parallel --jobs 5 "centrifuge_classify {} ${CENTRIFUGE_CLASSIFY} ${RAW_DATA_DIR}/{}${FASTQ_SUFFIX} NULL ${CENTRIFUGE_DB}/abv  10" :::  ${SAMPLES[*]}

# For paired-end reads
#parallel --jobs 5 "centrifuge_classify {} ${CENTRIFUGE_CLASSIFY} ${RAW_DATA_DIR}/{}${FORWARD_SUFFIX} ${RAW_DATA_DIR}/{}${REVERSE_SUFFIX} ${CENTRIFUGE_DB}/abv  10" :::  ${SAMPLES[*]}

find ${CENTRIFUGE_CLASSIFY}/ -type f -name "*.4krona" |sort  > ${CENTRIFUGE_CLASSIFY}/krona_files.txt
FILES=($(find ${CENTRIFUGE_CLASSIFY}/ -type f -name "*.4krona"))
basename -a -s '.4krona' ${FILES[*]} | sort | sed -i -E 's/\.centrifuge\.out//g'  > ${CENTRIFUGE_CLASSIFY}/sample_names.txt
KTAXONOMY_FILES=($(paste -d',' "${CENTRIFUGE_CLASSIFY}/krona_files.txt" "${CENTRIFUGE_CLASSIFY}/sample_names.txt"))


# Running ktImportTaxonomy to create a krona chart for samples
ktImportTaxonomy  \
    -tax ${KRONA_TAXONOMY} \
    -u  http://krona.sourceforge.net \
    -o ${CENTRIFUGE_CLASSIFY}/centrifuge_krona_report.html \
    -q 1 \
    -t 2 \
    ${KTAXONOMY_FILES[*]}



# Make a directory for Export

[ -d ${CENTRIFUGE_CLASSIFY}/Export ] || mkdir ${CENTRIFUGE_CLASSIFY}/Export

find ${CENTRIFUGE_CLASSIFY}/ -type f -name "*.report" -exec cp {} ${CENTRIFUGE_CLASSIFY}/Export \;

`
