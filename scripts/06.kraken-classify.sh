#!/usr/bin/env bash
#$ -S /bin/bash
#$ -q bioinfo.q
#$ -V
#$ -cwd
#$ -N classify_kraken
#$ -pe shared 72
#$ -M obadbotanist@yahoo.com


# AUTHOR: Olabiyi Aderemi Obayomi
# EMAIL: obadbotanist@yahoo.com


set -e

source /gpfs0/bioinfo/apps/Miniconda2/Miniconda_v4.3.21/bin/activate /gpfs0/bioinfo/apps/Miniconda2/Miniconda_v4.3.21/envs/Metagenomics
export PERL5LIB='/gpfs0/bioinfo/apps/Miniconda2/Miniconda_v4.3.21/envs/Metagenomics/lib/site_perl/5.26.2/x86_64-linux-thread-multi'
PROJECT_DIR='/gpfs0/bioinfo/projects/Amit_Gross/03.Metagenomics_stormwater_biofiltration/'
KRAKEN_DB='/gpfs0/bioinfo/apps/Miniconda2/Miniconda_v4.3.21/envs/Metagenomics/databases/kraken2/'
KRAKEN_CLASSIFY="${PROJECT_DIR}/04.read_based_analysis/03.kraken2-classify/"
KRAKEN_REPORT="${PROJECT_DIR}/04.read_based_analysis/04.kraken2-report/"
RAW_DATA_DIR="${PROJECT_DIR}/03.remove_host/fastq_files/"
PAIRED="False"
METADATA="${PROJECT_DIR}/metadata.tsv"
FASTQ_SUFFIX='_host_removed.fastq.gz'
SAMPLES=($(cat ${METADATA} | awk 'NR>1{print $1}'))
#forward=$(basename -a ${RAW_DATA_DIR}/*R1*)
#reverse=$(basename -a ${RAW_DATA_DIR}/*R2*)
KRONA_TAXONOMY='/gpfs0/bioinfo/apps/Miniconda2/Miniconda_v4.3.21/envs/Metagenomics/databases/krona/taxonomy'
COMBINE_REPORTS='/gpfs0/bioinfo/users/obayomi/KrakenTools/combine_kreports.py'


if [ PAIRED == "True" ]; then

parallel --jobs 5  "kraken2 \
  			--db $KRAKEN_DB \
  			--quick \
  			--threads 2 \
  			--use-names \
  			--output ${KRAKEN_CLASSIFY}/{}.tab \
  			--report ${KRAKEN_REPORT}/{}_tax.txt \
  			--unclassified-out ${KRAKEN_CLASSIFY}/{}.unclassified#.fastq \
  			--classified-out ${KRAKEN_CLASSIFY}/{}.classified#.fastq \
  			--paired \
  			${RAW_DATA_DIR}/{}*R1* ${RAW_DATA_DIR}/{}*R2*" ::: ${SAMPLES[*]}
  
  

else

parallel --jobs 5  "kraken2 \
  			--db $KRAKEN_DB \
  			--quick \
  			--threads 2 \
  			--use-names \
  			--output ${KRAKEN_CLASSIFY}/{}.tab \
  			--report ${KRAKEN_REPORT}/{}_tax.txt \
  			--unclassified-out ${KRAKEN_CLASSIFY}/{}.unclassified#.fastq \
  			--classified-out ${KRAKEN_CLASSIFY}/{}.classified#.fastq \
  			${RAW_DATA_DIR}/{}${FASTQ_SUFFIX}" ::: ${SAMPLES[*]}

fi


# Merge Kraken reports
${COMBINE_REPORTS} \
    -r ${KRAKEN_REPORT}/*.txt \
    -o ${KRAKEN_REPORT}/merged.tsv \
    --sample-names ${SAMPLES[*]}


# Make Krona plot

function kraken2krona(){

	local SAMPLE=$1
	local KRAKEN_CLASSIFY=$2

	awk 'BEGIN{FS="\t"}{printf("%s\t%s\n",$2,gensub(/.*taxid ([0-9]+).*/, "\\1", "g", $3))}' \
             ${KRAKEN_CLASSIFY}/${SAMPLE}.tab  \
             >  ${KRAKEN_CLASSIFY}/${SAMPLE}.4krona

}

export -f  kraken2krona
# Preprocess the kraken out in a version compatible with krona
parallel --jobs 0  "kraken2krona {} ${KRAKEN_CLASSIFY}" ::: ${SAMPLES[*]}

find ${KRAKEN_CLASSIFY}/ -type f -name "*.4krona" |sort  > ${KRAKEN_CLASSIFY}/krona_files.txt
FILES=($(find ${KRAKEN_CLASSIFY}/ -type f -name "*.4krona"))
basename -a -s '.4krona' ${FILES[*]} | sort  > ${KRAKEN_CLASSIFY}/sample_names.txt
KTAXONOMY_FILES=($(paste -d',' "${KRAKEN_CLASSIFY}/krona_files.txt" "${KRAKEN_CLASSIFY}/sample_names.txt"))


ktImportTaxonomy  \
	-tax ${KRONA_TAXONOMY} \
	-u http://krona.sourceforge.net \
	-o ${KRAKEN_REPORT}/kraken_krona_report.html \
	-q 1 \
	-t 2 \
	${KTAXONOMY_FILES[*]}
