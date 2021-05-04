#!/usr/bin/env bash
#$ -S /bin/bash
#$ -q bioinfo.q
#$ -V
#$ -cwd
#$ -N phyloflash-classify
#$ -pe shared 60
#$ -M obadbotanist@yahoo.com


# AUTHOR: Olabiyi Aderemi Obayomi
# EMAIL: obadbotanist@yahoo.com

set -e

source activate python2
export PERL5LIB="gpfs0/bioinfo/users/obayomi/miniconda3/envs/python2/lib/5.26.2/x86_64-linux-thread-multi"

PAIRED="False"
PROJECT_DIR='/gpfs0/bioinfo/projects/Amit_Gross/03.Metagenomics_stormwater_biofiltration/'
RAW_DATA_DIR="${PROJECT_DIR}/03.remove_host/fastq_files/"
OUT_DIR="${PROJECT_DIR}/04.read_based_analysis/06.phloFlash-classify/"
cd ${OUT_DIR}
SAMPLES=($(cat ${PROJECT_DIR}/metadata.tsv | awk 'NR>1{print $1}'))
SAMPLES_EDITED=($(echo ${SAMPLES[*]} | sed -E 's/\./_/g'))
FASTQ_SUFFIX='_host_removed.fastq.gz'
FORWARD_SUFFIX='_R1_host_removed.fastq.gz'
REVERSE_SUFFIX='_R2_host_removed.fastq.gz'

if [ "${PAIRED}" == "True" ]; then
	
	# Run both SPAdes and EMIRGE and produce all optional outputs
	parallel --jobs 5 "[ -d {} ] || mkdir {} && \
				cd {} && phyloFlash.pl -lib {} -almosteverything \
				-read1 ${RAW_DATA_DIR}/{}${FORWARD_SUFFIX} \
				-read2 ${RAW_DATA_DIR}/{}${REVERSE_SUFFIX}" ::: ${SAMPLES[*]}

else

	# Run both SPAdes and EMIRGE and produce all optional outputs - experimental for single-end
	parallel --jobs 5 --link  "[ -d {2} ] || mkdir {2} && \
				cd {2} && phyloFlash.pl -lib {2} -almosteverything \
				-read1 ${RAW_DATA_DIR}/{1}${FASTQ_SUFFIX}" \
				::: ${SAMPLES[*]} ::: ${SAMPLES_EDITED[*]}


fi

# Create an Export folder then move the tar.gz files there

[ -d ${OUT_DIR}/Export ] || mkdir ${OUT_DIR}/Export

find ${OUT_DIR} -type f -name "*.tar.gz" | xargs -I {} mv {} ${OUT_DIR}/Export

cd ${OUT_DIR}/Export

# Compare samples at different taxonomy levels
TAXON_LEVELS=($(seq 1 7))
TAXON_LEVEL_NAMES=('kingdom' 'phylum' 'class' 'order' 'family' 'genus' 'species')
parallel --jobs 0  --link \
	 "phyloFlash_compare.pl --allzip --task barplot,heatmap,matrix,ntu_table \
		--outfmt png --level {1} --out {2}  --displaytaxa 25 --min-ntu-count 100 \
		 > {2}.log 2>&1" \
	::: ${TAXON_LEVELS[*]} \
	::: ${TAXON_LEVEL_NAMES[*]}


