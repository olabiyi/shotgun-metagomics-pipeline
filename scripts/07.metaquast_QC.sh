#!/usr/bin/bash
#$ -q bioinfo.q
#$ -S /bin/bash
#$ -V
#$ -cwd
#$ -N QC-metaquast
#$ -pe shared 10

# A script check the quality of your Assembly using metaquast
# For multiple samples, it assumes that the assemblies of each 
# sample are store in separate directories

# AUTHOR: Olabiyi Aderemi Obayomi
# EMAIL: obadbotanist@yahoo.com

set -e

#Source my python 3 conda environment
#source /gpfs0/bioinfo/users/obayomi/miniconda3/bin/activate /gpfs0/bioinfo/users/obayomi/miniconda3/

PROJECT_DIR='/gpfs0/bioinfo/projects/Amit_Gross/03.Metagenomics_stormwater_biofiltration/'
METADATA="${PROJECT_DIR}/metadata.tsv"
SAMPLES=($(cat ${METADATA} | awk 'NR>1{print $1}'))
ASSEMBLY_DIR="${PROJECT_DIR}/05.assembly_based_analysis/01.assemble_reads/coassembly/"
QC_DIR="${PROJECT_DIR}/05.assembly_based_analysis/02.QC_assembly/coassembly/"
METAQUAST='/gpfs0/bioinfo/apps/Miniconda2/Miniconda_v4.3.21/envs/Metagenomics/lib/python2.7/site-packages/quast-5.0.2-py2.7.egg-info/scripts/metaquast.py'

[ -d ${QC_DIR} ] || mkdir ${QC_DIR}

COASSEMBLY='True' # Set to False if per-sample quality statistics is desired 

if [ ${COASSEMBLY} == 'True' ]; then

	CONTIGFILES="${ASSEMBLY_DIR}/final.contigs.fa"

else

	CONTIGFILES=''

	# Generate file names string

	for SAMPLE in ${SAMPLES[*]}; do

        	CONTIGFILES+=" ${ASSEMBLY_DIR}/${SAMPLE}/final.contigs.fa"

	done

fi

# After running megahit assess the quality of of your assembled contigs using metaquast
${METAQUAST} -t 10 -o ${QC_DIR} ${CONTIGFILES}


