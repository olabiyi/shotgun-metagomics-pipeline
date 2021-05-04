#!/usr/bin/bash

#NOTE: THIS SCRIPT MUST BE RUN DIRECTLY ON THE COMMANDLINE
# A script to generate html report when running metaerg
# Often metaerg fails to generate the final html report for some reason.
# This script shoul be run when the failure occurs to ensure that the report is generated

# AUTHOR: Olabiyi Aderemi Obayomi
# EMAIL: obadbotanist@yahoo.com


set -e 

source /gpfs0/bioinfo/users/obayomi/miniconda3/bin/activate /gpfs0/bioinfo/users/obayomi/miniconda3/envs/metaerg

export PERL5LIB="/gpfs0/bioinfo/users/obayomi/miniconda3/envs/metaerg/lib/5.26.2"


COASSEMBLY="False" #"True"
PROJECT_DIR='/gpfs0/bioinfo/projects/Amit_Gross/03.Metagenomics_stormwater_biofiltration/'
ANNOTATION_DIR="${PROJECT_DIR}/05.assembly_based_analysis/04.annotate_assembly"
METADATA="${PROJECT_DIR}/metadata.tsv"
PERL='/gpfs0/bioinfo/users/obayomi/miniconda3/envs/metaerg/bin/perl'
PROGRAM='/gpfs0/bioinfo/users/obayomi/metaerg/bin/output_reports.pl'

function generate_report(){

	local SAMPLE=$1
	local ANNOTATION_DIR=$2
	local OUT_DIR="${ANNOTATION_DIR}/${SAMPLE}"
	local CONTIGS="${OUT_DIR}/${SAMPLE}.fna"
	local GFF="${OUT_DIR}/data/all.gff"


	${PERL} \
		${PROGRAM} \
		-g ${GFF} \
		-o ${OUT_DIR} \
		-f ${CONTIGS}  \
		>  ${OUT_DIR}/${SAMPLE}.log 2>&1

}


if [ "${COASSEMBLY}" == 'True' ];then

	SAMPLE='coassembly'
	
	generate_report ${SAMPLE} ${ANNOTATION_DIR}



else

	SAMPLES=($(cat ${METADATA} | awk 'NR>1{print $1}'))

	for SAMPLE in ${SAMPLES[*]}; do

		generate_report ${SAMPLE} ${ANNOTATION_DIR}
	
	done

fi
