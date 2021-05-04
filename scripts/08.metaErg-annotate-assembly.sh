#!/usr/bin/bash
#$ -q bioinfo.q
#$ -S /bin/bash
#$ -V
#$ -cwd
#$ -N AnnotateContigs
#$ -pe shared 72


# A script to annotate an assembly of contigs using metaErg annotation pipeline program
# AUTHOR: Olabiyi Aderemi Obayomi
# EMAIL: obadbotanist@yahoo.com

set -e 

source /gpfs0/bioinfo/users/obayomi/miniconda3/bin/activate /gpfs0/bioinfo/users/obayomi/miniconda3/envs/metaerg
export PERL5LIB="/gpfs0/bioinfo/users/obayomi/miniconda3/envs/metaerg/lib/5.26.2"

PROJECT_DIR='/gpfs0/bioinfo/projects/Amit_Gross/03.Metagenomics_stormwater_biofiltration/'
COASSEMBLY="False" #"True"
# if Coassembly is true you also need to provide a depths file so that the abundances 
# can xbe calculate per sample
DEPTHS_FILE="${ASSEMBLY_DIR}/Alignments/global.depth.txt"
 
DATABASE="/gpfs0/bioinfo/users/obayomi/metaerg/db/"
MIN_CONTIG_LEN=200

if [ "${COASSEMBLY}" == 'True' ];then

	PREFIX="coassembly"
	THERADS=72
	ASSEMBLY_DIR="${PROJECT_DIR}/05.assembly_based_analysis/01.assemble_reads/coassembly/"
	CONTIGS=$(ls "${ASSEMBLY_DIR}/*contigs.fa*")
	OUT_DIR="${PROJECT_DIR}/05.assembly_based_analysis/04.annotate_assembly/${PREFIX}"
 
	perl /gpfs0/bioinfo/users/obayomi/metaerg/bin/metaerg.pl \
		--prefix ${PREFIX} \
        	--outdir ${OUT_DIR} --locustag ${PREFIX} \
        	--depth ${DEPTHS_FILE} \
		--dbdir ${DATABASE} \
		--mincontiglen ${MIN_CONTIG_LEN} \
		--cpus ${THREADS} \
		${CONTIGS}

else

	function annotate_contigs(){
	
		local SAMPLE=$1
		local ASSEMBLY_DIR=$2
		local CONTIGS=$(ls ${ASSEMBLY_DIR}/${SAMPLE}/*contigs.fa*)
		local OUT_DIR=$3
		local MIN_CONTIG_LEN=$4
		local DATABASE=$5
		local THREADS=${6:-10}
                

		perl /gpfs0/bioinfo/users/obayomi/metaerg/bin/metaerg.pl \
        		--prefix ${SAMPLE} \
        		--outdir ${OUT_DIR}/${SAMPLE} --locustag ${SAMPLE} \
        		--dbdir ${DATABASE} \
        		--mincontiglen ${MIN_CONTIG_LEN} \
       			--cpus ${THREADS} \
        		${CONTIGS}

	}

	export  -f annotate_contigs
	SAMPLES=($(cat ${PROJECT_DIR}/metadata.tsv | awk 'NR>1{print $1}'))
	OUT_DIR="${PROJECT_DIR}/05.assembly_based_analysis/04.annotate_assembly/"
	ASSEMBLY_DIR="${PROJECT_DIR}/05.assembly_based_analysis/01.assemble_reads/"
	parallel --jobs 10 "annotate_contigs {} ${ASSEMBLY_DIR} ${OUT_DIR} ${MIN_CONTIG_LEN} ${DATABASE} ${THREADS}" ::: ${SAMPLES[*]}

fi

# IF it fails to generate the html report - metaErg fails to generate the report sometimes
# run the commands below for each sample directly on the command line
# source activate metaerg
#/gpfs0/bioinfo/users/obayomi/miniconda3/envs/metaerg/bin/perl \
#	/gpfs0/bioinfo/users/obayomi/metaerg/bin/output_reports.pl \
#	-g ${OUT_DIR}/${SAMPLE}/data/all.gff \
#	-o ${OUT_DIR}/${SAMPLE} \
#	-f ${OUT_DIR}/${SAMPLE}.fna
