#!/usr/bin/bash
#$ -q bioinfo.q
#$ -S /bin/bash
#$ -V
#$ -cwd
#$ -N classifyBins
#$ -pe shared 72

# A script to assign taxonomy names to metagenomic bins using genome taxonomy database toolkit (GTDB-TK)
# AUTHOR: Olabiyi Aderemi Obayomi
# EMAIL: obadbotanist@yahoo.com

set -e

COASSEMBLY="False" #"True"
REMOVE_UNBINNED_SAMPLES="True"
PROJECT_DIR='/gpfs0/bioinfo/projects/Amit_Gross/03.Metagenomics_stormwater_biofiltration'
BIN_DIR="${PROJECT_DIR}/05.assembly_based_analysis/05.bin_assembly"
OUT_DIR="${PROJECT_DIR}/05.assembly_based_analysis/07.annotate_bins"
# Delete samples for which no bins were detected from the samples array
delete=('59.EX-B_17052019' '60.B-SC_22052019')
METADATA="${PROJECT_DIR}/metadata.tsv"


source /gpfs0/bioinfo/users/obayomi/miniconda3/bin/activate  /gpfs0/bioinfo/users/obayomi/miniconda3/envs/gtdbtk
export PERL5LIB="/gpfs0/bioinfo/users/obayomi/miniconda3/lib/5.32.0"


# Deleting low quality bins i.e bins with completness less than 50
#BIN_LOG_FILE="/gpfs0/bioinfo/projects/Amit_Gross/03.Metagenomics_stormwater_biofiltration/05.assembly_based_analysis/logs/BinContigs.o4125066"
#BIN_DIR="/gpfs0/bioinfo/projects/Amit_Gross/03.Metagenomics_stormwater_biofiltration/05.assembly_based_analysis/05.bin_assembly/coassembly/"
#CheckM_DIR="/gpfs0/bioinfo/projects/Amit_Gross/03.Metagenomics_stormwater_biofiltration/05.assembly_based_analysis/06.check_bins/coassembly/bins/"

function delete_empty_bins(){

	# Fasta files
	#grep -E "genome" ${BIN_LOG_FILE} | \
	#  grep -v "\[" | \
	#  awk 'NR>1&&!($13>=50){print $1}' |\
	#  xargs -I {} rm -rf ${BIN_DIR}/{}.fa
  
	# Directories
	#grep -E "genome" logs/BinContigs.o4125066 | \
	# grep -v "\[" | \
	# awk 'NR>1&&!($13>=50){print $1}' | \
	# xargs -I {} rm -rf ${CheckM_DIR}/{}
	#cd -

}

if [ ${COASSEMBLY} == "True" ];then
	
	SAMPLE="coassembly"
	BIN_DIR="${BIN_DIR}/${SAMPLE}/"
	OUT_DIR="${OUT_DIR}/${SAMPLE}/"
 	
	# Classify taxonomy of genome bins generated in the binning process for coassembly
	gtdbtk classify_wf --genome_dir ${BIN_DIR} -x fa --out_dir ${OUT_DIR} --prefix ${SAMPLE} --cpus 72

else

	SAMPLES=($(cat ${METADATA} | awk 'NR>1{print $1}'))
	

	if [ ${REMOVE_UNBINNED_SAMPLES} == "True" ];then
	
 		for d in ${delete[*]};do

			SAMPLES=(${SAMPLES[*]/$d})
		
		done
	fi

	# Classify taxonomy of genome bins generated in the binning process per sample
        parallel --jobs 2 "gtdbtk classify_wf --genome_dir ${BIN_DIR}/{} -x fa --out_dir ${OUT_DIR}/{} --prefix {} --cpus 36" ::: ${SAMPLES[*]}

fi





