#!/usr/bin/bash

set -eo pipefail

# Note: Run this script only after annotating / classifying your bins with gtdbtk by running "qsub scripts/run_gtdbtk.sh"
# It is important that this script is run directly on the command line because output_reports.pl fails
# when run via qsub
source /gpfs0/bioinfo/users/obayomi/miniconda3/bin/activate /gpfs0/bioinfo/users/obayomi/miniconda3/envs/metaerg
export PERL5LIB="/gpfs0/bioinfo/users/obayomi/miniconda3/envs/metaerg/lib/5.26.2"


#Filter out contig sequences shorter than 500bp
#perl $HOME/metaerg/bin/filterContigByLength.pl contig.fasta 500

function annotate_bin(){
	
	local BIN=$1
	local BIN_DIR=$2
	local CONTIGS_ANNOTATION_DIR=$3
        local GFF="${CONTIGS_ANNOTATION_DIR}/data/all.gff"
	local BIN_CONTIGS="${BIN_DIR}/${BIN}.fa"
        local BIN_ANNOTATION_DIR=$4

	# Step1, extracting the gff-format annotations for the contigs 
	# included in "${BIN_CONTIGS}.fa" from the total metaerg dataset annotation:
	[ -d ${BIN_ANNOTATION_DIR}/${BIN} ] || mkdir ${BIN_ANNOTATION_DIR}/${BIN}
	perl /gpfs0/bioinfo/users/obayomi/metaerg/bin/fastaContig2Gff.pl \
                        -c ${BIN_CONTIGS}  -g ${GFF} \
                        > ${BIN_ANNOTATION_DIR}/${BIN}/${BIN}.gff 


	# Step 2, generating the html reports for the extracted contig subset
	perl /gpfs0/bioinfo/users/obayomi/metaerg/bin/output_reports.pl  \
		-g ${BIN_ANNOTATION_DIR}/${BIN}/${BIN}.gff \
		-f ${BIN_CONTIGS} \
		-o ${BIN_ANNOTATION_DIR}/${BIN}/ \
		> ${BIN_ANNOTATION_DIR}/${BIN}/${BIN}.log  2>&1;

}


function add_bin_id(){

	local BIN_DIR=$1
	local ANNOTATION_DIR=$2

	# Let's assume mybindir contains many nucleotide fasta files,
	#  one for each bin: Bin.1.fa", "Bin.2.fa", "Bin.3.fa"... files.
	#  The following commands will:
	# Add bin id to the fasta format of the protein coding sequence 
	# and protein coding sequence id will be in the format of "binid_geneid"
	perl /gpfs0/bioinfo/users/obayomi/metaerg/bin/add_binid2cds.pl \
		-p 'genome' \
		-d ${BIN_DIR} \
		-c ${ANNOTATION_DIR}/data/cds.faa \
		-g ${ANNOTATION_DIR}/data/all.gff > ${ANNOTATION_DIR}/data/cds_with_bin_id.faa

	# Add bin ids to master.tsv file, as the first column
	perl /gpfs0/bioinfo/users/obayomi/metaerg/bin/add_binid2master_dot_tsv.pl \
		-d ${BIN_DIR} \
		-t ${ANNOTATION_DIR}/data/master.tsv.txt > ${ANNOTATION_DIR}/data//master_with_bin_id.tsv.txt

}




REMOVE_UNBINNED_SAMPLES="False" #"True"
COASSEMBLY="True" #"False"
PROJECT_DIR='/gpfs0/bioinfo/projects/Amit_Gross/03.Metagenomics_stormwater_biofiltration'
BIN_DIR="${PROJECT_DIR}/05.assembly_based_analysis/05.bin_assembly"
CONTIGS_ANNOTATION_DIR="${PROJECT_DIR}/05.assembly_based_analysis/04.annotate_assembly"
BIN_ANNOTATION_DIR="${PROJECT_DIR}/05.assembly_based_analysis/07.annotate_bins"
METADATA="${PROJECT_DIR}/metadata.tsv"
# Samples not to include because there were no bins generated from the binning step
DELETE=('59.EX-B_17052019' '60.B-SC_22052019')

if [ ${COASSEMBLY} == "True" ]; then

	SAMPLES=('coassembly')

else
	SAMPLES=($(cat ${METADATA} | awk 'NR>1{print $1}'))

	# Delete samples for which no bins were detected from the samples array

	if [ ${REMOVE_UNBINNED_SAMPLES} == "True" ];then
	
 		for d in ${DELETE[*]};do

			SAMPLES=(${SAMPLES[*]/$d})
		
		done

	fi
fi
		
for SAMPLE in ${SAMPLES[*]}; do

     	SAMPLE_BIN_DIR="${BIN_DIR}/${SAMPLE}"
     	SAMPLE_BINS=($(ls -1 ${SAMPLE_BIN_DIR}/))
	SAMPLE_BINS=($(basename -a -s '.fa' ${SAMPLE_BINS[*]}))
     	SAMPLE_BIN_ANNOTATION_DIR="${BIN_ANNOTATION_DIR}/${SAMPLE}"

	if [ ${COASSEMBLY} == "True" ]; then

		SAMPLE_CONTIGS_ANNOTATION_DIR="${CONTIGS_ANNOTATION_DIR}/${SAMPLE}"
	else

		SAMPLE_CONTIGS_ANNOTATION_DIR="${CONTIGS_ANNOTATION_DIR}/${SAMPLE}"	
	
	fi

	for BIN in ${SAMPLE_BINS[*]}; do

		# Generate a bine specific gff file and html report annotating the bin
                annotate_bin ${BIN} ${SAMPLE_BIN_DIR} ${SAMPLE_CONTIGS_ANNOTATION_DIR} ${SAMPLE_BIN_ANNOTATION_DIR}

	done
	
	# Add bin id to the sample cds.faa file and to the master.tsv file
	add_bin_id ${SAMPLE_BIN_DIR} ${SAMPLE_CONTIGS_ANNOTATION_DIR}       
done


