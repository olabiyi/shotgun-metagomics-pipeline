#!/usr/bin/bash
#$ -q bioinfo.q
#$ -S /bin/bash
#$ -V
#$ -cwd
#$ -N calculate-depth-per-sample
#$ -pe shared 72

# A script to calculate the depth i.e. the number of sequence that mapped to a contig
# It has two functions that can perform the calculation on a per sample or coassembly 
# level

# AUTHOR: Olabiyi Aderemi Obayomi
# EMAIL: obadbotanist@yahoo.com

set -e
#metaBAT
export PATH=/gpfs0/bioinfo/users/obayomi/metabat/:$PATH	
source /gpfs0/bioinfo/users/obayomi/miniconda3/bin/activate /gpfs0/bioinfo/users/obayomi/miniconda3/
export PERL5LIB="/gpfs0/bioinfo/users/obayomi/miniconda3/lib/5.32.0"

PROJECT_DIR='/gpfs0/bioinfo/projects/Amit_Gross/03.Metagenomics_stormwater_biofiltration/'
READS_DIR="${PROJECT_DIR}/03.remove_host/fastq_files/"

SAMPLES=($(cat ${PROJECT_DIR}/metadata.tsv | awk 'NR>1{print $1}'))
PAIRED="False"
COASSEMBLY="False" #"True"

function coassembly_map_reads_to_contigs(){
	
	local SAMPLE=$1
	local PAIRED=$2
	local CONTIGS=$3
	local ASSEMBLY_DIR=$4
	local ALIGNMENT_DIR=$5
	local READS_DIR=$6

	[ -d  ${ALIGNMENT_DIR}/${SAMPLE}/ ] || mkdir  ${ALIGNMENT_DIR}/${SAMPLE}/

		bowtie2-build \
			-f ${CONTIGS} \
        		${CONTIGS%%.fa*}.idx



	if [  ${PAIRED} == "True" ];then
		forward=$(basename ${READS_DIR}/${SAMPLE}*R1*)
		reverse=$(basename ${READS_DIR}/${SAMPLE}*R2*)
        	# Map the reads against the assembled contigs
       		 bowtie2 \
                	-x ${CONTIGS%%.fa*}.idx \
                	-1 ${READS_DIR}/${forward} \
                	-2 ${READS_DIR}/${reverse}  \
                	-S ${ALIGNMENT_DIR}/${SAMPLE}/${SAMPLE}.sam \
                	-p 10
		
	else
	 	bowtie2 \
                	-x ${CONTIGS%%.fa*}.idx \
                	-U ${READS_DIR}/${SAMPLE}* \
                	-S ${ALIGNMENT_DIR}/${SAMPLE}/${SAMPLE}.sam \
                	-p 10
	
	
	fi
	# if this error: samtools: error while loading shared libraries: libcrypto.so.1.0.0: cannot open shared object file: No such file or directory
	# Occurs while running samtools, solve it like so:
	# source activate base
	# conda install -c bioconda samtools=1.9 --force-reinstall
	samtools sort \
             	${ALIGNMENT_DIR}/${SAMPLE}/${SAMPLE}.sam \
            	-o ${ALIGNMENT_DIR}/${SAMPLE}/${SAMPLE}_sorted.bam


	# Convert the BAM files of each sample into a “depth.txt” file
	jgi_summarize_bam_contig_depths \
                     --outputDepth ${ALIGNMENT_DIR}/${SAMPLE}/${SAMPLE}.depth.txt \
                     ${ALIGNMENT_DIR}/${SAMPLE}/${SAMPLE}_sorted.bam
  
}


function sample_map_reads_to_contigs(){

        local SAMPLE=$1
        local PAIRED=$2
        local ASSEMBLY_DIR=$3
        local READS_DIR=$4

	local ALIGNMENT_DIR="${ASSEMBLY_DIR}/${SAMPLE}/Alignments"

        [ -d  ${ALIGNMENT_DIR} ] || mkdir  ${ALIGNMENT_DIR}


                CONTIGS=$(ls ${ASSEMBLY_DIR}/${SAMPLE}/*contigs.fa*)
                bowtie2-build \
                        -f ${CONTIGS} \
                        ${CONTIGS%%.fa*}.idx



        if [  ${PAIRED} == "True" ];then
                forward=$(basename ${READS_DIR}/${SAMPLE}*R1*)
                reverse=$(basename ${READS_DIR}/${SAMPLE}*R2*)
                # Map the reads against the assembled contigs
                 bowtie2 \
                        -x ${CONTIGS%%.fa*}.idx \
                        -1 ${READS_DIR}/${forward} \
                        -2 ${READS_DIR}/${reverse}  \
                        -S ${ALIGNMENT_DIR}/${SAMPLE}.sam \
                        -p 10

        else
                bowtie2 \
                        -x ${CONTIGS%%.fa*}.idx \
                        -U ${READS_DIR}/${SAMPLE}* \
                        -S ${ALIGNMENT_DIR}/${SAMPLE}.sam \
                        -p 10


        fi
        # if this error: samtools: error while loading shared libraries: libcrypto.so.1.0.0: cannot open shared object file: No such file or directory
        # Occurs while running samtools, solve it like so:
        # source activate base
        # conda install -c bioconda samtools=1.9 --force-reinstall
        samtools sort \
                ${ALIGNMENT_DIR}/${SAMPLE}.sam \
                -o ${ALIGNMENT_DIR}/${SAMPLE}_sorted.bam


        # Convert the BAM files of each sample into a “depth.txt
        jgi_summarize_bam_contig_depths \
                     --outputDepth ${ALIGNMENT_DIR}/${SAMPLE}.depth.txt \
                     ${ALIGNMENT_DIR}/${SAMPLE}_sorted.bam

}




if [ "${COASSEMBLY}" == 'True' ];then

	ASSEMBLY_DIR="${PROJECT_DIR}/05.assembly_based_analysis/01.assemble_reads/coassembly/"
	CONTIGS=$(ls "${ASSEMBLY_DIR}/*contigs.fa*")
	ALIGNMENT_DIR="${ASSEMBLY_DIR}/Alignments/"
	[ -d  ${ALIGNMENT_DIR} ] || mkdir  ${ALIGNMENT_DIR}


	# Build database index for bowtie2
	bowtie2-build \
		-f ${CONTIGS} \
       	${CONTIGS%%.fa*}.idx

	export -f  coassembly_map_reads_to_contigs
	parallel --jobs 10 "coassembly_map_reads_to_contigs {} ${COASSEMBLY} ${PAIRED} ${CONTIGS} ${ASSEMBLY_DIR} ${ALIGNMENT_DIR} ${READS_DIR}" ::: ${SAMPLES[*]} 

	SORTED_BAMS=($(find ${ALIGNMENT_DIR}/ -type f -name "*_sorted.bam"))
	# Generate a matrix of depth for all the samples together
	jgi_summarize_bam_contig_depths \
                      --outputDepth ${ALIGNMENT_DIR}/global.depth.txt \
                      ${SORTED_BAMS[*]}


else 
	ASSEMBLY_DIR="${PROJECT_DIR}/05.assembly_based_analysis/01.assemble_reads/"
	export -f  sample_map_reads_to_contigs
	parallel --jobs 10 "sample_map_reads_to_contigs {} ${PAIRED} ${ASSEMBLY_DIR} ${READS_DIR}" ::: ${SAMPLES[*]}

fi


