#!/usr/bin/bash
#$ -q bioinfo.q
#$ -S /bin/bash
#$ -V
#$ -cwd
#$ -N assemble-reads
#$ -pe shared 72
#$ -M obadbotanist@yahoo.com

#  Ascript to 

set -e

source /gpfs0/bioinfo/users/obayomi/miniconda3/bin/activate /gpfs0/bioinfo/users/obayomi/miniconda3/
export PERL5LIB="/gpfs0/bioinfo/users/obayomi/miniconda3/lib/5.32.0"

PROJECT_DIR='/gpfs0/bioinfo/projects/Amit_Gross/03.Metagenomics_stormwater_biofiltration/'
RAW_DATA_DIR="${PROJECT_DIR}/03.remove_host/fastq_files/"

# For single-end reads
SUFFIX='_host_removed.fastq.gz'
FILES=($(find ${RAW_DATA_DIR}/ -type f -name "*${SUFFIX}"))
# For paired-end reads
FOWARD_SUFFIX='R1_host_removed.fastq.gz'
REVERSE_SUFFIX='R2_host_removed.fastq.gz'
F_FILES=($(find ${RAW_DATA_DIR}/ -type f -name "*${FOWARD_SUFFIX}"))
R_FILES=($(find ${RAW_DATA_DIR}/ -type f -name "*${REVERSE_SUFFIX}"))

ASSEMBLY_DIR="${PROJECT_DIR}/05.assembly_based_analysis/01.assemble_reads/"
COASSEMBLE='False' # Set top True for a coassembly
PAIRED='False'



function assemble_reads(){


	#USAGE:
	# PAIRED-END
	# assemble_reads ${ASSEMBLY_DIR}/sample_A -1  10 A_R1.fastq.gz A_R2.fastq.gz
        # SINGLE-END
        # assemble_reads ${ASSEMBLY_DIR}/sample_A -1  10 A.fastq.gz  

        local OUTDIR=$1
        local THREADS=$2
	local FORWARD=$3
        local REVERSE=${4:-False}


        if [ ${REVERSE} == 'False' ]; then
		
		megahit \
			--presets meta-sensitive \
			-t ${THREADS} \
			-r ${FORWARD} \
			-o ${OUTDIR}
        else

                megahit \
                        --presets meta-sensitive \
                        -t ${THREADS} \
                        -1 ${FORWARD} \
			-2 ${REVERSE} \
                        -o ${OUTDIR}

	fi


}


export -f assemble_reads

if [ ${PAIRED} == 'True' ];then


        if [ ${COASSEMBLE} == 'True' ]; then

                THREADS=72
                # Create a comma separated list of file names
                F_FILES=$(echo ${F_FILES[*]} |sed 's/ /,/g')
		R_FILES=$(echo ${R_FILES[*]} |sed 's/ /,/g')
                # Co-assembly - The output directory must not be in existence
		assemble_reads ${ASSEMBLY_DIR}/coassembly ${THREADS} ${F_FILES} ${R_FILES}


        else
                # Get a list of sample names
                SAMPLES=($(basename -a -s "${SUFFIX}" ${FILES[*]}))
		THREADS=10
                # Assembly reads per sample - The output directory must not be in existence
                parallel --jobs 10 \
			--link \
			--keep-order \
			"assemble_reads ${ASSEMBLY_DIR}/{3} ${THREADS} {1} {2}" \
			::: ${F_FILES[*]} :::  ${R_FILES[*]} ::: ${SAMPLES[*]}


	fi

else

	if [ ${COASSEMBLE} == 'True' ]; then
	
		THREADS=72	
		# Create a comma separated list of file names
		FILES=$(echo ${FILES[*]} |sed 's/ /,/g')
		# Co-assembly - The output directory must not be in existence
		assemble_reads ${ASSEMBLY_DIR}/coassembly ${THREADS} ${FILES}
	

	else
		# Get a list of sample names
		SAMPLES=($(basename -a -s "${SUFFIX}" ${FILES[*]}))
		THREADS=10
		# Assembly reads per sample - The output directory must not be in existence
		parallel --jobs 10 \
			--link  \
			--keep-order \
                        "assemble_reads ${ASSEMBLY_DIR}/{2} ${THREADS} {1}" \
			 ::: ${FILES[*]} ::: ${SAMPLES[*]}

	fi

fi
