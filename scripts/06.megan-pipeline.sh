#!/usr/bin/env bash
#$ -S /bin/bash
#$ -q bioinfo.q
#$ -V
#$ -cwd
#$ -N megan-pipeline
#$ -pe shared 72
#$ -M obadbotanist@yahoo.com

# MEGAN taxonomy and functional annotation pipeline

set -eo pipefail
source activate /gpfs0/bioinfo/users/obayomi/miniconda3/envs/bioinfo
export PERL5LIB='/gpfs0/bioinfo/users/obayomi/miniconda3/envs/bioinfo/lib/site_perl/5.26.2/x86_64-linux-thread-multi'
export XDG_RUNTIME_DIR=/run/user/210093
export XDG_SESSION_ID=246490
# VARIABLES
PROJECT_DIR='/gpfs0/bioinfo/projects/Amit_Gross/03.Metagenomics_stormwater_biofiltration/'
PAIRED_READS=False
RAW_DATA_DIR="${PROJECT_DIR}/03.remove_host/fastq_files/"
DATABASE='/gpfs0/bioinfo/users/obayomi/databases/nr'
METADATA="${PROJECT_DIR}/metadata.tsv"
SAMPLES=($(cat ${METADATA} | awk 'NR>1{print $1}'))
DAA_DIR="${PROJECT_DIR}/04.read_based_analysis/01.diamond-blast-sequences"
OUT_DIR="{PROJECT_DIR}/04.read_based_analysis/02.megan-classify/"
ACC2TAXA='/gpfs0/bioinfo/users/obayomi/megan/prot_acc2tax-Jul2019X1.abin'
ACC2SEED='/gpfs0/bioinfo/users/obayomi/megan/acc2seed-May2015XX.abin'
ACC2INTERPRO='/gpfs0/bioinfo/users/obayomi/megan/acc2interpro-Jul2019X.abin'
ACC2EGGNOG='/gpfs0/bioinfo/users/obayomi/megan/acc2eggnog-Jul2019X.abin'
SUFFIX='_host_removed.fastq.gz'
THREADS=10

# Taxonomy and function annotation using Diamond and megan pipeline --------------------------------------------------------------------------------------

function run_diamond_and_meganize(){

	local SAMPLE=$1
	local DATABASE=$2
	local RAW_DATA_DIR=$3
	local SUFFIX=$4
	local DAA_DIR=$5
	local THREADS=$6
	local ACC2TAXA=$7
	local ACC2SEED=$8
	local ACC2INTERPRO=$9
	local ACC2EGGNOG=${10}


		diamond blastx \
                      -d ${DATABASE} \
                       -q ${RAW_DATA_DIR}/${SAMPLE}${SUFFIX} \
                       --daa ${DAA_DIR}/${SAMPLE}.daa \
                       --threads ${THREADS} \
                       --unal 1 


		daa-meganizer --in ${DAA_DIR}/${SAMPLE}.daa  \
                		--minScore 50 \
                		--maxExpected 0.01 \
                		--minPercentIdentity 0  \
                		--topPercent 5 \
                		--minSupportPercent 0.05 \
                		--minSupport 0 \
                		--lcaAlgorithm weighted \
                		--lcaCoveragePercent 80  \
                		--readAssignmentMode readCount \
                		--acc2taxa ${ACC2TAXA} \
                		--acc2seed ${ACC2SEED} \
                		--acc2interpro2go ${ACC2INTERPRO} \
                		--acc2eggnog ${ACC2EGGNOG}


}



export -f run_diamond_and_meganize

parallel --jobs 5  " run_diamond_and_meganize {} ${DATABASE} ${RAW_DATA_DIR} ${SUFFIX} ${DAA_DIR} \
					${THREADS} ${ACC2TAXA} ${ACC2SEED} ${ACC2INTERPRO} ${ACC2EGGNOG} " ::: ${SAMPLES[*]}



# After the run completes do the following
# 1. Launch MEGAN by typing MEGAN at the command prompt
# 2. SELECT the file -> compare samples
# 3. import all your .rma/daa files that have been meganized above
# 4. explore the appropriate TAXA or SEED OR INTERPRO
# 5. important - export the files in stamp format which is very similar to the usual count tables
#  it is tab delimited with samples and taxonomic assignment levels as columns

