#!/usr/bin env bash

# Author: Olabiyi Aderemi Obayomi
# Email: obadbotanist@yahoo.com
# To Download the host genome or any sequence from NCBI for that matter
# Visit NCBI to find the accession number(s) of the sequence(s) for download
# Assign the accession numbers to the ACCESSIONS variable below
# See the help for efetch for other interesting arguements


set -eo pipefail

# Export PERL5LIB becuase the efetch utility from NCBI is written in perl
# We set it in order to avoid conflict between different perl versions
# that may be installed

# The is the PERL5LIB of my bioinfo conda environment 
# Set it to something different when using another conda environment
export PERL5LIB="/gpfs0/bioinfo/users/obayomi/miniconda3/envs/bioinfo/lib/5.26.2"

# Activate the conda environment here bioinfo with NCBI's efectch utility installed
# By the way, this environment contains many interesting bioinfomatics tools
source activate bioinfo 

# set genomes directory where the sequence will be downloaded to
GENOMES_DIR="/gpfs0/bioinfo/projects/Amit_Gross/03.Metagenomics_stormwater_biofiltration/genomes/"
PHIX='NC_001422.1'
PARTIALS_SEQUENCES=(DQ998899.1 DQ998952.1 DQ998846.1 AY699637.1 AF148726.1 MN763569.1)
ACCESSIONS=(${PHIX} ${PARTIALS_SEQUENCES[*]})
OUT_FILE='PHIX_and_host_partial_sequences.fna'
GENOME='GCA_011114355.1'
PROJECT_ID='PRJNA540106'



function download_genome(){

        local PROJECT_ID=$1

	esearch -db bioproject -query ${PROJECT_ID} \
  	| elink -target assembly \
  	| efetch -format docsum \
  	| xtract -pattern DocumentSummary -element FtpPath_GenBank \
  	| xargs -n 1 sh -c 'wget "$0"/*fna.gz'

}


cd ${GENOMES_DIR}

# Download the Genome contained in PROJECT_ID from NCBI
download_genome ${PROJECT_ID}

# Nucleotide database
# Partial sequences - partial sequences are downloaded because the Host genome is not available
# Get the partial sequences of the host from NCBI
for ACC in ${ACCESSIONS[*]}; do efetch -db=nuccore -format=fasta -id=${ACC} >> ${OUT_FILE};done
