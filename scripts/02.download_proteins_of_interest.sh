#!/usr/bin/env bash


# Author: Olabiyi Aderemi Obayomi
# Email: obadbotanist@yahoo.com
# To Download the proteins of interest or any sequence from NCBI for that matter
# Visit NCBI to find the accession numbers for download
# Assign the accession numbers to the ACCESSIONS variable below
# See the help for efetch for other interesting arguments

set -euo pipefail

export PERL5LIB="/gpfs0/bioinfo/users/obayomi/miniconda3/envs/bioinfo/lib/5.26.2"

source activate bioinfo 
	

# Set the project directory and acession numbers to download
PROJECT_DIR="/gpfs0/bioinfo/projects/Amit_Gross/03.Metagenomics_stormwater_biofiltration"
ACCESSIONS=(ABM10554.1 AAK50270.1 AFS77481.1 AFR30904.1 WP_003538537.1 \
		AGO53779.1 AGO54289.1 AKR54998.1 QKJ72710.1 AFD26914.1 \
		NP_001169343.2 WP_013393031.1 NP_249210.1 NP_249214.1 \
		CBE68939.1 GAX62860.1)

OUT_FILE='proteins_of_interest.ffa'



cd ${PROJECT_DIR}
[ -d proteins_of_interest ] || mkdir proteins_of_interest
cd proteins_of_interest

	
# Download proteins of interest from NCBIs protein database using efetch			
for ACC in ${ACCESSIONS[*]}; do efetch -db=protein -format=fasta -id=${ACC} >> ${OUT_FILE}; done
