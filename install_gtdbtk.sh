#!/usr/bin/env bash

# Installing and running gtdbtk for MAGs bin taxonomy assignment
# https://ecogenomics.github.io/GTDBTk/installing/bioconda.html
# Create the GTDB-Tk environment
bash
# latest version with a new conda environment
conda create -n gtdbtk -c conda-forge -c bioconda gtdbtk

# Download and alias the GTDB-Tk reference data
# GTDB-Tk requires an environment variable named 
# GTDBTK_DATA_PATH to be set to the directory 
# containing the unarchived GTDB-Tk reference data.

#  Automatically download, and extract the GTDB-Tk reference data
source activate gtdbtk
download-db.sh


# #move into "binning" directory
#classify taxonomy of genome bins generated in the binning process
gtdbtk classify_wf --genome_dir . -x fa --out_dir gtdbtk_out --prefix your_sample_id.bin --cpus 8
