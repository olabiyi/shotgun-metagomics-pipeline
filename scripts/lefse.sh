#!/usr/bin/env bash

# The tutorial is located here: https://github.com/biobakery/biobakery/wiki/lefse
# install lefse with conda
#conda install -c biobakery lefse

# input file format
# two rows of metadata, one row of sample names and the corresponding microbial abundance table

INPUT="hmp_aerobiosis_small.txt"
BASE=${INPUT%.*}

# Format input file
# Run the following command to format the input file (hmp_small_aerobiosis.txt). This will generate a file (hmp_aerobiosis_small.in)
format_input.py ${INPUT} ${BASE}.in -c 1 -s 2 -u 3 -o 1000000

# RUN Lefse
# Run the following command, passing the file generated in the previous step as input.
# This will generate a file (hmp_aerobiosis_small.res) consisting of LEfSe analysis results.:
run_lefse.py ${BASE}.in  ${BASE}.res

# Visualization
# To visualize the results, LEfSe provides a couple of options. 
# For all the options you will need the output from run_lefse.py (in this case: hmp_aerobiosis_small.res)

# BAR plot of taxa against LDA score
# To plot the results of the LEfSe analysis generated from the previous step, run the following command.:
plot_res.py ${BASE}.res ${BASE}.png

# CLADOGRAM
# Run the following command to generate the Cladogram figure. This will use the LEfSe results file generated previously
plot_cladogram.py ${BASE}.res ${BASE}.cladogram.png --format png
