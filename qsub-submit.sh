#!/usr/bin/env bash
#$ -S /bin/bash
#$ -q bioinfo.q
#$ -V
#$ -cwd
#$ -N submit-jobs
#$ -pe shared 72

set -e

# Activate the main conda environment
source /gpfs0/bioinfo/apps/Miniconda2/Miniconda_v4.3.21/bin/activate /gpfs0/bioinfo/apps/Miniconda2/Miniconda_v4.3.21/envs/Metagenomics
export PERL5LIB='/gpfs0/bioinfo/apps/Miniconda2/Miniconda_v4.3.21/envs/Metagenomics/lib/site_perl/5.26.2/x86_64-linux-thread-multi'

# Generate the rule graph on the commadline
# Rule graph
# snakemake -s snakefile --rulegraph | dot -Tpng > rulegraph.png
# Directed Acyclic Graph (DAG)
# snakemake -s snakefile --dag | dot -Tpng > dag.png

# Run snmakemake on the cluster
# submit a maximum 100 jobs
#  # wait for 60 seconds before declaring that a job has failed
snakemake \
	--cluster-config config/config.yaml \
	--cluster "qsub -q {cluster.queue} -S {cluster.shell} -cwd -V -N {rule} -pe shared {threads}" \
	--jobs 100 \ 
	--latency-wait 60 

