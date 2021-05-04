rule humann2_sample_pe:
    """
    Runs HUMAnN2 pipeline using general defaults.
    Other HUMAnN2 parameters can be specified as a quoted string in 
    PARAMS: HUMANN2: OTHER. 
    Going to do just R1 reads for now. Because of how I've split PE vs SE
    processing and naming, still will need to make a separate rule for PE. 
    """
    input:
        paired_f  = "data/{sample}/{run}/host_filtered/{sample}_R1.trimmed.host_filtered.fq.gz",
        unpaired_f = "data/{sample}/{run}/host_filtered/{sample}_U1.trimmed.host_filtered.fq.gz",
        metaphlan_in = "data/combined_analysis/{run}/humann2/joined_taxonomic_profile_max.tsv"
    output:
        genefamilies = "data/{sample}/{run}/humann2/{sample}_genefamilies.biom",
        pathcoverage = "data/{sample}/{run}/humann2/{sample}_pathcoverage.biom",
        pathabundance = "data/{sample}/{run}/humann2/{sample}_pathabundance.biom"
    params:
        metaphlan_dir = config['PARAMS']['HUMANN2']["METAPHLAN_DIR"],
        humann2_nt_db = config['PARAMS']['HUMANN2']["HUMANN2_NT_DB"],
        humann2_aa_db = config['PARAMS']['HUMANN2']["HUMANN2_AA_DB"],
        other = config['PARAMS']['HUMANN2']['OTHER']
    threads:
        8
    log:
        "logs/{run}/analysis/humann2_sample_pe_{sample}.log"
    benchmark:
        "benchmarks/{run}/analysis/humann2_sample_pe_{sample}.json"
    run:
        with tempfile.TemporaryDirectory(dir=TMP_DIR_ROOT) as temp_dir:
            shell("""
                  set +u; {HUMANN2_ENV}; set -u
                  zcat {input.paired_f} {input.unpaired_f} > {temp_dir}/input.fastq
                  humann2 --input {temp_dir}/input.fastq \
                  --output {temp_dir}/{wildcards.sample} \
                  --output-basename {wildcards.sample} \
                  --nucleotide-database {params.humann2_nt_db} \
                  --protein-database {params.humann2_aa_db} \
                  --taxonomic-profile {input.metaphlan_in} \
                  --metaphlan {params.metaphlan_dir} \
                  --o-log {log} \
                  --threads {threads} \
                  --output-format biom {params.other} 2> {log} 1>&2
                  scp {temp_dir}/{wildcards.sample}/{wildcards.sample}_genefamilies.biom {output.genefamilies}
                  scp {temp_dir}/{wildcards.sample}/{wildcards.sample}_pathcoverage.biom {output.pathcoverage}
                  scp {temp_dir}/{wildcards.sample}/{wildcards.sample}_pathabundance.biom {output.pathabundance}
                  """)

rule humann2_renorm_tables:
    """
    Renormalizes HUMAnN2 per-sample tables, per recommendation in the HUMAnN2
    website. 
    Counts-per-million (cpm) or Relative Abundance (Relabund) can be specified
    as a list in the PARAMS: HUMANN2: NORMS variable in the config file.
    """
    input:
        genefamilies = "data/{sample}/{run}/humann2/{sample}_genefamilies.biom",
        pathcoverage = "data/{sample}/{run}/humann2/{sample}_pathcoverage.biom",
        pathabundance = "data/{sample}/{run}/humann2/{sample}_pathabundance.biom"
    output:
        genefamilies = "data/{sample}/{run}/humann2/{sample}_genefamilies_{norm}.biom",
        pathcoverage = "data/{sample}/{run}/humann2/{sample}_pathcoverage_{norm}.biom",
        pathabundance = "data/{sample}/{run}/humann2/{sample}_pathabundance_{norm}.biom"
    threads:
        1
    log:
        "logs/{run}/analysis/humann2_renorm_tables_{sample}.log"
    benchmark:
        "benchmarks/{run}/analysis/humann2_renorm_tables_{sample}.json"
    run:
        shell("""
              set +u; {HUMANN2_ENV}; set -u
              humann2_renorm_table --input {input.genefamilies} \
              --output {output.genefamilies} \
              --units {wildcards.norm} 2> {log} 1>&2
              humann2_renorm_table --input {input.pathcoverage} \
              --output {output.pathcoverage} \
              --units {wildcards.norm} 2>> {log} 1>&2
              humann2_renorm_table --input {input.pathabundance} \
              --output {output.pathabundance} \
              --units {wildcards.norm} 2>> {log} 1>&2
              """)


rule humann2_combine_tables:
    """
    Combines the per-sample normalized tables into a single run-wide table. 
    Because HUMAnN2 takes a directory as input, first copies all the individual
    tables generated in this run to a temp directory and runs on that.
    """
    input:
        lambda wildcards: expand("data/{sample}/{run}/humann2/{sample}_genefamilies_{norm}.biom",
               sample=SAMPLES_PE, run=RUN, norm=wildcards.norm),
        lambda wildcards: expand("data/{sample}/{run}/humann2/{sample}_pathcoverage_{norm}.biom",
               sample=SAMPLES_PE, run=RUN, norm=wildcards.norm),
        lambda wildcards: expand("data/{sample}/{run}/humann2/{sample}_pathabundance_{norm}.biom",
               sample=SAMPLES_PE, run=RUN, norm=wildcards.norm)
    output:
        genefamilies = "data/combined_analysis/{run}/humann2/combined_genefamilies_{norm}_all.biom",
        pathcoverage = "data/combined_analysis/{run}/humann2/combined_pathcoverage_{norm}_all.biom",
        pathabundance = "data/combined_analysis/{run}/humann2/combined_pathabundance_{norm}_all.biom"
    log:
        "logs/{run}/analysis/humann2_combine_tables_{norm}.log"
    benchmark:
        "benchmarks/{run}/humann2_combine_tables_{norm}.log"
    run:
        with tempfile.TemporaryDirectory(dir='data/combined_analysis') as temp_dir:
            for file in input:
                shell("cp {0} {1}/.".format(file, temp_dir))
            shell("""
                  set +u; {HUMANN2_ENV}; set -u
                  humann2_join_tables --input {temp_dir} \
                  --output {output.genefamilies} \
                  --file_name genefamilies 2> {log} 1>&2
                  humann2_join_tables --input {temp_dir} \
                  --output {output.pathcoverage} \
                  --file_name pathcoverage 2>> {log} 1>&2
                  humann2_join_tables --input {temp_dir} \
                  --output {output.pathabundance} \
                  --file_name pathabundance 2>> {log} 1>&2
                  """)

rule humann2_remove_unmapped:
    """
    By default, HUMAnN2 includes the un-annoated reads (either unmapped in the
    first step of the pipeline or not matched in the translated alignment step)
    in the output files. In my experience, this causes relatively small
    differences in run quality (e.g. different read lengths) to have huge
    effects on the evaluated outcome, as the overall proportion of unmatched
    reads varies by run, and the compositionality of the data then causes large
    fluctuations in the count estimates of the annotated genes and pathways.
    To remove this problem, this rule renormalizes the combined tables after
    extracting unmatched read categories.
    """
    input:
        genefamilies = "data/combined_analysis/{run}/humann2/combined_genefamilies_{norm}_all.biom",
        pathcoverage = "data/combined_analysis/{run}/humann2/combined_pathcoverage_{norm}_all.biom",
        pathabundance = "data/combined_analysis/{run}/humann2/combined_pathabundance_{norm}_all.biom"
    output:
        genefamilies = "data/combined_analysis/{run}/humann2/combined_genefamilies_{norm}_mapped.biom",
        pathcoverage = "data/combined_analysis/{run}/humann2/combined_pathcoverage_{norm}_mapped.biom",
        pathabundance = "data/combined_analysis/{run}/humann2/combined_pathabundance_{norm}_mapped.biom"
    threads:
        1
    log:
        "logs/{run}/analysis/humann2_remove_un.log"
    benchmark:
        "benchmarks/{run}/analysis/humann2_remove_un.json"
    run:
        shell("""
              set +u; {HUMANN2_ENV}; set -u
              humann2_renorm_table --input {input.genefamilies} \
              --output {output.genefamilies} \
              --units {wildcards.norm} -s n 2> {log} 1>&2
              humann2_renorm_table --input {input.pathcoverage} \
              --output {output.pathcoverage} \
              --units {wildcards.norm} -s n 2>> {log} 1>&2
              humann2_renorm_table --input {input.pathabundance} \
              --output {output.pathabundance} \
              --units {wildcards.norm} -s n 2>> {log} 1>&2
              """)


rule humann2_split_stratified_tables:
    """
    Splits the grouped tables into separate grouped taxonomy-stratified and un-
    stratified tables for downstream analysis. Does this for both the combined
    tables including unmatched reads (*_all) and those excluding unmatched
    reads (*_mapped).
    The un-stratified tables should then be directly useful for downstream
    analysis in e.g. beta diversity. 
    """
    input:
        genefamilies = "data/combined_analysis/{run}/humann2/combined_genefamilies_{norm}_{mapped}.biom",
        pathcoverage = "data/combined_analysis/{run}/humann2/combined_pathcoverage_{norm}_{mapped}.biom",
        pathabundance = "data/combined_analysis/{run}/humann2/combined_pathabundance_{norm}_{mapped}.biom"
    output:
        genefamilies = "data/combined_analysis/{run}/humann2/stratified/combined_genefamilies_{norm}_{mapped}_stratified.biom",
        pathcoverage = "data/combined_analysis/{run}/humann2/stratified/combined_pathcoverage_{norm}_{mapped}_stratified.biom",
        pathabundance = "data/combined_analysis/{run}/humann2/stratified/combined_pathabundance_{norm}_{mapped}_stratified.biom",
        genefamilies_unstrat = "data/combined_analysis/{run}/humann2/stratified/combined_genefamilies_{norm}_{mapped}_unstratified.biom",
        pathcoverage_unstrat = "data/combined_analysis/{run}/humann2/stratified/combined_pathcoverage_{norm}_{mapped}_unstratified.biom",
        pathabundance_unstrat = "data/combined_analysis/{run}/humann2/stratified/combined_pathabundance_{norm}_{mapped}_unstratified.biom"
    threads:
        1
    log:
        "logs/{run}/analysis/humann2_split_stratified_tables.log"
    benchmark:
        "benchmarks/{run}/analysis/humann2_split_stratified_tables.json"
    run:
        shell("""
              set +u; {HUMANN2_ENV}; set -u
              humann2_split_stratified_table --input {input.genefamilies} \
              --output data/combined_analysis/{wildcards.run}/humann2/stratified 2> {log} 1>&2
              humann2_split_stratified_table --input {input.pathcoverage} \
              --output data/combined_analysis/{wildcards.run}/humann2/stratified 2>> {log} 1>&2
              humann2_split_stratified_table --input {input.pathabundance} \
              --output data/combined_analysis/{wildcards.run}/humann2/stratified 2>> {log} 1>&2
              """)

