import pandas as pd
from os import path

onsuccess:
    print("Workflow completed without any error")

onerror:
    print("An error occurred")
    shell("mail -s 'An error occurred' obadbotanist@yahoo.com < {log}")

configfile: "config/config.yaml"

sample_file = config["sample_file"]

metadata = pd.read_table(sample_file)

# Generate the rule graph on the commadline
# Rule graph
# snakemake -s snakefile --rulegraph | dot -Tpng > rulegraph.png
# Directed Acyclic Graph (DAG)
# snakemake -s snakefile --dag | dot -Tpng > dag.png

# Run snmakemake on the cluster
#snakemake --cluster "qsub -q bioinfo.q -S /bin/bash -cwd -V -N {rule} " --jobs 100

# used for generated merged kaiju tables
TAXON_LEVELS=['phylum', 'class', 'order', 'family', 'genus', 'species']

ruleorder: Diamond_blastx > Megan_classify


# Pseudo rule that defines the targets for the pipeline
rule all:
    input:
        expand(["01.raw_data/{sample}/{sample}_1.fq.gz", "01.raw_data/{sample}/{sample}_2.fq.gz"], sample=config['samples']), 
        "02.QC/pre_trim/multiqc_report.html",
        "04.QC/post_trim/multiqc_report.html",
        expand("06.Map_reads_to_Host/{sample}/{sample}.sam", sample=config['samples']),
        "08.QC/unmapped_reads/multiqc_report.html",
        config['databases']['proteins_of_interest'].replace('ffa', 'dmnd'),
        expand("09.Concat_reads/{sample}/{sample}.fastq.gz", sample=config['samples']),
        expand("10.reads_find_proteins/{sample}/{sample}_count.txt", sample=config['samples']),
        "10.Kaiju_classify/{project}.html".format(project=config['project_name']),
        expand("10.Kaiju_classify/merged.{taxon_level}.tsv", taxon_level=TAXON_LEVELS),
        expand("10.Kaiju_classify/{sample}/{sample}.phylum.tsv",sample=config['samples']),
        "11.Kraken_classify/merged.tsv",
        "11.Kraken_classify/merged.biom.tsv",
        "11.Kraken_classify/{project}.html".format(project=config['project_name']),
        expand("12.PhyloFlash_classify/{taxon_level}.ntu_table.tsv", taxon_level=TAXON_LEVELS),
        "13.Centrifuge_classify/{project}.html".format(project=config['project_name']),
        expand("14.Megan_classify/{sample}/{sample}.tkn", sample=config['samples']),
        expand("15.mOtus_classify/{sample}/{sample}.mOTU.motus", sample=config['samples']),
        "15.mOtus_classify/merged.mOTU.motus",
        expand('15.mOtus_classify/merged.{taxon_level}.motus', taxon_level=TAXON_LEVELS[:-1]),
        "16.Humann2_classify/Export/metaphlan2_merged.spf",
        "16.Humann2_classify/{project}.html".format(project=config['project_name']),
        expand("16.Humann2_classify/{sample}/{sample}_genefamilies_relab.tsv", sample=config['samples']),
        "18.QC_assembly_per_sample/report.html",
        expand("18.MetaErg_Annotation_per_sample/{sample}/index.html", sample=config['samples']),
        expand("19.assembly_find_proteins/{sample}/{sample}_count.txt", sample=config['samples']),
        expand("22.Convert_Bam2Depth_per_sample/{sample}/{sample}.depth.txt", sample=config['samples']),
        "29.MetaErg_CO_Annotate/coassembly/index.html",
        expand("31.Check_bins_per_sample/{sample}/", sample=config['samples']),
        expand("32.classify_bins_per_sample/{sample}/", sample=config['samples'])

       

# Rename files so that the file names are easy to work with 
rule rename_files:
    output:
        expand(["01.raw_data/{sample}/{sample}_1.fq.gz", "01.raw_data/{sample}/{sample}_2.fq.gz"],
                 sample=config['samples'])
    threads: 5
    run:
        for old,new in zip(metadata.Old_name,metadata.New_name):
            shell("mv {old} {new}".format(old=old, new=new))


        
# ------------------------- Quality check, trimming and contaminants removal  ----------------------------------

rule QC_pre_trim:
    input:
        forward="01.raw_data/{sample}/{sample}_1.fq.gz",
        rev="01.raw_data/{sample}/{sample}_2.fq.gz"
    output:
        forward_html="02.QC/pre_trim/{sample}/{sample}_1_fastqc.html",
        reverse_html="02.QC/pre_trim/{sample}/{sample}_2_fastqc.html"
    threads: 5
    params:
        program=config['programs_path']['fastqc'],
        out_dir=lambda w, output: path.dirname(output.forward_html)
    shell:
        "{params.program} --outdir {params.out_dir}/ "
        "--threads {threads} {input.forward} {input.rev}"

rule SummarizeQC_pre_trim:
    input:
        expand(["02.QC/pre_trim/{sample}/{sample}_1_fastqc.html",
                "02.QC/pre_trim/{sample}/{sample}_2_fastqc.html"],
                 sample=config['samples'])
    output:
        "02.QC/pre_trim/multiqc_report.html"
    threads: 1
    params:
        program=config['programs_path']['multiqc'],
        out_dir=lambda w, output: path.dirname(output[0])
    shell:
        "{params.program} --interactive -f {params.out_dir} -o {params.out_dir}"

            
adaptors=config['parameters']['trimmomatic']['adaptors']
min_length=config['parameters']['trimmomatic']['min_len']
rule Trim_reads:
    input:
        forward="01.raw_data/{sample}/{sample}_1.fq.gz",
        rev="01.raw_data/{sample}/{sample}_2.fq.gz"
    output:
        r1="03.trimmed/{sample}/{sample}_1.fq.gz",
        r2="03.trimmed/{sample}/{sample}_2.fq.gz",
        # reads where trimming entirely removed the mate
        r1_unpaired="03.trimmed/{sample}/{sample}_1.unpaired.fq.gz",
        r2_unpaired="03.trimmed/{sample}/{sample}_2.unpaired.fq.gz"
    log:
        "logs/trimmomatic/{sample}/{sample}.log"
    params:
        program=config['programs_path']['trimmomatic'],
        trimmer="ILLUMINACLIP:{adaptors}:2:30:10"
                " LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20"
                " MINLEN:{min_length}".format(adaptors=adaptors,
                                          min_length=min_length)
    threads: 20
    resources:
        mem_mb=1024
    shell:
        "{params.program} "
        "-threads {threads} "
        "{input.forward} {input.rev} "
        "{output.r1} {output.r1_unpaired} "
        "{output.r2} {output.r2_unpaired} "
        "{params.trimmer} {log}"


rule QC_post_trim:
    input:
        forward="03.trimmed/{sample}/{sample}_1.fq.gz",
        rev="03.trimmed/{sample}/{sample}_2.fq.gz"
    output:
        forward_html="04.QC/post_trim/{sample}/{sample}_1_fastqc.html",
        reverse_html="04.QC/post_trim/{sample}/{sample}_2_fastqc.html"
    threads: 5
    params:
        program=config['programs_path']['fastqc'],
        out_dir=lambda w, output: path.dirname(output.forward_html) 
    shell:
        "{params.program} --outdir {params.out_dir} "
        "--threads {threads} {input.forward} {input.rev}"

rule SummarizeQC_post_trim:
    input:
        expand(["04.QC/post_trim/{sample}/{sample}_1_fastqc.html",
                "04.QC/post_trim/{sample}/{sample}_2_fastqc.html"],
                 sample=config['samples'])
    output:
        "04.QC/post_trim/multiqc_report.html"
    params:
        program=config['programs_path']['multiqc'],
        out_dir=lambda w, output: path.dirname(output[0])
    threads: 1
    shell:
        "{params.program} --interactive -f {params.out_dir} -o {params.out_dir}"

#---- Build index and map reads to index using bwa mem
rule Build_Host_index:
    input:
        config['host_genome']
    output:
        "05.Build_Host_index/{project}_bwa_index".format(project=config['project_name'])
    params:
        program=config['programs_path']['bwa']
    threads: 5
    shell:
        "{params.program} index -p {output} {input}"


rule Map_reads_to_Host:
    input:
        index="05.Build_Host_index/{project}_bwa_index".format(project=config['project_name']),
        r1="03.trimmed/{sample}/{sample}_1.fq.gz",
        r2="03.trimmed/{sample}/{sample}_2.fq.gz"
    output:
        "06.Map_reads_to_Host/{sample}/{sample}.sam"
    log:
        "logs/bwa/{sample}.log"
    params:
        program=config['programs_path']['bwa'],
        read_group="@RG\\tID:{sample}\\tSM:{sample}"
    threads: 10
    shell:
        "{params.program} mem "
        "-t {threads} -R '{params.read_group}' "
        "{input.index} {input.r1} {input.r2} > {output} 2> {log}"

rule Filter_unmapped_reads:
    input:
        sam="06.Map_reads_to_Host/{sample}/{sample}.sam"
    output:
        sorted_bam="07.Filter_unmapped_reads/{sample}/{sample}.sorted.bam",
        flag_stat="07.Filter_unmapped_reads/{sample}/{sample}.flagstat.txt",
        index="07.Filter_unmapped_reads/{sample}/{sample}.sorted.bam.bai",
        stats="07.Filter_unmapped_reads/{sample}/{sample}.stats.txt",
        idxstats="07.Filter_unmapped_reads/{sample}/{sample}.idxstats.txt",
        fastqF="07.Filter_unmapped_reads/{sample}/{sample}.F.fastq",
        fastqR="07.Filter_unmapped_reads/{sample}/{sample}.R.fastq",
        fastqS="07.Filter_unmapped_reads/{sample}/{sample}.S.fastq"
    params:
        program=config['programs_path']['samtools'],
        flags="-f 12 -t -F 256"
    threads: 10
    shell:
        """ 
            {params.program} sort -@ {threads} -o {output.sorted_bam} {input.sam}
            {params.program} flagstat {output.sorted_bam} > {output.flag_stat}
            {params.program} index {output.sorted_bam}
            {params.program} stats --remove-dups {output.sorted_bam} > {output.stats}
            {params.program} idxstats {output.sorted_bam} > {output.idxstats}
            {params.program} fastq {params.flags} -1 {output.fastqF} -2 {output.fastqR} \
		-0 {output.fastqS}  {output.sorted_bam}
        """

rule QC_unmapped_reads:
    input:
        forward="07.Filter_unmapped_reads/{sample}/{sample}.F.fastq",
        rev="07.Filter_unmapped_reads/{sample}/{sample}.R.fastq"
    output:
        forward_html="08.QC/unmapped_reads/{sample}/{sample}_1_fastqc.html",
        reverse_html="08.QC/unmapped_reads/{sample}/{sample}_2_fastqc.html"
    params:
        program = config['programs_path']['fastqc'],
        out_dir= lambda w, output: path.dirname(output.forward_html)
    threads: 5
    shell:
        "{params.program} --outdir {params.out_dir} "
        "--threads {threads} {input.forward} {input.rev}"

rule SummarizeQC_unmapped_reads_and_alignment:
    input:
        expand(["08.QC/unmapped_reads/{sample}/{sample}_1_fastqc.html", 
               "08.QC/unmapped_reads/{sample}/{sample}_2_fastqc.html"],
               sample=config['samples'])
    output:
        "08.QC/unmapped_reads/multiqc_report.html"
    threads: 1
    params:
        program = config['programs_path']['multiqc'],
        out_dir = lambda w, output: path.dirname(output[0]),
        filter_dir= "07.Filter_unmapped_reads/" 
    shell:
        "{params.program} --interactive -f {params.filter_dir} "
        "{params.out_dir} -o {params.out_dir}"

rule Concat_reads:
    input:
        forward="07.Filter_unmapped_reads/{sample}/{sample}.F.fastq",
        rev="07.Filter_unmapped_reads/{sample}/{sample}.R.fastq"
    output:
        "09.Concat_reads/{sample}/{sample}.fastq.gz"
    params:
        temp="09.Concat_reads/{sample}/{sample}.fastq"
    threads: 2
    shell:
        "cat {input.forward} {input.rev} > {params.temp} && gzip {params.temp}"



# Make diamond database of proteins of interest
rule make_diamond_DB:
    input: 
        config['databases']['proteins_of_interest']
    output: 
        path.splitext(config['databases']['proteins_of_interest'])[0] + '.dmnd'
    params: 
        program=config['programs_path']['diamond'] 
    threads: 10
    shell:
        "{params.program} makedb --in {input} --db {output}"


#-----------------Read-based finding proteins of interest --------------------------------
rule reads_find_proteins:
    input:
        query="09.Concat_reads/{sample}/{sample}.fastq.gz",
        database=path.splitext(config['databases']['proteins_of_interest'])[0] + '.dmnd'
    output:
        tsv="10.reads_find_proteins/{sample}/{sample}_matches.tsv",
        counts="10.reads_find_proteins/{sample}/{sample}_count.txt"
    threads: 10
    params:
        program=config['programs_path']['diamond'],
        percent_identity=config['parameters']['find_protein']['percent_identity'],
        query_cover=config['parameters']['find_protein']['query_coverage'],
        evalue=config['parameters']['find_protein']['evalue'],
        outfmt=config['parameters']['find_protein']['out_format']
    threads: 10
    shell:
        """ 
        # Query protein database
        {params.program} blastx \
            --db {input.database} \
            --query {input.query} \
            --out {output.tsv} \
            --outfmt {params.outfmt} \
            --max-target-seqs 1 \
            --evalue {params.evalue} \
            --id {params.percent_identity} \
            --threads {threads} \
            --query-cover {params.query_cover} \
            --more-sensitive 
        
        # Counts the sequences per hit
        awk -F "\\t" '{{print $3}}' {output.tsv} |sort | uniq -c |sort -nr > {output.counts}
        """


# --------------------------- Read-based Taxonomy classification ---------------------------------------------

# ------------------------------ Kaiju -----------------------------------------------------

rule Kaiju_classify:
    input:
        forward="07.Filter_unmapped_reads/{sample}/{sample}.F.fastq",
        rev="07.Filter_unmapped_reads/{sample}/{sample}.R.fastq"
    output:
        "10.Kaiju_classify/{sample}/{sample}.kaiju.out"
    params:
        program=config['programs_path']['kaiju']['kaiju'],
        fmi=config['databases']['kaiju']['fmi'],
        nodes=config['databases']['kaiju']['nodes'],
        e_value=config['parameters']['kaiju']['evalue'],
        conda_activate=config['conda']['metagenomics']['env']
    threads: 10
    shell:
        """ 
        set +u
       {params.conda_activate}  
        set -u
       {params.program} -f {params.fmi} -t {params.nodes} -z {threads} -E {params.e_value} \
        -i  {input.forward} -j {input.rev} -o {output}
        """ 

rule Kaiju2krona:
    input:
        "10.Kaiju_classify/{sample}/{sample}.kaiju.out"
    output:
        "10.Kaiju_classify/{sample}/{sample}.krona.txt"
    params:
        program=config['programs_path']['kaiju']['kaiju2krona'],
        names=config['databases']['kaiju']['names'],
        nodes=config['databases']['kaiju']['nodes']
    threads: 2
    shell:
        "{params.program} -i {input} -o {output} -t {params.nodes} -n {params.names}"

rule Kaiju_ktImportText:
    input:
        files=expand("10.Kaiju_classify/{sample}/{sample}.krona.txt", sample=config['samples'])
    output:
        "10.Kaiju_classify/{project}.html".format(project=config['project_name'])
    params:
        program=config['programs_path']['kt_import_text']
    threads: 10
    run:
        files = ["{file},{sample}".format(file=file, sample=sample) 
                  for file,sample in zip(input.files,config['samples'])]
        shell("{params.program} -o {output} {files}") 

rule Kaiju2table:
    input:
        "10.Kaiju_classify/{sample}/{sample}.kaiju.out"
    output:
        "10.Kaiju_classify/{sample}/{sample}.phylum.tsv",
        "10.Kaiju_classify/{sample}/{sample}.class.tsv",
        "10.Kaiju_classify/{sample}/{sample}.order.tsv",
        "10.Kaiju_classify/{sample}/{sample}.family.tsv",
        "10.Kaiju_classify/{sample}/{sample}.genus.tsv",
        "10.Kaiju_classify/{sample}/{sample}.species.tsv"
    params:
        program=config['programs_path']['kaiju']['kaiju2table'],
        names=config['databases']['kaiju']['names'],
        nodes=config['databases']['kaiju']['nodes'],
        prefix="10.Kaiju_classify/{sample}/{sample}"
    threads: 2
    run:
        for taxon_level in TAXON_LEVELS:
            shell("{params.program} -t {params.nodes} -n {params.names} "
                  "-p -r {taxon_level} -o {params.prefix}.{taxon_level}.tsv  {input}")

rule Kaiju2_merge_tables:
    input:
        expand("10.Kaiju_classify/{sample}/{sample}.kaiju.out", sample=config['samples'])
    output:
        "10.Kaiju_classify/merged.phylum.tsv",
        "10.Kaiju_classify/merged.class.tsv",
        "10.Kaiju_classify/merged.order.tsv",
        "10.Kaiju_classify/merged.family.tsv",
        "10.Kaiju_classify/merged.genus.tsv",
        "10.Kaiju_classify/merged.species.tsv"
    params:
        program=config['programs_path']['kaiju']['kaiju2table'],
        names=config['databases']['kaiju']['names'],
        nodes=config['databases']['kaiju']['nodes'],
        prefix="10.Kaiju_classify/merged"
    threads: 5
    run:
        for taxon_level in TAXON_LEVELS:
            shell("{params.program} -t {params.nodes} -n {params.names} "
                  "-p -r {taxon_level} -o {params.prefix}.{taxon_level}.tsv  {input}")
            shell("sed -i -E 's/.+\/(.+).kaiju\.out/\1/g' {params.prefix}.{taxon_level}.tsv  && " 
                  "sed -i -E 's/file/sample/' {params.prefix}.{taxon_level}.tsv")


# ------------------------ Kraken2 -----------------------------------------------

rule Kraken_classify:
    input:
        forward="07.Filter_unmapped_reads/{sample}/{sample}.F.fastq",
        rev="07.Filter_unmapped_reads/{sample}/{sample}.R.fastq"
    output:
        tab="11.Kraken_classify/{sample}/{sample}.tab",
        tax="11.Kraken_classify/{sample}/{sample}.tax.txt",
        unclassified="11.Kraken_classify/{sample}/{sample}.unclassified#.fastq",
        classified="11.Kraken_classify/{sample}/{sample}.classified#.fastq"
    params:
        database=config['databases']['kraken'],
        program=config['programs_path']['kraken']['kraken']
    threads: 10
    shell:
        "{params.program} --db {params.database} "
        "--quick --threads {threads} --use-names "
        "--output {output.tab} --report {output.tax} "
        "--unclassified-out {output.unclassified} --unclassified-out {output.unclassified} "

rule Kraken_merge_reports:
    input:
        expand("11.Kraken_classify/{sample}/{sample}.tax.txt",  sample=config['samples'])
    output:
        "11.Kraken_classify/merged.tsv"
    params:
        program=config['programs_path']['kraken']['combine_reports'],
        samples=config['samples']
    threads: 2
    shell:
        "{params.program} -r {input} -o {output} --sample-names {params.samples}"

rule Kraken2krona:
    input:
        "11.Kraken_classify/{sample}/{sample}.tab"
    output:
        "11.Kraken_classify/{sample}/{sample}.4krona"
    threads: 1
    shell:
        """ awk 'BEGIN{{FS="\\t"}}{{printf("%s\\t%s\\n",$2,gensub(/.*taxid ([0-9]+).*/, "\\1", "g", $3))}}' {input} > {output} """

rule Kraken_ImportTaxonomy:
    input:
        files = expand("11.Kraken_classify/{sample}/{sample}.4krona", sample=config['samples'])
    output:
        "11.Kraken_classify/{project}.html".format(project=config['project_name'])
    params:
        program=config['programs_path']['kt_import_taxonomy'],
        taxonomy=config['databases']['krona_taxonomy'],
        url="http://krona.sourceforge.net",
        tax_column=2,
        query_column=1
    threads: 10
    run:
        files = ["{file},{sample}".format(file=file, sample=sample) for file,sample in zip(input.files,config['samples'])]
        shell("{params.program} -tax {params.taxonomy} -u {params.url} " 
              "-o {output} -q {params.query_column} -t {params.tax_column}  {files}") 

rule Kraken_biom:
    input:
        expand("11.Kraken_classify/{sample}/{sample}.tab", sample=config['samples'])
    output:
        gzip="11.Kraken_classify/merged.biom.gz",
        tsv="11.Kraken_classify/merged.biom.tsv"
    params:
        program_kraken=config['programs_path']['kraken']['kraken_biom'],
        program_biom=config['programs_path']['biom_convert'],
        highest_taxonomy='D',
        lowest_taxonomy='S'
    threads: 5
    shell:
        """ 
        {params.program_kraken} --gzip \
        --max {params.highest_taxonomy} --min {params.lowest_taxonomy} --output_fp {output.gzip} \
        {input}

        {params.program_biom} \
        -i {output.gzip} \
        -o {output.tsv} \
        --to-tsv --header-key taxonomy --output-metadata-id "Consensus Lineage" 
        """

# ---------------- Phyloflash ---------------------------------------

rule PhyloFlash_classify:
    input:
        forward="07.Filter_unmapped_reads/{sample}/{sample}.F.fastq",
        rev="07.Filter_unmapped_reads/{sample}/{sample}.R.fastq"
    output:
        html="12.PhyloFlash_classify/{sample}/{sample}.phyloFlash.html",
        gzip="12.PhyloFlash_classify/{sample}/{sample}.phyloFlash.tar.gz"
    log:
        "12.PhyloFlash_classify/{sample}/{sample}.phyloFlash.log"
    params:
        program=config['programs_path']['phyloflash']['main'],
        conda_activate=config['conda']['python2']['env'],
        out_dir=lambda w, output: path.dirname(output.html),
        summarise_at=config['parameters']['phyloflash']['summarize_at']
    threads: 20
    shell:
        """
        set +u
        {params.conda_activate}
        set -u
        [ -d {params.out_dir} ] || mkdir {params.out_dir} && \
        cd {params.out_dir} && {params.program} -taxlevel {params.summarise_at} -lib {wildcards.sample} \
        -almosteverything -read1 {input.forward} -read2 {input.rev}
        """


rule PhyloFlash_tables:
    input:
        expand("12.PhyloFlash_classify/{sample}/{sample}.phyloFlash.tar.gz", sample=config['samples'])
    output:
        expand("12.PhyloFlash_classify/{taxon_level}.matrix.tsv", taxon_level=TAXON_LEVELS),
        expand("12.PhyloFlash_classify/{taxon_level}.ntu_table.tsv", taxon_level=TAXON_LEVELS),
        expand("12.PhyloFlash_classify/{taxon_level}.barplot.png", taxon_level=TAXON_LEVELS),
        expand("12.PhyloFlash_classify/{taxon_level}.heatmap.png", taxon_level=TAXON_LEVELS)
    log:
        expand("12.PhyloFlash_classify/{taxon_level}.log", taxon_level=TAXON_LEVELS)
    params:
        program=config['programs_path']['phyloflash']['compare'],
        min_ntu_count=config['parameters']['phyloflash']['min_ntu_count'],
        taxa2display=config['parameters']['phyloflash']['taxa2display'],
        out_dir="12.PhyloFlash_classify"
    threads: 10
    shell:
        """
        TAXON_INDICES=($(seq 1 7))
        TAXON_LEVELS=('kingdom' 'phylum' 'class' 'order' 'family' 'genus' 'species')
        ZIP_FILES=$(echo {input} | sed -E 's/ /,/g')

        for i in ${{TAXON_INDICES[*]}}; do
            {params.program} --zip ${{ZIP_FILES}}  --task barplot,heatmap,matrix,ntu_table \
                  --level ${{i}} --out ${{TAXON_LEVELS[$i-1]}}  --displaytaxa {params.taxa2display} \
                  --min-ntu-count {params.min_ntu_count} \
                   > {params.out_dir}/${{TAXON_LEVELS[$i-1]}}.log 2>&1 
        done

        """

# ------------------------ Centrifuge ------------------------------------

rule Centrifuge_classify:
    input:
        forward="07.Filter_unmapped_reads/{sample}/{sample}.F.fastq",
        rev="07.Filter_unmapped_reads/{sample}/{sample}.R.fastq"
    output:
        tax="13.Centrifuge_classify/{sample}/{sample}.centrifuge.out",
        report="13.Centrifuge_classify/{sample}/{sample}.centrifuge.report"
    params:
        conda_activate=config['conda']['metagenomics']['env'],
        program_tax=config['programs_path']['centrifuge']['classify'],
        program_report=config['programs_path']['centrifuge']['report'],
        database=config['databases']['centrifuge'],
        random_seed=4
    threads: 10
    shell:
        """
        set +u
        {params.conda_activate}
        set -u
        {params.program_tax} -x {params.database} \
          --seed {params.random_seed} -q -S {output.tax} \
          -1 {input.forward} -2 {input.rev} 

        {params.program_report} -x {params.database} {output.tax} > {output.report} 
        """ 

rule Centrifuge2krona:
    input:
        "13.Centrifuge_classify/{sample}/{sample}.centrifuge.out"
    output:
        "13.Centrifuge_classify/{sample}/{sample}.centrifuge.4krona"
    threads: 1
    shell:
        "cut -f 1,3 {input} > {output}"


rule Centrifuge_ImportTaxonomy:
    input:
        files=expand("13.Centrifuge_classify/{sample}/{sample}.centrifuge.4krona",
                sample=config['samples'])
    output:
        "13.Centrifuge_classify/{project}.html".format(project=config['project_name'])
    params:
        program=config['programs_path']['kt_import_taxonomy'],
        taxonomy=config['databases']['krona_taxonomy'],
        url="http://krona.sourceforge.net",
        tax_column=2,
        query_column=1
    threads: 10
    run:
        files = ["{file},{sample}".format(file=file, sample=sample) for file,sample in zip(input.files,config['samples'])]
        shell("{params.program} -tax {params.taxonomy} -u {params.url} " 
              "-o {output} -q {params.query_column} -t {params.tax_column}  {files}")


# ------------------- Megan pipeline ------------------------------------

rule Diamond_blastx:
    input:
        "09.Concat_reads/{sample}/{sample}.fastq.gz"
    output:
        "14.Megan_classify/{sample}/{sample}.daa"
    params:
        database=config['databases']['diamond']['nr'],
        program=config['programs_path']['diamond']
    #threads: workflow.cores * 0.75 # set maximum threads to 75% of the avaliable cores
    threads: 10
    resources:
        # remember to set the --restart-times option to how many attempy you 
        # would like sankemake to try example --restart-times 3 for 3 attempts.
        # this will automatically increament the amount of memmory required 
        # by 500MB for every failed attempt
        mem_mb=lambda wildcards, attempt: attempt * 500 
    shell:
        "{params.program} blastx --query {input} --daa {output} --db {params.database} --threads {threads} --unal 1"

rule Megan_classify:
    input:
        # This just says that this rule is dependent on Diamond_blastx
        rules.Diamond_blastx.output
    output:
        temp("14.Megan_classify/{sample}/{sample}.tkn")
    params:
        program=config['programs_path']['megan']['meganizer'],
        ACC2TAXA=config['databases']['megan']['acc2taxa'],
        ACC2SEED=config['databases']['megan']['acc2seed'],
        ACC2INTERPRO=config['databases']['megan']['acc2interpro'],
        ACC2EGGNOG=config['databases']['megan']['acc2eggnog']
    threads: 10
    shell:
        """
        {params.program} \ 
        --in {input} \
        --minScore 50 \ 
        --maxExpected 0.01 \ 
        --minPercentIdentity 0 \
        --topPercent 5 \
        --minSupportPercent 0.05 \ 
        --minSupport 0 \
        --lcaAlgorithm weighted \
        --lcaCoveragePercent 80 \
        --readAssignmentMode readCount \
        --acc2taxa {params.ACC2TAXA} \
        --acc2seed {params.ACC2SEED} \
        --acc2interpro2go {params.ACC2INTERPRO} \
        --acc2eggnog {params.ACC2EGGNOG} && touch {output}
        """

rule mOtus_classify:
    input:
        forward="07.Filter_unmapped_reads/{sample}/{sample}.F.fastq",
        rev="07.Filter_unmapped_reads/{sample}/{sample}.R.fastq"        
    output:
        '15.mOtus_classify/{sample}/{sample}.phylum.motus',
        '15.mOtus_classify/{sample}/{sample}.class.motus',
        '15.mOtus_classify/{sample}/{sample}.order.motus',
        '15.mOtus_classify/{sample}/{sample}.family.motus',
        '15.mOtus_classify/{sample}/{sample}.genus.motus',
        '15.mOtus_classify/{sample}/{sample}.mOTU.motus'
    params:
        program=config['programs_path']['motus']
    threads: 10
    run:
        TAXON_LEVELS=('phylum', 'class', 'order', 'family', 'genus', 'mOTU')
        for taxon_level in  TAXON_LEVELS:
            shell("""
                {params.program} profile  -k {taxon_level} -c -t {threads} \
                -n {wildcards.sample} -f {input.forward} -r {input.rev} \
                -o 15.mOtus_classify/{wildcards.sample}/{wildcards.sample}.{taxon_level}.motus
                """.format(taxon_level=taxon_level))


rule mOtus_merge_phylum:
    input:
        expand('15.mOtus_classify/{sample}/{sample}.phylum.motus',
                sample=config['samples'])
    output:
        '15.mOtus_classify/merged.phylum.motus'
    params:
        program=config['programs_path']['motus']
    threads: 1
    shell:
        """
        # make a comma separated list
        FILES=$(echo {input} |sed 's/ /,/g')
        {params.program} merge -i ${{FILES}} -o {output}
       """

rule mOtus_merge_class:
    input:
        expand('15.mOtus_classify/{sample}/{sample}.class.motus',
                sample=config['samples'])
    output:
        '15.mOtus_classify/merged.class.motus'
    params:
        program=config['programs_path']['motus']
    threads: 1
    shell:
        """
        # make a comma separated list
        FILES=$(echo {input} |sed 's/ /,/g')
        {params.program} merge -i ${{FILES}} -o {output}
       """


rule mOtus_merge_order:
    input:
        expand('15.mOtus_classify/{sample}/{sample}.order.motus',
                sample=config['samples'])
    output:
        '15.mOtus_classify/merged.order.motus'
    params:
        program=config['programs_path']['motus']
    threads: 1
    shell:
        """
        # make a comma separated list
        FILES=$(echo {input} |sed 's/ /,/g')
        {params.program} merge -i ${{FILES}} -o {output}
       """


rule mOtus_merge_family:
    input:
        expand('15.mOtus_classify/{sample}/{sample}.family.motus',
                sample=config['samples'])
    output:
        '15.mOtus_classify/merged.family.motus'
    params:
        program=config['programs_path']['motus']
    threads: 1
    shell:
        """
        # make a comma separated list
        FILES=$(echo {input} |sed 's/ /,/g')
        {params.program} merge -i ${{FILES}} -o {output}
       """


rule mOtus_merge_genus:
    input:
        expand('15.mOtus_classify/{sample}/{sample}.genus.motus',
                sample=config['samples'])
    output:
        '15.mOtus_classify/merged.genus.motus'
    params:
        program=config['programs_path']['motus']
    threads: 1
    shell:
        """
        # make a comma separated list
        FILES=$(echo {input} |sed 's/ /,/g')
        {params.program} merge -i ${{FILES}} -o {output}
        """


rule mOtus_merge_mOTU:
    input:
        expand('15.mOtus_classify/{sample}/{sample}.mOTU.motus',
               sample=config['samples'])
    output:
        '15.mOtus_classify/merged.mOTU.motus'
    params:
        program=config['programs_path']['motus']
    threads: 1
    shell:
        """
        # make a comma separated list
        FILES=$(echo {input} |sed 's/ /,/g')
        {params.program} merge -i ${{FILES}} -o {output}
       """

analyses = ['ec', 'eggnog', 'genefamilies', 'go', 'infogo1000', 'kegg-orthology',
            'ko', 'level4ec', 'metacyc-rxn', 'pathabundance', 'pfam', 'rxn']

# ---------- Funtional classification and taxonomy assignment with Metaphlan2 -----
rule Humann2_classify:
    input:
        "09.Concat_reads/{sample}/{sample}.fastq.gz"
    output:
        bugs_list="16.Humann2_classify/{sample}/{sample}_humann2_temp/"
                  "{sample}_metaphlan_bugs_list.tsv",
        gene_families="16.Humann2_classify/{sample}/{sample}_genefamilies.tsv",
        pathabundance="16.Humann2_classify/{sample}/{sample}_pathabundance.tsv",
        pathcoverage="16.Humann2_classify/{sample}/{sample}_pathcoverage.tsv"
    params:
        mpa_dir=config['programs_path']['humann2']['mpa_dir'],
        UNIREF_DB=config['databases']['humann2']['uniref90'],
        CHOCOPHLAN_DB=config['databases']['humann2']['chocophlan'],
        conda_activate=config['conda']['metagenomics']['env'],
        program=config['programs_path']['humann2']['humann2']
    threads: 10
    shell: 
        """
        set +u
        {params.conda_activate}
        set -u
        {params.program} \
         --gap-fill on \
         --input-format fastq.gz \
         --minpath on \
         --nucleotide-database {params.CHOCOPHLAN_DB} \
         --protein-database {params.UNIREF_DB} \
         --threads {threads} \
         --output 16.Humann2_classify/{wildcards.sample} \
         --input {input}
        """

rule Metaphlan_merge:
    input:
        expand("16.Humann2_classify/{sample}/{sample}_humann2_temp/"
               "{sample}_metaphlan_bugs_list.tsv",
                sample=config['samples'])
    output:
        tsv="16.Humann2_classify/Export/metaphlan2_merged.tsv",
        # STAMP format required for further downstream analyses
        spf="16.Humann2_classify/Export/metaphlan2_merged.spf"
    params:
        conda_activate=config['conda']['metagenomics']['env'],
        program=config['programs_path']['metaphlan']['merge_metaphlan'],
        metaplan2stamp=config['programs_path']['metaphlan']['metaphlan2stamp'],
        PERL5LIB=config['conda']['metagenomics']['perl5lib']
    threads: 5
    shell:
        """
        set +u
        {params.conda_activate}
        {params.PERL5LIB}
        set -u
        [ -d  16.Humann2_classify/Export/ ] || mkdir 16.Humann2_classify/Export/
        {params.program} {input}  > {output.tsv}
        {params.metaplan2stamp} {output.tsv}  > {output.spf}
        """

rule Metaphlan2krona:
    input:
        "16.Humann2_classify/{sample}/{sample}_humann2_temp/"
        "{sample}_metaphlan_bugs_list.tsv"
    output:
        "16.Humann2_classify/{sample}/{sample}_humann2_temp/{sample}.krona.txt"
    params:
        conda_activate=config['conda']['metagenomics']['env'],
        program=config['programs_path']['metaphlan']['metaphlan2krona']
    threads: 1
    shell:
        """
        set +u
        {params.conda_activate}
        set -u
        {params.program} -p {input} -k {output}
        """

rule Metaphlan_ktextImport:
    input:
        files=expand("16.Humann2_classify/{sample}/{sample}_humann2_temp/{sample}.krona.txt",
                      sample=config['samples'])
    output:
        "16.Humann2_classify/{project}.html".format(project=config['project_name'])
    params:
        program=config['programs_path']['kt_import_text']
    threads: 1
    run:
        files = ["{file},{sample}".format(file=file, sample=sample)
                  for file,sample in zip(input.files,config['samples'])]
        shell("{params.program} -o {output} {files}")

 
rule Humann2_normalize:
    input:
        gene_families="16.Humann2_classify/{sample}/{sample}_genefamilies.tsv",
        pathabundance="16.Humann2_classify/{sample}/{sample}_pathabundance.tsv" 
    output:
        gene_families="16.Humann2_classify/{sample}/{sample}_genefamilies_relab.tsv",
        pathabundance="16.Humann2_classify/{sample}/{sample}_pathabundance_relab.tsv" 
    params:
        conda_activate=config['conda']['metagenomics']['env'],
        program=config['programs_path']['humann2']['renorm'] 
    threads: 1
    shell:
        """
        set +u
        {params.conda_activate}
        set -u
        {params.program} \
         -i {input.gene_families} \
         -o {output.gene_families} \ 
         --units relab

        {params.program} \
         -i {input.pathabundance} \
         -o {output.pathabundance} \ 
         --units relab
        """

#rule Humann2_export:
#    input:
#    output:
#    params:
#    shell:

# ------------------------ Assembly based analysis ---------------------------------

# Assembly and annotation per sample

rule Assembly_per_sample:
    input:
        forward="07.Filter_unmapped_reads/{sample}/{sample}.F.fastq",
        rev="07.Filter_unmapped_reads/{sample}/{sample}.R.fastq"        
    output:
        contigs="17.Assembly_per_sample/{sample}/{sample}.contigs.fa"
    log:
        "17.Assembly_per_sample/{sample}/{sample}.log"
    params:
        out_dir=lambda w, output: path.dirname(output.contigs),
        program=config['programs_path']['megahit']
    threads: 10
    shell:
        "{params.program} "
        "--presets meta-sensitive "
        "--continue "
        "--num-cpu-threads {threads} "
        "--out-dir {params.out_dir} "
        "--out-prefix {wildcards.sample} "
        "-1 {input.forward} "
        "-2 {input.rev} "


rule QC_assembly_per_sample:
    input:
        expand("17.Assembly_per_sample/{sample}/{sample}.contigs.fa",
               sample=config['samples'])
    output:
        "18.QC_assembly_per_sample/report.html",
        "18.QC_assembly_per_sample/icarus.html"
    log:
        "18.QC_assembly_per_sample/metaquast.log"
    params:
        out_dir=lambda w, output: path.dirname(output[0]),
        program=config['programs_path']['metaquast'],
        conda_activate=config['conda']['metagenomics']['env']
    threads: 10
    shell:
        """
        set +u
        {params.conda_activate}
        set -u
        {params.program} --output-dir {params.out_dir} {input}
        """

rule MetaErg_Annotation_per_sample:
    input:
        "17.Assembly_per_sample/{sample}/{sample}.contigs.fa"
    output:
        contigs="18.MetaErg_Annotation_per_sample/{sample}/{sample}.fna",
        zip_file="18.MetaErg_Annotation_per_sample/{sample}.tar.gz",
        html="18.MetaErg_Annotation_per_sample/{sample}/index.html",
        coding_sequences="18.MetaErg_Annotation_per_sample/{sample}/data/cds.faa"
    log:
        "18.MetaErg_Annotation_per_sample/{sample}/{sample}.log"
    threads: 10
    params:
        program="{bin}/metaerg.pl".format(bin=config['programs_path']['metaerg']),
        PERL5LIB=config['conda']['metaerg']['perl5lib'],
        conda_activate=config['conda']['metaerg']['env'],
        out_dir= lambda w, output: path.dirname(output.contigs),
        db_dir=config['databases']['metaerg'],
        min_contig_len=config['parameters']['metaerg']['min_contig_len']
    shell:
        """
        set +u
        {params.conda_activate}
        {params.PERL5LIB}
        set -u
        {params.program} \
         --dbdir {params.db_dir} \
         --mincontiglen {params.min_contig_len} \
         --cpus {threads} \
         --outdir {params.out_dir} \
         --force \
         --prefix {wildcards.sample} \
         --locustag {wildcards.sample} \
         {input} 2> {log}
        """

# -----------------Assembly-based finding proteins of interest -------------------------------
rule assembly_find_proteins:
    input:
        # coding sequences
        query="18.MetaErg_Annotation_per_sample/{sample}/data/cds.faa",
        database=path.splitext(config['databases']['proteins_of_interest'])[0] + '.dmnd'
    output:
        tsv="19.assembly_find_proteins/{sample}/{sample}_matches.tsv",
        counts="19.assembly_find_proteins/{sample}/{sample}_count.txt"
    threads: 10
    params:
        program=config['programs_path']['diamond'],
        percent_identity=config['parameters']['find_protein']['percent_identity'],
        query_cover=config['parameters']['find_protein']['query_coverage'],
        evalue=config['parameters']['find_protein']['evalue'],
        outfmt=config['parameters']['find_protein']['out_format']
    shell:
        """ 
        # Query protein database
        {params.program} blastx \
            --db {input.database} \
            --query {input.query} \
            --out {output.tsv} \
            --outfmt {params.outfmt} \
            --max-target-seqs 1 \
            --evalue {params.evalue} \
            --id {params.percent_identity} \
            --threads {threads} \
            --query-cover {params.query_cover} \
            --more-sensitive 
        
        # Counts the sequences per hit
        awk -F "\\t" '{{print $3}}' {output.tsv} |sort | uniq -c |sort -nr > {output.counts}
        """


# ----------------- Calculate depth/coverage per sample

# Bowtie build index per sample
rule Build_index_per_sample:
    input:
        "17.Assembly_per_sample/{sample}/{sample}.contigs.fa"
    output:
        "19.Build_index_per_sample/{sample}/{sample}.index" 
    log:
        "19.Build_index_per_sample/{sample}/{sample}.log"
    params:
        program=config['programs_path']['bowtie']['build'],
        conda_activate=config['conda']['metagenomics']['env']
    threads: 5
    shell:
        """
        set +u 
        {params.conda_activate}
        set -u
         {params.program} {input} {output} 2> {log}
        """


rule Map_reads2index_per_sample:
    input:
        index="19.Build_index_per_sample/{sample}/{sample}.index",
        fastqF="07.Filter_unmapped_reads/{sample}/{sample}.F.fastq",
        fastqR="07.Filter_unmapped_reads/{sample}/{sample}.R.fastq",
        fastqS="07.Filter_unmapped_reads/{sample}/{sample}.S.fastq"    
    output:
        stats="20.Map_reads2index_per_sample/{sample}/{sample}.stats",
        sam="20.Map_reads2index_per_sample/{sample}/{sample}.sam"
    log:
        "20.Map_reads2index_per_sample/{sample}/{sample}.log"
    params:
        program=config['programs_path']['bowtie']['bowtie'],
        conda_activate=config['conda']['metagenomics']['env']
    threads: 20        
    shell:
        """
        set +u
        {params.conda_activate}
        set -u
        {params.program} \
         -p {threads} \
         --rg-id {wildcards.sample} \
         --rg SM:{wildcards.sample} \
         -x {input.index} \
         -1 {input.fastqF} \
         -2 {input.fastqR} \
         -U {input.fastqS} \
         --met-file {output.stats} \
         -S {output.sam} \
          2> {log}
        """

rule Sort_sam_per_sample:
    input:
        sam="20.Map_reads2index_per_sample/{sample}/{sample}.sam"
    output:
        sorted_bam="21.Sort_sam_per_sample/{sample}/{sample}.sorted.bam",
        index="21.Sort_sam_per_sample/{sample}/{sample}.sorted.bam.bai"
    log:
        "21.Sort_sam_per_sample/{sample}/{sample}.sorted.log"
    threads: 10
    params:
        program=config['programs_path']['samtools']
    shell:
        """
            {params.program} sort -@ {threads} -o {output.sorted_bam} {input.sam}
            {params.program} index {output.sorted_bam}
        """

rule Convert_Bam2Depth_per_sample:
    input:
        "21.Sort_sam_per_sample/{sample}/{sample}.sorted.bam"
    output:
        "22.Convert_Bam2Depth_per_sample/{sample}/{sample}.depth.txt"
    log:
        "22.Convert_Bam2Depth_per_sample/{sample}/{sample}.log"
    threads: 10
    params:
        program=config['programs_path']['jgi_summarize'],
        conda_activate=config['conda']['base']['env'],
        path=config['programs_path']['metabat']['path']
    shell:
        """
        set +u
        {params.conda_activate}
        set -u
        # export metabat path
        {params.path}
        {params.program} --outputDepth {output} {input}
        """


# ----------- Co-Assembly and annotation----------------------------------------

rule Co_Assembly:
    input:
        forward=expand("07.Filter_unmapped_reads/{sample}/{sample}.F.fastq", 
                       sample=config['samples']),
        rev=expand("07.Filter_unmapped_reads/{sample}/{sample}.R.fastq", 
                    sample=config['samples']),
        single=expand("07.Filter_unmapped_reads/{sample}/{sample}.R.fastq",
                       sample=config['samples'])        
    output:
        contigs="23.Co_Assembly/coassembly/coassembly.contigs.fa"
    log:
        "23.Co_Assembly/coassembly/coassembly.log"
    params:
        out_dir= lambda w, output: path.dirname(output.contigs),
        program=config['programs_path']['megahit'],
        prefix="coassembly"
    threads: 40
    shell:
        """
        # Substitute spaces with commas
        FORWARD=$(echo {input.forward} | sed -E 's/ /,/g')
        REVERSE=$(echo {input.rev} | sed -E 's/ /,/g')
        SINGLE=$(echo {input.single} | sed -E 's/ /,/g')

        {params.program} \
        --presets meta-sensitive \
        --continue \
        --num-cpu-threads {threads} \
        --out-dir {params.out_dir} \
        --out-prefix {params.prefix} \
        -1 ${{FORWARD}} \
        -2 ${{REVERSE}} \
        -r ${{SINGLE}}

        """

rule QC_Coassembly:
    input:
        "23.Co_Assembly/coassembly/coassembly.contigs.fa"
    output:
        "24.QC_Coassembly/coassembly/report.html",
        "24.QC_Coassembly/coassembly/icarus.html"
    log:
        "24.QC_Coassembly/coassembly/metaquast.log"
    params:
        out_dir= lambda w, output: path.dirname(output[0]),
        program=config['programs_path']['metaquast'],
        conda_activate=config['conda']['metagenomics']['env']
    shell:
        """
        set +u
        {params.conda_activate}
        set -u
        {params.program} --output-dir {params.out_dir} {input}
        """

#------- Calculate depth/coverage per sample of the Coassembly

# Bowtie build index for the coassembly
rule Build_index_of_CoAssembly:
    input:
        "23.Co_Assembly/coassembly/coassembly.contigs.fa"
    output:
        "25.Build_index_of_CoAssembly/coassembly/coassembly.index" 
    log:
        "25.Build_index_of_CoAssembly/coassembly/coassembly.log"
    threads: 10
    params:
        program=config['programs_path']['bowtie']['build'],
        conda_activate=config['conda']['metagenomics']['env']
    shell:
        """
        set +u
        {params.conda_activate}
        set -u
        {params.program} {input} {output}
        """


rule Map_reads2index_of_CoAssembly:
    input:
        index="25.Build_index_of_CoAssembly/coassembly/coassembly.index",
        fastqF="07.Filter_unmapped_reads/{sample}/{sample}.F.fastq",
        fastqR="07.Filter_unmapped_reads/{sample}/{sample}.R.fastq",
        fastqS="07.Filter_unmapped_reads/{sample}/{sample}.S.fastq"    
    output:
        stats="26.Map_reads2index_of_CoAssembly/{sample}/{sample}.stats",
        sam="26.Map_reads2index_of_CoAssembly/{sample}/{sample}.sam"
    log:
        "26.Map_reads2index_of_CoAssembly/{sample}/{sample}.log"
    params:
        program=config['programs_path']['bowtie']['bowtie'],
        conda_activate=config['conda']['metagenomics']['env']
    threads: 20        
    shell:
        """
        set +u
        {params.conda_activate}
        set -u
        {params.program} \
         -p {threads} \
         --rg-id {wildcards.sample} \
         --rg SM:{wildcards.sample} \
         -x {input.index} \
         -1 {input.fastqF} \
         -2 {input.fastqR} \
         -U {input.fastqS} \
         --met-file {output.stats} \
         -S {output.sam} \
          2> {log}
        """

rule Sort_sam_of_CoAssembly:
    input:
        sam="26.Map_reads2index_of_CoAssembly/{sample}/{sample}.sam"
    output:
        sorted_bam="27.Sort_sam_of_CoAssembly/{sample}/{sample}.sorted.bam",
        index="27.Sort_sam_of_CoAssembly/{sample}/{sample}.sorted.bam.bai"
    log:
        "27.Sort_sam_of_CoAssembly/{sample}/{sample}.sorted.log"
    threads: 10
    params:
        program=config['programs_path']['samtools']
    shell:
        """
            {params.program} sort -@ {threads} -o {output.sorted_bam} {input.sam}
            {params.program} index {output.sorted_bam}
        """

rule Convert_Bam2Depth_of_CoAssembly:
    input:
        expand("27.Sort_sam_of_CoAssembly/{sample}/{sample}.sorted.bam",
               sample=config['samples'])
    output:
        "28.Convert_Bam2Depth_of_CoAssembly/coassembly/coassembly.depth.txt"
    log:
        "28.Convert_Bam2Depth_of_CoAssembly/coassembly/coassembly.log"
    params:
        program=config['programs_path']['jgi_summarize'],
        conda_activate=config['conda']['base']['env'],
        path=config['programs_path']['metabat']['path']
    threads: 10
    shell:
        """
        set +u
        {params.conda_activate}
        set -u
        # export metabat path
        {params.path}
        {params.program} --outputDepth {output} {input}
        """ 

rule MetaErg_CO_Annotate:
    input:
        contigs="23.Co_Assembly/coassembly/coassembly.contigs.fa",
        depth="28.Convert_Bam2Depth_of_CoAssembly/coassembly/coassembly.depth.txt"
    output:
        contigs="29.MetaErg_CO_Annotate/coassembly/coassembly.fna",
        zip_file="29.MetaErg_CO_Annotate/coassembly.tar.gz",
        html="29.MetaErg_CO_Annotate/coassembly/index.html"
    log:
        "29.MetaErg_CO_Annotate/coassembly/coassembly.log"
    params:
        program="{bin}/metaerg.pl".format(bin=config['programs_path']['metaerg']),
        PERL5LIB=config['conda']['metaerg']['perl5lib'],
        conda_activate=config['conda']['metaerg']['env'],
        out_dir= lambda w, output: path.dirname(output.contigs),
        db_dir=config['databases']['metaerg'],
        prefix="coassembly",
        min_contig_len=config['parameters']['metaerg']['min_contig_len']
    threads: 40
    shell:
        """
        set +u
        {params.conda_activate}
        {params.PERL5LIB}
        set -u
        {params.program} \
         --dbdir {params.db_dir} \
         --mincontiglen {params.min_contig_len} \
         --cpus {threads} \
         --outdir {params.out_dir} \
         --force \
         --prefix {params.prefix} \
         --locustag {params.prefix} \
         --depth {input.depth} \
         {input.contigs} 2> {log}
        """

## ---------------- Binning, refining and Qality Check of bins --------------------------------

# Binning per sample using metabat
rule Bin_assembly_per_sample:
    input:
        # Bin only contigs with length >= 200bp
        contigs="18.MetaErg_Annotation_per_sample/{sample}/{sample}.fna",
        depth="22.Convert_Bam2Depth_per_sample/{sample}/{sample}.depth.txt"
    output:
        directory("30.Bin_assembly_per_sample/{sample}/")
    params:
        program=config['programs_path']['metabat']['metabat'],
        conda_activate=config['conda']['base']['env'],
        #out_dir="30.Bin_assembly_per_sample/{sample}/"
    threads: 10
    shell:
        """
        set +u
        {params.conda_activate}
        set -u
        {params.program} -i {input.contigs} -a {input.depth} -o {output}/genome 
        """


#rule Refine_bins_per_sample:
#    input:
#    output:
#    log:
#    params:
#    shell:

rule Check_bins_per_sample:
    input:
        "30.Bin_assembly_per_sample/{sample}/"
    output:
        directory("31.Check_bins_per_sample/{sample}/")
    params:
        program=config['programs_path']['checkm'],
        conda_activate=config['conda']['base']['env'],
        fasta_extension='fa'
    threads: 10     
    shell:
        """
           set +u
           {params.conda_activate}
           set -u
           {params.program} \
            -t {threads} \
            -x {params.fasta_extension} \
            {input} {output} 
        """

# Assign taxonomy to MAGS using GTDB-TK
rule classify_bins_per_sample:
    input:
        "30.Bin_assembly_per_sample/{sample}/"
    output:
        directory("32.classify_bins_per_sample/{sample}/")
    threads: 72
    log:
        "logs/classify_bins_per_sample/{sample}.log"
    params:
        fasta_extension='fa',
        threads=36,
        conda_activate=config['conda']['gtdb-tk']['env']
    shell:
        """
        set +u
        {params.conda_activate}
        set -u
        gtdbtk classify_wf \
        --genome_dir {input} \
        -x {params.fasta_extension} \
        --out_dir {output} \
        --prefix {wildcards.sample} \
        --cpus {params.threads} 2> {log}
        """

localrules: metaerg_annotate_bins_per_sample

# This step should not be run on the cluste but run locally
rule metaerg_annotate_bins_per_sample:
    input:
        "30.Bin_assembly_per_sample/{sample}/"
    output:
        directory("33.metaerg_annotate_bins_per_sample/{sample}/")
    threads: 10
    log: 
        "logs/metaerg_annotate_bins_per_sample/{sample}.log"
    params:
        metaerg_bin=config['programs_path']['metaerg'],
        conda_activate=config['conda']['metaerg']['env'],
        PERL5LIB=config['conda']['metaerg']['perl5lib']
    shell:
        """
	function annotate_bin(){

        	local BIN=$1
        	local BIN_DIR=$2
        	local CONTIGS_ANNOTATION_DIR=$3
        	local GFF="${{CONTIGS_ANNOTATION_DIR}}/data/all.gff"
        	local BIN_CONTIGS="${{BIN_DIR}}/${{BIN}}.fa"
        	local BIN_ANNOTATION_DIR=$4

        	# Step1, extracting the gff-format annotations for the contigs
        	# included in "${{BIN_CONTIGS}}.fa" from the total metaerg dataset annotation:
        	[ -d ${{BIN_ANNOTATION_DIR}}/${{BIN}} ] || mkdir ${{BIN_ANNOTATION_DIR}}/${{BIN}}
        	{params.metaerg_bin}/fastaContig2Gff.pl \
                        -c ${{BIN_CONTIGS}}  -g ${{GFF}} \
                        > ${{BIN_ANNOTATION_DIR}}/${{BIN}}/${{BIN}}.gff


        	# Step 2, generating the html reports for the extracted contig subset
        	{params.metaerg_bin/output_reports.pl  \
                	-g ${{BIN_ANNOTATION_DIR}}/${{BIN}}/${{BIN}}.gff \
                	-f ${{BIN_CONTIGS}} \
                	-o ${{BIN_ANNOTATION_DIR}}/${{BIN}}/ \
                	> ${{BIN_ANNOTATION_DIR}}/${{BIN}}/${{BIN}}.log  2>&1;

	}


	function add_bin_id(){

        	local BIN_DIR=$1
       		local ANNOTATION_DIR=$2

        	# Let's assume mybindir contains many nucleotide fasta files,
        	#  one for each bin: Bin.1.fa", "Bin.2.fa", "Bin.3.fa"... files.
        	#  The following commands will:
        	# Add bin id to the fasta format of the protein coding sequence
        	# and protein coding sequence id will be in the format of "binid_geneid"
        	{params.metaerg_bin}/add_binid2cds.pl \
                	-p 'genome' \
                	-d ${{BIN_DIR}} \
                	-c ${{ANNOTATION_DIR}}/data/cds.faa \
                	-g ${{ANNOTATION_DIR}}/data/all.gff \
			> ${{ANNOTATION_DIR}}/data/cds_with_bin_id.faa

        	# Add bin ids to master.tsv file, as the first column
        	{params.metaerg_bin}/add_binid2master_dot_tsv.pl \
                	-d ${{BIN_DIR}} \
                	-t ${{ANNOTATION_DIR}}/data/master.tsv.txt \
			> ${{ANNOTATION_DIR}}/data//master_with_bin_id.tsv.txt

	}

	for BIN in ${SAMPLE_BINS[*]}; do

                # Generate a bin specific gff file and html report annotating the bin
                annotate_bin ${{BIN}} ${{SAMPLE_BIN_DIR}} ${{SAMPLE_CONTIGS_ANNOTATION_DIR}} ${{SAMPLE_BIN_ANNOTATION_DIR}}

        done

        # Add bin id to the sample cds.faa file and to the master.tsv file
        add_bin_id ${{SAMPLE_BIN_DIR}} ${{SAMPLE_CONTIGS_ANNOTATION_DIR}}


        """
