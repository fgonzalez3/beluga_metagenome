rule all:
    input:
        expand("results/{genera}/2_Taxonomic_Assignment/1_Taxonomic_Classification/Kraken2/{sample}/classified/csseqs_#.fq", sample=SAMPLES, genera=config["genera"]),
        expand("results/{genera}/2_Taxonomic_Assignment/1_Taxonomic_Classification/Kraken2/{sample}/unclassified/ucseqs_#.fq", sample=SAMPLES, genera=config["genera"]),
        expand("results/{genera}/2_Taxonomic_Assignment/1_Taxonomic_Classification/Kraken2/{sample}/k2_output.txt", sample=SAMPLES, genera=config["genera"]),
        expand("results/{genera}/2_Taxonomic_Assignment/1_Taxonomic_Classification/Kraken2/{sample}/k2_report.txt", sample=SAMPLES, genera=config["genera"])


        expand("results/{genera}/GTDB-TK/{sample}/gtdbtk_summary.tsv", sample=SAMPLES, genera=config["genera"]),
        expand("results/{genera}/unbinned_contigs_taxonomy_assignment/{sample}/taxonomyResult.tsv", sample=SAMPLES, genera=config["genera"]),
        expand("results/{genera}/unbinned_contigs_taxonomy_assignment/{sample}/taxonomyResult_report.report", sample=SAMPLES, genera=config["genera"]),
        expand("results/{genera}/unbinned_contigs_taxonomy_assignment/{sample}/report.html", sample=SAMPLES, genera=config["genera"]),
        expand("results/{genera}/bracken/{sample}/{sample}_bracken_level_{level}.txt", sample=SAMPLES, genera=config["genera"], level=["P", "C", "O", "F", "G", "S"]),
        expand("results/{genera}/bracken/{sample}/{sample}_bracken_combined.txt", sample=SAMPLES, genera=config["genera"]),
        expand("results/{genera}/bracken/all_combined_bracken.txt", genera=config["genera"])


checkpoint Kraken2:
    """
    Classify taxonomy for reads with Kraken2
    """
    input:
        r1 = "results/{genera}/1_pre_processing/dedup_reads/{sample}/{sample}_host_removed_dedup_R1.fastq",
        r2 = "results/{genera}/1_pre_processing/dedup_reads/{sample}/{sample}_host_removed_dedup_R2.fastq"
    output:
        classified = "results/{genera}/2_Taxonomic_Assignment/1_Taxonomic_Classification/Kraken2/{sample}/classified/csseqs_#.fq",
        unclassified = "results/{genera}/2_Taxonomic_Assignment/1_Taxonomic_Classification/Kraken2/{sample}/unclassified/ucseqs_#.fq",
        stdout = "results/{genera}/2_Taxonomic_Assignment/1_Taxonomic_Classification/Kraken2/{sample}/k2_output.txt",
        report = "results/{genera}/2_Taxonomic_Assignment/1_Taxonomic_Classification/Kraken2/{sample}/k2_report.txt"
    params:
        threads = 4,
        db = ""
    log:
        stdout = "logs/{genera}/2_Taxonomic_Assignment/1_Taxonomic_Classification/Kraken2/{sample}/Kraken2_Tax.out",
        stderr = "logs/{genera}/2_Taxonomic_Assignment/1_Taxonomic_Classification/Kraken2/{sample}/Kraken2_Tax.err"
    shell:
        """
        /vast/palmer/pi/turner/flg9/conda_envs/kraken2/kraken2 \
        --paired \
        --db {params.db} \
        --classified-out {output.classified} \
        --unclassified-out {output.unclassified} \
        --output {output.stdout} \
        --report {output.report} \
        {input.r1} {input.r2} \
        1>> {log.stdout} 2>> {log.stderr}
        """

rule MetaPhlaAn2:
    """
    """
    input:
    output:
    params:
    log:
    shell:
        """
        """

# MAG-Based Taxonomic Classification Steps

rule reconstruct_16S:
    """
    Reconstruct 16S rRNA genes and compare taxonomic composition between this output and those produced later in this pipeline
    """
    input:
    output:
    params:
    log:
    shell:
        """
        """

rule taxonomy_assignment_bins:
    """
    Assign taxonomy to bins using GTDB-Tk
    """
    input:
        bins=lambda wildcards: sorted(glob.glob(f"results/{wildcards.genera}/binning/{wildcards.sample}/MAXBIN.*.fasta"))
    output:
        "results/{genera}/GTDB-TK/{sample}/gtdbtk_summary.tsv"
    params:
        genera=config["genera"],
        bin_dir = "results/{genera}/binning/{sample}",
        outdir = "results/{genera}/GTDB-TK/{sample}", 
        threads = 4,
        prefix = "GTB-TK_"
    log:
        stdout = "logs/{genera}/gtdb-tk/{sample}/gtdb.out",
        stderr = "logs/{genera}/gtdb/{sample}/gtdb.err"
    shell:
        """
        module unload miniconda 
        source activate /home/flg9/.conda/envs/gtdbtk-2.4.1

        # Before running, activate your conda env and 'run conda env config vars set GTDBTK_DATA_PATH="/path/to/unarchived/gtdbtk/data";'

        gtdbtk classify_wf --genome_dir {params.bin_dir} --out_dir {params.outdir} \
        --cpus {params.threads} --debug --prefix {params.prefix} 1> {log.stdout} 2> {log.stderr}
        """

rule pseudo_binning_and_taxonomy_assignment:
    """
    Assign un-binned contigs to 'pseudobins' via protein clustering and assign taxonomy using MMseqs2
    """
    input:
        unbinned_contigs =lambda wildcards: sorted(glob.glob(f"results/{wildcards.genera}/binning/{wildcards.sample}/MAXBIN.*.noclass"))
    output:
        tsv = "results/{genera}/unbinned_contigs_taxonomy_assignment/{sample}/taxonomyResult.tsv",
        report = "results/{genera}/unbinned_contigs_taxonomy_assignment/{sample}/taxonomyResult_report.report",
        krona = "results/{genera}/unbinned_contigs_taxonomy_assignment/{sample}/report.html",
    params:
        genera=config["genera"],
        outdir = "results/{genera}/unbinned_contigs_taxonomy_assignment/{sample}",
        tmp = "results/{genera}/unnbinned_contigs_taxonomy_assignment/tmp",
        ranks = "genus,family,order,superkingdom"
    log:
        stdout = "logs/{genera}/unbinned_contigs_taxonomy_assignment/{sample}/mmseqs.out",
        stderr = "logs/{genera}/unbinned_contigs_taxonomy_assignment/{sample}/mmseqs.err"
    shell:
        """
        module unload miniconda 
        module load MMseqs2/14-7e284-gompi-2022b

        mkdir -p {params.outdir} {params.tmp}

        ###### Cluster ######

        # 1. Create MMseqs2 database using unbinned contigs as input
        mmseqs createdb {input.unbinned_contigs} {params.outdir}/DB

        # 2. Cluster the database and generate tsv formatted output file of clusters
        mmseqs cluster {params.outdir}/DB {params.outdir}/DB_clu {params.tmp}
        mmseqs createtsv {params.outdir}/DB {params.outdir}/DB {params.outdir}/DB_clu {params.outdir}/DB_clu/DB_clu.tsv

        # 3. Extract representative sequences from clustering
        mmseqs createsubdb {params.outdir}/DB_clu {params.outdir}/DB {params.outdir}/DB_clu_reps
      
        ###### Taxonomy Assignment ######

        # 1. Taxonomy Assignment
        mmseqs taxonomy {params.outdir}/DB_clu_reps path/to/referenceDB {params.outdir}/taxonomyResult {params.tmp}

        # 2. Reformat output into tsv
        mmseqs createtsv {params.outdir}/DB_clu_reps {params.outdir}/taxonomyResult {output.tsv}.tsv --lca-ranks {params.ranks}

        # 3. Reformat tsv into Kraken style for Bracken input
        mmseqs taxonomyreport {params.outdir}/DB_clu_reps {params.outdir}/TaxonomyResult {params.outdir}/taxonomyResult_report.report

        # 4. Create Krona interactive taxonomy report map
        mmseqs taxonomyreport {params.outdir}/DB_clu_reps {params.outdir}/taxonomyResult {output.krona} --report-mode 1

        1> {log.stdout} 2> {log.stderr}
        """

rule abundance_estimation:
    """
    Estimate relative abundance using mmseq taxonomic classification output from unbinned contigs
    """
    input:
        "results/{genera}/unbinned_contigs_taxonomy_assignment/{sample}/taxonomyResult_report.report"
    output:
        out1 = expand("results/{{genera}}/relative_abundance_calc/{{sample}}/{{sample}}_bracken_level_{level}.txt", level=["P", "C", "O", "F", "G", "S"]),
        out2 = "results/{genera}/relative_abundance_calc/{sample}/{sample}_bracken_combined.txt",
        out3 = "results/{genera}/relative_abundance_calc/all_combined_bracken.txt"
    params:
        genera=config["genera"]
        db = "/vast/palmer/pi/turner/data/db/kraken2_GTDBv220",
        threads = 10, 
        kmer_len = 35,
        read_len = 150
    shell:
        """
        module unload miniconda
        source activate /home/flg9/.conda/envs/bracken

        # Generate Bracken database file
        bracken-build -d {params.db} -t {params.threads} -k {params.kmer_len} -l {params.read_len}
        
        # Define taxonomic levels that we want abundance counts for 
        levels = ("P" "C" "O" "F" "G" "S")

        # Loop through each of the taxon levels and output a separate file for each 
        for level in "${{levels[@]}}"; do
            bracken -d {params.db} -i {input} -o {output.out1} -l $level -r {params.read_len}

        done

        # Combine all levels for the current sample
        cat {results/{wildcards.genera}/relative_abundance_calc/{wildcards.sample}/{wildcards.sample}_bracken_level_*.txt} > {output.out2}

        # Merge all of the files into one
        find results/{wildcards.genera}/relative_abundance_calc/ -name "*_bracken_level_*.txt" -exec cat {{}} + > {output.out3}

        # Old code
        #bracken -d KRAKEN2_DB -i contigs.report -o bracken_level_1.txt -l P -t 10
        #bracken -d KRAKEN2_DB -i contigs.report -o bracken_level_2.txt -l C -t 10
        #bracken -d KRAKEN2_DB -i contigs.report -o bracken_level_3.txt -l O -t 10
        #bracken -d KRAKEN2_DB -i contigs.report -o bracken_level_4.txt -l F -t 10
        #bracken -d KRAKEN2_DB -i contigs.report -o bracken_level_5.txt -l G -t 10
        #bracken -d KRAKEN2_DB -i contigs.report -o bracken_level_6.txt -l S -t 10

        # Concatenate the results into a single output file
        #cat bracken_level_*.txt > bracken.txt
        """
