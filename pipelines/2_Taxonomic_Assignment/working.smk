rule all:
    input:
        expand("results/{genera}/bracken/{sample}/{sample}_bracken_level_{level}.txt", sample=SAMPLES, genera=config["genera"], level=["P", "C", "O", "F", "G", "S"])

rule bracken_build: # bracken k-mer files are already built for this db so no need to run this 
    """
    Generate bracken database file necessary for abundance estimation
    """
    input:
        db = "/vast/palmer/pi/turner/data/db/kraken2_GTDBv220"
    output:
        check = "results/{genera}/2_Taxonomic_Assignment/1_Taxonomic_Classification/Bracken/check.txt"
    params:
        threads = 4,
        kmer = 35,
        readlen = 150,
        executables = "/vast/palmer/pi/turner/flg9/conda_envs/kraken2"
    log:
        stdout = "logs/{genera}/2_Taxonomic_Assignment/1_Taxonomic_Classification/Bracken/Bracken_build.out",
        stderr = "logs/{genera}/2_Taxonomic_Assignment/1_Taxonomic_Classification/Bracken/Bracken_build.err"
    shell:
        """
        module unload miniconda
        source activate /home/flg9/.conda/envs/bracken

        bracken-build \
        -d {input.db} \
        -t {params.threads} \
        -k {params.kmer} \
        -l {params.readlen} \
        -x {params.executables} \
        1>> {log.stdout} 2>> {log.stderr}

        # Write to an empty file once this process wraps up
        touch {output}
        """

rule MetaPhlaAn2:
    """
    Run taxonomic assignment on short reads w/ MetaPhlan
    """
    input:
        r1 = "results/{genera}/1_pre_processing/dedup_reads/{sample}/{sample}_host_removed_dedup_R1.fastq",
        r2 = "results/{genera}/1_pre_processing/dedup_reads/{sample}/{sample}_host_removed_dedup_R2.fastq"
    output:
        profile = "profiled_metagenome.txt",
        mapout = "metagenome.bowtie2.bz2",
        threads = 5,
        db = "/vast/palmer/pi/turner/data/db/MetaPhlAn"
    params:
        type = "fastq"
    log:
        stdout = "logs/{genera}/2_Taxonomic_Assignment/1_Taxonomic_Classification/MetaPhlan/{sample}/MetaPhlan_Tax.out",
        stderr = "logs/{genera}/2_Taxonomic_Assignment/1_Taxonomic_Classification/MetaPhlan/{sample}/MetaPhlan_Tax.err"
    shell:
        """
        module unload miniconda
        source activate /vast/palmer/pi/turner/flg9/conda_envs/metaphlan
        
        metaphlan \
        {input.r1},{input.r2} \
        --input_type {params.type} \
        -o {output.profile} \
        --mapout {output.mapout} \
        --nproc {params.threads}
        --ignore_eukaryotes \
        --ignore_archaea \
        --db_dir {params.db} \
        -v \
        1>> {log.stdout} 2>> {log.stderr}
        """

# MAG-Based Taxonomic Classification Steps

rule GTDB-Tk:
    """
    Assign taxonomy to bins using GTDB-Tk
    """
    input:
        dastool_bins="results/{genera}/1_metagenome_assembly/6_binning/DASTool/{assembler}_individual_assembly/refined_bins/{sample}/_DASTool_bins"
    output:
        "results/{genera}/2_Taxonomic_Assignment/1_Taxonomic_Classification/GTDB-Tk/{sample}/gtdbtk_summary.tsv"
    params:
        threads = 4,
        ext = ".fa",
        prefix = "GTB-TK_",
        outdir = "results/{genera}/2_Taxonomic_Assignment/1_Taxonomic_Classification/GTDB-Tk/{sample}/"
    log:
        stdout = "logs/{genera}/2_Taxonomic_Assignment/1_Taxonomic_Classification/GTDB-Tk/{sample}/GTDB-Tk_Tax.out",
        stderr = "logs/{genera}/2_Taxonomic_Assignment/1_Taxonomic_Classification/GTDB-Tk/{sample}/GTDB-Tk_Tax.err"
    shell:
        """
        module unload miniconda 
        source activate /home/flg9/.conda/envs/gtdbtk-2.4.1
        export GTDBTK_DB=/vast/palmer/pi/turner/data/db/gtdbtk-2.4.1

        gtdbtk classify_wf \
        --genome_dir {input.dastool_bins} \
        --out_dir {params.outdir} \
        --cpus {params.threads} \
        --prefix {params.prefix} \
        -x {params.ext}
        --debug \
        1> {log.stdout} 2> {log.stderr}
        """



# misc code

rule bracken_abundance:
    """
    Estimate relative abundance using mmseq taxonomic classification output from unbinned contigs
    """
    input:
        ""
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