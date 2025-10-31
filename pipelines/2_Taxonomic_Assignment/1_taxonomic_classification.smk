import os
import glob
import pandas as pd

configfile: "config/assembly.yaml"

samples_df = pd.read_csv("tsv/test_beluga_raw_reads.tsv", sep="\t")
SAMPLES = samples_df["sample_id"].tolist()
READS = {row.sample_id: {"r1": row.r1, "r2": row.r2} for row in samples_df.itertuples()}

# This workflow will largely be split into two parts - :
    # 1. Read-based taxonomic classification
    # 2. Taxonomic classification of individual MAGs

# Read-Based Taxonomic Classification Steps

rule Kraken2:
    """
    Classify taxonomy for reads with Kraken2
    """
    input:
        r1 = "results/{genera}/1_metagenome_assembly/1_pre_processing/dedup_reads/{sample}/{sample}_host_removed_dedup_R1.fastq",
        r2 = "results/{genera}/1_metagenome_assembly/1_pre_processing/dedup_reads/{sample}/{sample}_host_removed_dedup_R2.fastq"
    output:
        stdout = "results/{genera}/2_Taxonomic_Assignment/1_Taxonomic_Classification/Kraken2/{sample}/{sample}.kraken2",
        report = "results/{genera}/2_Taxonomic_Assignment/1_Taxonomic_Classification/Kraken2/{sample}/{sample}.kreport2"
    params:
        threads = 4,
        db = "/vast/palmer/pi/turner/data/db/kraken2_GTDBv220"
    log:
        stdout = "logs/{genera}/2_Taxonomic_Assignment/1_Taxonomic_Classification/Kraken2/{sample}/Kraken2_Tax.out",
        stderr = "logs/{genera}/2_Taxonomic_Assignment/1_Taxonomic_Classification/Kraken2/{sample}/Kraken2_Tax.err"
    shell:
        """
        /vast/palmer/pi/turner/flg9/conda_envs/kraken2/kraken2 \
        --paired \
        --db {params.db} \
        --output {output.stdout} \
        --report {output.report} \
        {input.r1} {input.r2} \
        1>> {log.stdout} 2>> {log.stderr}
        """

rule bracken_build:
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
        executables = "/vast/palmer/pi/turner/flg9/conda_envs/kraken2/kraken2"
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

rule bracken_abund:
    """
    Use Kraken2's report file to estimate microbial abundance
    """
    input:
        "results/{genera}/2_Taxonomic_Assignment/1_Taxonomic_Classification/Kraken2/{sample}/{sample}.kreport2"
    output:
        out1 = expand("results/{{genera}}/2_Taxonomic_Assignment/1_Taxonomic_Classification/Bracken/{{sample}}/{{sample}}_bracken_level_{level}.txt", level=["P", "C", "O", "F", "G", "S"]),
        out2 = "results/{genera}/2_Taxonomic_Assignment/1_Taxonomic_Classification/Bracken/{sample}/{sample}_bracken_combined.txt",
        out3 = "results/{genera}/2_Taxonomic_Assignment/1_Taxonomic_Classification/Bracken/all_combined_bracken.txt"
    params:
        readlen = 150,
        threshold = 10,
        lvls = "P,C,O,F,G,S",
        genera=config["genera"],
        db = "/vast/palmer/pi/turner/data/db/kraken2_GTDBv220"

    log:
        stdout = "logs/{genera}/2_Taxonomic_Assignment/1_Taxonomic_Classification/Bracken/{sample}/Bracken_abund.out",
        stderr = "logs/{genera}/2_Taxonomic_Assignment/1_Taxonomic_Classification/Bracken/{sample}/Bracken_abund.err"
    shell:
        """
        module unload miniconda
        source activate /home/flg9/.conda/envs/bracken

        # Define taxonomic levels that we want abundance counts for 
        levels = ("P" "C" "O" "F" "G" "S")

        # Loop through each of the taxon levels and output a separate file for each 
        for level in "${{levels[@]}}"; do
            bracken \
            -d {params.db} \
            -i {input} \
            -o {output.out1} \
            -l $level \
            -r {params.read_len} \
            1>> {log.stdout} 2>> {log.stderr}
        done

        # Combine all levels for the current sample
        cat \
        {results/{wildcards.genera}/2_Taxonomic_Assignment/1_Taxonomic_Classification/Bracken/{wildcards.sample}/{wildcards.sample}_bracken_level_*.txt} > {output.out2}

        # Merge all of the files into one
        find results/{wildcards.genera}/2_Taxonomic_Assignment/1_Taxonomic_Classification/Bracken/ -name "*_bracken_level_*.txt" -exec cat {{}} + > {output.out3}
        """