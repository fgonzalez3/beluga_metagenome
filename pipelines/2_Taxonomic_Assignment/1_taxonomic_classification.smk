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

checkpoint Kraken2:
    """
    Classify taxonomy for reads with Kraken2
    """
    input:
        r1 = "results/{genera}/1_metagenome_assembly/1_pre_processing/dedup_reads/{sample}/{sample}_host_removed_dedup_R1.fastq",
        r2 = "results/{genera}/1_metagenome_assembly/1_pre_processing/dedup_reads/{sample}/{sample}_host_removed_dedup_R2.fastq"
    output:
        stdout = "results/{genera}/2_Taxonomic_Assignment/1_Taxonomic_Classification/Kraken2/{sample}/k2_output.txt",
        report = "results/{genera}/2_Taxonomic_Assignment/1_Taxonomic_Classification/Kraken2/{sample}/k2_report.txt"
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