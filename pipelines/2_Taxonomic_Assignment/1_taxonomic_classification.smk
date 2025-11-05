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

rule bracken_abund:
    """
    Get abundance from Kraken2 report file
    """
    input:
        "results/{genera}/2_Taxonomic_Assignment/1_Taxonomic_Classification/Kraken2/{sample}/{sample}.kreport2"
    output:
        "results/{genera}/2_Taxonomic_Assignment/1_Taxonomic_Classification/Bracken/{sample}/bracken_level_{level}.txt"
    params:
        readlen = 150,
        db = "/vast/palmer/pi/turner/data/db/kraken2_GTDBv220"
    log:
        stdout = "logs/{genera}/2_Taxonomic_Assignment/1_Taxonomic_Classification/Kraken2/{sample}/Bracken_{level}.out",
        stderr = "logs/{genera}/2_Taxonomic_Assignment/1_Taxonomic_Classification/Kraken2/{sample}/Bracken_{level}.err"
    shell:
        """
        module unload miniconda
        source activate /home/flg9/.conda/envs/bracken

        bracken \
            -d {params.db} \
            -i {input} \
            -o {output} \
            -l {wildcards.level} \
            -r {params.readlen} \
            1>> {log.stdout} 2>> {log.stderr}
        """

