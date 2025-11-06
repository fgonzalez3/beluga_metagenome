# This workflow will largely be split into two parts - :
    # 1. Read-based taxonomic classification
    # 2. Taxonomic classification of individual MAGs

# Read-Based Taxonomic Classification Steps

rule Kraken2:
    """
    Classify taxonomy for reads with Kraken2
    The GTDB-Tk database for this is ~500gb, so >600gb of RAM are needed to run this
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

rule bracken_merge:
    """
    Merge Bracken abundance files into one master file per sample
    """
    input:
        expand("results/{genera}/2_Taxonomic_Assignment/1_Taxonomic_Classification/Bracken/{sample}/bracken_level_{level}.txt",
        sample=SAMPLES,
        genera=config["genera"],
        level=["P", "C", "O", "F", "G", "S"])
    output:
        "results/{genera}/2_Taxonomic_Assignment/1_Taxonomic_Classification/Bracken/{sample}/bracken_merged.csv"
    shell:
        """
        cat {input} > {output}
        """

rule Kaiju_Taxonomy:
    """
    Classify taxonomy for reads with Kaiju
    The nr database used here is ~180gb, so at least 200gb RAM are needed to run this
    """
    input:
        r1 = "results/{genera}/1_metagenome_assembly/1_pre_processing/dedup_reads/{sample}/{sample}_host_removed_dedup_R1.fastq",
        r2 = "results/{genera}/1_metagenome_assembly/1_pre_processing/dedup_reads/{sample}/{sample}_host_removed_dedup_R2.fastq"
    output:
        "results/{genera}/2_Taxonomic_Assignment/1_Taxonomic_Classification/Kaiju/{sample}/kaiju.out"
    params:
        mode = "nr",
        node = "nodes.dmp",
        refseq_index = "/vast/palmer/pi/turner/data/db/kaiju/nr/kaiju_db_nr.fmi",
        max_exact_matches = 12, # conservative params that result in closer precision to Kraken
        min_score = 70, # conservative params that result in closer precision to Kraken
        mismatches = 5 # run this on greedy-5 mode for highest sensitivity at tradeoff of slightly lower precision
    log:
        stdout = "logs/{genera}/2_Taxonomic_Assignment/1_Taxonomic_Classification/Kaiju/{sample}/Kaiju_Tax.out",
        stderr = "logs/{genera}/2_Taxonomic_Assignment/1_Taxonomic_Classification/Kaiju/{sample}/Kaiju_Tax.err"
    shell:
        """
        module unload miniconda
        source activate /vast/palmer/pi/turner/flg9/conda_envs/kaiju

        # Run kaiju-makedb -s nr to download nr database for reference
        # Other interesting databases include plasmids and viruses dbs
        
        kaiju \
        -t {params.node} \
        -f {params.refseq_index} \
        -i {input.r1} \
        -j {input.r2} \
        -o {output} \
        -m {params.max_exact_matches} \
        -s {params.min_score} \
        -e {params.mismatches} \
        -v \
        1>> {log.stdout} 2>> {log.stderr}
        """