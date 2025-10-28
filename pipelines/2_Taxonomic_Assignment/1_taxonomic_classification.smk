# This workflow will largely be split into two parts - :
    # 1. Read-based taxonomic classification
    # 2. Taxonomic classification of individual MAGs
    
rule all:
    input:
        expand("results/{genera}/2_Taxonomic_Assignment/1_Taxonomic_Classification/{sample}/Kaiju/kaiju.out", 
        genera=config["genera"],
        sample=SAMPLES),
        expand("results/{genera}/2_Taxonomic_Assignment/1_Taxonomic_Classification/Kaiju/{sample}/kaiju_summary.tsv",
        genera=config["genera"],
        sample=SAMPLES)

# Read-Based Taxonomic Classification Steps

rule Kaiju_Taxonomy:
    """
    Classify taxonomy for reads with Kaiju
    """
    input:
        r1 = "results/{genera}/1_pre_processing/dedup_reads/{sample}/{sample}_host_removed_dedup_R1.fastq",
        r2 = "results/{genera}/1_pre_processing/dedup_reads/{sample}/{sample}_host_removed_dedup_R2.fastq"
    output:
        "results/{genera}/2_Taxonomic_Assignment/1_Taxonomic_Classification/Kaiju/{sample}/kaiju.out"
    params:
        mode = "nr"
        node = "nodes.dmp",
        refseq_index = "refseq/kaiju_db_nr.fmi",
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

rule Kaiju_Summary:
    """
    Summarize findings from Kaiju into a tsv report
    """
    input:
        "results/{genera}/2_Taxonomic_Assignment/1_Taxonomic_Classification/Kaiju/{sample}/kaiju.out"
    output:
        "results/{genera}/2_Taxonomic_Assignment/1_Taxonomic_Classification/Kaiju/{sample}/kaiju_summary.tsv"
    params:
        nodes = "nodes.dmp",
        names = "names.dmp",
        ranks = "superkingdom,phylum,class,order,family,genus,species"
    log:
        stdout = "logs/{genera}/2_Taxonomic_Assignment/1_Taxonomic_Classification/Kaiju/{sample}/Kaiju_Summ.out",
        stderr = "logs/{genera}/2_Taxonomic_Assignment/1_Taxonomic_Classification/Kaiju/{sample}/Kaiju_Summ.err"
    shell:
        """
        module unload miniconda
        source activate /vast/palmer/pi/turner/flg9/conda_envs/kaiju

        kaiju2table \
        -t {params.node} \
        -n {params.names} \
        -l {params.ranks} \
        -o {output} \
        1>> {log.stdout} 2>> {log.stderr}
        """