# This pipeline goes through the following steps - 
    # 1-3. De-duplication of contigs assembled from each of our assemblies
       # 1. Individual assemblies from each PE read file from both SPAdes and Megahit assemblies
       # 2. Co-assemblies produced from concatenating PE read files for each corresponding whale
       # 3. Co-assemblies produced from concatenating all PE read file regardless of individual

rule deduplicate_contigs_single_assemblies_spades: # test 
    """
    Run CD-HIT to deduplicate overrepresented contigs that may have been introduced during metagenome assembly
    """
    input:
        SPAdes_single = "results/{genera}/2_assembly/SPAdes/individual_metagenome_assembly/{sample}/contigs.fasta"
    output:
        SPAdes_single_dedups = "results/{genera}/3_dedup_contigs/SPAdes/individual_metagenome_assembly/{sample}/{sample}_DEDUP95.fasta"
    params:
        genera=config["genera"],
        out = "results/{genera}/3_dedup_contigs/SPAdes/individual_metagenome_assembly/{sample}",
    log:
        stdout = "logs/{genera}/3_dedup_contigs/SPAdes/individual_metagenome_assembly/{sample}/dedup_contigs.out",
        stderr = "logs/{genera}/3_dedup_contigs/SPAdes/individual_metagenome_assembly/{sample}/dedup_contigs.err"
    shell:
        """
        module unload miniconda
        module load CD-HIT/4.8.1-GCC-10.2.0

        # Create output directories 
        mkdir -p {params.out}

        # Deduplicate contigs from SPAdes single assemblies
        cd-hit-est \
        -i {input.SPAdes_single} -o {output.SPAdes_single_dedups} -T 4 -c 0.95 \
        1>> {log.stdout} 2>> {log.stderr}
        """

rule deduplicate_contigs_single_assemblies_megahit: # test 
    """
    Run CD-HIT to deduplicate overrepresented contigs that may have been introduced during metagenome assembly
    """
    input:
        megahit_single = "results/{genera}/2_assembly/megahit/individual_metagenome_assembly/{sample}/final.contigs.fa"
    output:
        megahit_single_dedups = "results/{genera}/3_dedup_contigs/megahit/individual_metagenome_assembly/{sample}/{sample}_DEDUP95.fasta"
    params:
        genera=config["genera"],
        out = "results/{genera}/3_dedup_contigs/megahit/individual_metagenome_assembly/{sample}"
    log:
        stdout = "logs/{genera}/3_dedup_contigs/megahit/individual_metagenome_assembly/{sample}/dedup_contigs.out",
        stderr = "logs/{genera}/3_dedup_contigs/megahit/individual_metagenome_assembly/{sample}/dedup_contigs.err"
    shell:
        """
        module unload miniconda
        module load CD-HIT/4.8.1-GCC-10.2.0

        # Create output directories 
        mkdir -p {params.out}

        # Deduplicate contigs from Megahit single assemblies
        cd-hit-est \
        -i {input.megahit_single} -o {output.megahit_single_dedups} -T 4 -c 0.95 \
        1>> {log.stdout} 2>> {log.stderr}
        """