# This pipeline goes through the following steps - 
    # 1-3. De-duplication of contigs assembled from each of our assemblies
       # 1. Individual assemblies from each PE read file from both SPAdes and Megahit assemblies
       # 2. Co-assemblies produced from concatenating PE read files for each corresponding whale
       # 3. Co-assemblies produced from concatenating all PE read file regardless of individual

rule deduplicate_contigs_single_assemblies: # test 
    """
    Run CD-HIT to deduplicate overrepresented contigs that may have been introduced during metagenome assembly
    """
    input:
        SPAdes_single = "results/{genera}/2_assembly/SPAdes/individual_metagenome_assembly/{sample}/contigs.fasta",
        megahit_single = "results/{genera}/2_assembly/megahit/individual_metagenome_assembly/{sample}/final.contigs.fa"
    output:
        SPAdes_single_dedups = "results/{genera}/3_dedup_contigs/SPAdes_single/{sample}/{sample}_DEDUP95.fasta",
        megahit_single_dedups = "results/{genera}/3_dedup_contigs/megahit_single/{sample}/{sample}_DEDUP95.fasta"
    params:
        genera=config["genera"],
        out1 = "results/{genera}/3_dedup_contigs/SPAdes_single/{sample}",
        out2 = "results/{genera}/3_dedup_contigs/megahit_single/{sample}"
    log:
        stdout1 = "logs/{genera}/3_dedup_contigs/SPAdes_single/{sample}/dedup_contigs.out",
        stderr1 = "logs/{genera}/3_dedup_contigs/SPAdes_single/{sample}/dedup_contigs.err",
        stdout2 = "logs/{genera}/3_dedup_contigs/megahit_single/{sample}/dedup_contigs.out",
        stderr2 = "logs/{genera}/3_dedup_contigs/megahit_single/{sample}/dedup_contigs.err"
    shell:
        """
        module unload miniconda
        module load CD-HIT/4.8.1-GCC-10.2.0

        # Create output directories 
        mkdir -p {params.out1}
        mkdir -p {params.out2}

        # Deduplicate contigs from SPAdes single assemblies
        cd-hit-est \
        -i {input.SPAdes_single} -o {output.SPAdes_single_dedups} -T 4 -c 0.95 \
        1>> {log.stdout1} 2>> {log.stderr1}

        # Deduplicate contigs from Megahit single assemblies
        cd-hit-est \
        -i {input.megahit_single} -o {output.megahit_single_dedups} -T 4 -c 0.95 \
        1>> {log.stdout2} 2>> {log.stderr2}
        """

rule deduplicate_contigs_per_individidual_coassemblies: # test 
    """
    Run CD-HIT to deduplicate overrepresented contigs that may have been introduced during metagenome assembly
    """
    input:
        SPAdes_whales = "results/{genera}/2_assembly/SPAdes/individual_metagenome_coassembly/{individual}/contigs.fasta",
        megahit_whales = "results/{genera}/2_assembly/megahit/indvidual_metagenome_coassembly/{individual}/final.contigs.fa",
    output:
        SPAdes_whales_dedups = "results/{genera}/3_dedup_contigs/dedup_contigs/SPAdes_whales/{individual}/{individual}_DEDUP95.fasta",
        megahit_whales_dedups = "results/{genera}/3_dedup_contigs/dedup_contigs/megahit_whales/{individual}/{individual}_DEDUP95.fasta",
    params:
        genera=config["genera"],
        out1 = "results/{genera}/3_dedup_contigs/dedup_contigs/SPAdes_whales/{individual}",
        out2 = "results/{genera}/3_dedup_contigs/dedup_contigs/megahit_whales/{individual}"
    log:
        stdout1 = "logs/{genera}/3_dedup_contigs/dedup_contigs/SPAdes_whales/{individual}/dedup_contigs.out",
        stderr1 = "logs/{genera}/3_dedup_contigs/dedup_contigs/SPAdes_whales/{individual}/dedup_contigs.err",
        stdout2 = "logs/{genera}/3_dedup_contigs/dedup_contigs/megahit_whales/{individual}/dedup_contigs.out",
        stderr2 = "logs/{genera}/3_dedup_contigs/dedup_contigs/megahit_whales/{individual}/dedup_contigs.err"
    shell:
        """
        module unload miniconda
        module load CD-HIT/4.8.1-GCC-10.2.0

        # Create output directories 
        mkdir -p {params.out1}
        mkdir -p {params.out2}


        # Deduplicate contigs from SPAdes per individual assemblies
        cd-hit-est \
        -i {input.SPAdes_whales} -o {output.SPAdes_whales_dedups} -T 4 -c 0.95 \
        1>> {log.stdout1} 2>> {log.stderr1}

        # Megahit de-duplication 

        # Deduplicate contigs from Megahit per individual assemblies
        cd-hit-est \
        -i {input.megahit_whales} -o {output.megahit_whales_dedups} -T 4 -c 0.95 \
        1>> {log.stdout2} 2>> {log.stderr2}
        """

rule deduplicate_contigs_master_coassemblies: # test 
    """
    Run CD-HIT to deduplicate overrepresented contigs that may have been introduced during metagenome assembly
    """
    input:
        SPAdes_master = "results/{genera}/2_assembly/SPAdes/master_metagenome_coassembly/contigs.fasta",
        megahit_master = "results/{genera}/2_assembly/megahit/master_metagenome_coassembly/final.contigs.fa"
    output:
        SPAdes_master_dedups = "results/{genera}/3_dedup_contigs/dedup_contigs/SPAdes_master/DEDUP95.fasta",
        megahit_master_dedups = "results/{genera}/3_dedup_contigs/dedup_contigs/megahit_master/DEDUP95.fasta"
    params:
        genera=config["genera"],
        out1 = "results/{genera}/3_dedup_contigs/dedup_contigs/SPAdes_master",
        out2 = "results/{genera}/3_dedup_contigs/dedup_contigs/megahit_master"
    log:
        stdout1 = "logs/{genera}/3_dedup_contigs/dedup_contigs/SPAdes_master/dedup_contigs.out",
        stderr1 = "logs/{genera}/3_dedup_contigs/dedup_contigs/SPAdes_master/dedup_contigs.err",
        stdout2 = "logs/{genera}/3_dedup_contigs/dedup_contigs/megahit_master/dedup_contigs.out",
        stderr2 = "logs/{genera}/3_dedup_contigs/dedup_contigs/megahit_master/dedup_contigs.err"
    shell:
        """
        module unload miniconda
        module load CD-HIT/4.8.1-GCC-10.2.0

        # Create output directories 
        mkdir -p {params.out1}
        mkdir -p {params.out2}

        # Deduplicate contigs from SPAdes master assemblies
        cd-hit-est \
        -i {input.SPAdes_master} -o {output.SPAdes_master_dedups} -T 4 -c 0.95 \
        1>> {log.stdout1} 2>> {log.stderr1}

        # Deduplicate contigs from Megahit master assemblies
        cd-hit-est \
        -i {input.megahit_master} -o {output.megahit_master_dedups} -T 4 -c 0.95 \
        1>> {log.stdout2} 2>> {log.stderr2}