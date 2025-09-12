import os
import glob
import pandas as pd

# This pipeline contains 3 assemblies - 
    # 1. Individual assemblies where each sample undergoes its own assembly
        # individual_metagenome_assembly_spades & individual_metagenome_assembly_megahit
    # 2. Individual co-assemblies where reads are assembled by individual
        # individual_metagenome_coassembly_spades & individual_metagenome_coassembly_megahit
    # 3. Master co-assembly where all of our reads from all individuals are grouped and co-assembled
        # master_metagenome_coassembly_spades & master_metagenome_coassembly_megahit

rule individual_metagenome_assembly_spades: # done
    """
    Assemble contigs from individual samples using SPAdes
    """
    input:
        r1 = "results/{genera}/1_pre_processing/dedup_reads/{sample}/{sample}_host_removed_dedup_R1.fastq",
        r2 = "results/{genera}/1_pre_processing/dedup_reads/{sample}/{sample}_host_removed_dedup_R2.fastq"
    output:
        "results/{genera}/2_assembly/SPAdes/individual_metagenome_assembly/{sample}/contigs.fasta"
    params:
        genera=config["genera"],
        outdir = "results/{genera}/2_assembly/SPAdes/individual_metagenome_assembly/{sample}"
    log:
        stdout = "logs/{genera}/2_assembly/SPAdes/individual_metagenome_assembly/{sample}/assembly.out",
        stderr = "logs/{genera}/2_assembly/SPAdes/individual_metagenome_assembly{sample}/assembly.err"
    resources:
        mem_mb=200000,
        threads=4
    shell:
        """
        module unload miniconda
        module load SPAdes/3.15.5-GCC-12.2.0

        spades.py \
        --meta --threads {resources.threads} \
        -1 {input.r1} -2 {input.r2} -o {params.outdir} \
        1> {log.stdout} 2> {log.stderr}
        """

rule individual_metagenome_assembly_megahit: # test
    """
    Assemble contigs from individual samples with Megahit
    """
    input:
        r1 = "results/{genera}/1_pre_processing/dedup_reads/{sample}/{sample}_host_removed_dedup_R1.fastq",
        r2 = "results/{genera}/1_pre_processing/dedup_reads/{sample}/{sample}_host_removed_dedup_R2.fastq"
    output:
        "results/{genera}/2_assembly/megahit/individual_metagenome_assembly/{sample}/final.contigs.fa"
    params:
        genera=config["genera"],
        preset="meta-large",
        prefix = "final",
        outdir = "results/{genera}/2_assembly/megahit/individual_metagenome_assembly/{sample}"
    log:
        stdout = "logs/{genera}/2_assembly/megahit/individual_metagenome_assembly/{sample}/assembly.out",
        stderr = "logs/{genera}/2_assembly/megahit/individual_metagenome_assembly/{sample}/assembly.err"
    shell:
        """
        module unload miniconda
        source activate /home/flg9/.conda/envs/megahit

        megahit \
        -1 {input.r1} -2 {input.r2} -o {params.outdir} --verbose \
        --presets {params.preset} --out-prefix {params.prefix} --continue \
        1> {log.stdout} 2> {log.stderr}
        """