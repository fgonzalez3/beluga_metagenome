import os
import glob
import pandas as pd

configfile: "config/assembly.yaml"

samples_df = pd.read_csv("tsv/test_beluga_raw_reads.tsv", sep="\t")
SAMPLES = samples_df["sample_id"].tolist()
READS = {row.sample_id: {"r1": row.r1, "r2": row.r2} for row in samples_df.itertuples()}
ASSEMBLERS=config["assembler"]

        ############################# Metagenome Assembly Outputs #############################

def assembly_outputs():
    outputs = []
    for assembler in config["assembler"]:
        for sample in SAMPLES:
            if assembler == "spades":
                outputs.append(f"results/{config["genera"]}/1_metagenome_assembly/2_assembly/spades/individual_metagenome_assembly/{sample}/contigs.fasta")
            elif assembler == "megahit":
                outputs.append(f"results/{config["genera"]}/1_metagenome_assembly/2_assembly/megahit/individual_metagenome_assembly/{sample}/final.contigs.fa")
    return outputs

def contig_deduplication_outputs():
    outputs=[]
    for assembler in config["assembler"]:
        for sample in SAMPLES:
            if assembler == "spades":
                outputs.append(f"results/{config["genera"]}/1_metagenome_assembly/3_dedup_contigs/spades/individual_metagenome_assembly/{sample}/{sample}_DEDUP95.fasta")
            elif assembler == "megahit":
                outputs.append(f"results/{config["genera"]}/1_metagenome_assembly/3_dedup_contigs/megahit/individual_metagenome_assembly/{sample}/{sample}_DEDUP95.fasta")
    return outputs

def align_reads_to_individual_assemblies():
    outputs=[]
    for assembler in config["assembler"]:
        for sample in SAMPLES:
            if assembler == "spades":
                outputs.append(f"results/{config["genera"]}/1_metagenome_assembly/4_align_reads_to_contigs/spades/contig_index_individual_assemblies/{sample}/{sample}_indexed_contig.1.bt2"),
                outputs.append(f"results/{config["genera"]}/1_metagenome_assembly/4_align_reads_to_contigs/spades/contig_index_individual_assemblies/{sample}/{sample}_indexed_contig.2.bt2"),
                outputs.append(f"results/{config["genera"]}/1_metagenome_assembly/4_align_reads_to_contigs/spades/contig_index_individual_assemblies/{sample}/{sample}_indexed_contig.3.bt2"),
                outputs.append(f"results/{config["genera"]}/1_metagenome_assembly/4_align_reads_to_contigs/spades/contig_index_individual_assemblies/{sample}/{sample}_indexed_contig.4.bt2"),
                outputs.append(f"results/{config["genera"]}/1_metagenome_assembly/4_align_reads_to_contigs/spades/contig_index_individual_assemblies/{sample}/{sample}_indexed_contig.rev.1.bt2"),
                outputs.append(f"results/{config["genera"]}/1_metagenome_assembly/4_align_reads_to_contigs/spades/contig_index_individual_assemblies/{sample}/{sample}_indexed_contig.rev.2.bt2"),
                outputs.append(f"results/{config["genera"]}/1_metagenome_assembly/4_align_reads_to_contigs/spades/contig_read_alignment_individual_assemblies/{sample}_aligned_sorted.bam"),
                outputs.append(f"results/{config["genera"]}/1_metagenome_assembly/4_align_reads_to_contigs/spades/contig_read_alignment_individual_assemblies/{sample}_aligned_sorted.bam.bai")

            elif assembler == "megahit":
                outputs.append(f"results/{config["genera"]}/1_metagenome_assembly/4_align_reads_to_contigs/megahit/contig_index_individual_assemblies/{sample}/{sample}_indexed_contig.1.bt2"),
                outputs.append(f"results/{config["genera"]}/1_metagenome_assembly/4_align_reads_to_contigs/megahit/contig_index_individual_assemblies/{sample}/{sample}_indexed_contig.2.bt2"),
                outputs.append(f"results/{config["genera"]}/1_metagenome_assembly/4_align_reads_to_contigs/megahit/contig_index_individual_assemblies/{sample}/{sample}_indexed_contig.3.bt2"),
                outputs.append(f"results/{config["genera"]}/1_metagenome_assembly/4_align_reads_to_contigs/megahit/contig_index_individual_assemblies/{sample}/{sample}_indexed_contig.4.bt2"),
                outputs.append(f"results/{config["genera"]}/1_metagenome_assembly/4_align_reads_to_contigs/megahit/contig_index_individual_assemblies/{sample}/{sample}_indexed_contig.rev.1.bt2"),
                outputs.append(f"results/{config["genera"]}/1_metagenome_assembly/4_align_reads_to_contigs/megahit/contig_index_individual_assemblies/{sample}/{sample}_indexed_contig.rev.2.bt2"),
                outputs.append(f"results/{config["genera"]}/1_metagenome_assembly/4_align_reads_to_contigs/megahit/contig_read_alignment_individual_assemblies/{sample}_aligned_sorted.bam"),
                outputs.append(f"results/{config["genera"]}/1_metagenome_assembly/4_align_reads_to_contigs/megahit/contig_read_alignment_individual_assemblies/{sample}_aligned_sorted.bam.bai")
    return outputs

def evaluate_individual_assemblies():
    outputs=[]
    for assembler in config["assembler"]:
        for sample in SAMPLES:
            if assembler == "spades":
                outputs.append(f"results/{config["genera"]}/1_metagenome_assembly/5_evaluate_assemblies/spades/filter_individual_assemblies/{sample}/assembly_DEDUP95_m1500.fasta"),
                outputs.append(f"results/{config["genera"]}/1_metagenome_assembly/5_evaluate_assemblies/QUAST/report.html")
            elif assembler == "megahit":
                outputs.append(f"results/{config["genera"]}/1_metagenome_assembly/5_evaluate_assemblies/megahit/filter_individual_assemblies/{sample}/assembly_DEDUP95_m1500.fasta"),
                outputs.append(f"results/{config["genera"]}/1_metagenome_assembly/5_evaluate_assemblies/QUAST/report.html")
    return outputs

def binning_individual_assemblies():
    outputs=[]
    for assembler in config["assembler"]:
        for sample in SAMPLES:
            if assembler == "spades":

                # Concoct Outputs
                outputs.append(f"results/{config["genera"]}/1_metagenome_assembly/6_binning/concoct/spades_individual_assembly/{sample}/contigs_20k.fa"),
                outputs.append(f"results/{config["genera"]}/1_metagenome_assembly/6_binning/concoct/spades_individual_assembly/{sample}/contigs_20k.bed"),
                outputs.append(f"results/{config["genera"]}/1_metagenome_assembly/6_binning/concoct/spades_individual_assembly/{sample}/coverage_table.tsv"),
                outputs.append(f"results/{config["genera"]}/1_metagenome_assembly/6_binning/concoct/spades_individual_assembly/{sample}/bins/clustering.csv"),
                outputs.append(f"results/{config["genera"]}/1_metagenome_assembly/6_binning/concoct/spades_individual_assembly/{sample}/bins/clustering_merged.csv"),
                outputs.append(f"results/{config["genera"]}/1_metagenome_assembly/6_binning/concoct/spades_individual_assembly/{sample}/fasta_bins/manifest.txt"),

                # Metabat Outputs
                outputs.append(f"results/{config["genera"]}/1_metagenome_assembly/6_binning/metabat/spades_individual_assembly/{sample}/METABAT.txt"),
                outputs.append(f"results/{config["genera"]}/1_metagenome_assembly/6_binning/metabat/spades_individual_assembly/{sample}/bins/check.txt"),

                # MaxBin Outputs
                outputs.append(f"results/{config["genera"]}/1_metagenome_assembly/6_binning/maxbin/spades_individual_assembly/{sample}/maxbin.txt"),
                outputs.append(f"results/{config["genera"]}/1_metagenome_assembly/6_binning/maxbin/spades_individual_assembly/{sample}/bins/check.txt"),

                # SemiBin Outputs
                outputs.append(f"results/{config["genera"]}/1_metagenome_assembly/6_binning/semibin2/spades_individual_assembly/single_binning/{sample}/check.tsv"),

                # DASTool Outputs
                outputs.append(f"results/{config["genera"]}/1_metagenome_assembly/6_binning/DASTool/spades_individual_assembly/contigs2bin/{sample}.metabat.contigs2bin.tsv"),
                outputs.append(f"results/{config["genera"]}/1_metagenome_assembly/6_binning/DASTool/spades_individual_assembly/contigs2bin/{sample}.concoct.contigs2bin.tsv"),
                outputs.append(f"results/{config["genera"]}/1_metagenome_assembly/6_binning/DASTool/spades_individual_assembly/contigs2bin/{sample}.maxbin.contigs2bin.tsv"),
                outputs.append(f"results/{config["genera"]}/1_metagenome_assembly/6_binning/DASTool/spades_individual_assembly/contigs2bin/{sample}.metabat.contigs2bin.tsv"),
                outputs.append(f"results/{config["genera"]}/1_metagenome_assembly/6_binning/DASTool/spades_individual_assembly/refined_bins/{sample}/check.tsv"),

                # MetaWRAP Outputs
                #outputs.append(f"results/{config["genera"]}/1_metagenome_assembly/6_binning/MetaWRAP/spades_individual_assembly/{sample}/bin_refinement/Binning_refiner.stats"),

                # CheckM Outputs
                outputs.append(f"results/{config["genera"]}/1_metagenome_assembly/6_binning/CheckM/spades_individual_assembly/{sample}/lineage.ms"),
                outputs.append(f"results/{config["genera"]}/1_metagenome_assembly/6_binning/CheckM/spades_individual_assembly/{sample}/lineage_results.tsv"),
                outputs.append(f"results/{config["genera"]}/1_metagenome_assembly/6_binning/CheckM/spades_individual_assembly/{sample}/qa_results.tsv"),
                outputs.append(f"results/{config["genera"]}/1_metagenome_assembly/6_binning/CheckM2/spades_individual_assembly/{sample}/quality_report.tsv"),
                outputs.append(f"results/{config["genera"]}/1_metagenome_assembly/6_binning/GUNC/spades_individual_assembly/{sample}/check.txt")

            elif assembler == "megahit":
                # Concoct Outputs
                outputs.append(f"results/{config["genera"]}/1_metagenome_assembly/6_binning/concoct/megahit_individual_assembly/{sample}/contigs_20k.fa"),
                outputs.append(f"results/{config["genera"]}/1_metagenome_assembly/6_binning/concoct/megahit_individual_assembly/{sample}/contigs_20k.bed"),
                outputs.append(f"results/{config["genera"]}/1_metagenome_assembly/6_binning/concoct/megahit_individual_assembly/{sample}/coverage_table.tsv"),
                outputs.append(f"results/{config["genera"]}/1_metagenome_assembly/6_binning/concoct/megahit_individual_assembly/{sample}/bins/clustering.csv"),
                outputs.append(f"results/{config["genera"]}/1_metagenome_assembly/6_binning/concoct/megahit_individual_assembly/{sample}/bins/clustering_merged.csv"),
                outputs.append(f"results/{config["genera"]}/1_metagenome_assembly/6_binning/concoct/megahit_individual_assembly/{sample}/fasta_bins/manifest.txt"),

                # Metabat Outputs
                outputs.append(f"results/{config["genera"]}/1_metagenome_assembly/6_binning/metabat/megahit_individual_assembly/{sample}/METABAT.txt"),
                outputs.append(f"results/{config["genera"]}/1_metagenome_assembly/6_binning/metabat/megahit_individual_assembly/{sample}/bins/check.txt"),

                # MaxBin Outputs
                outputs.append(f"results/{config["genera"]}/1_metagenome_assembly/6_binning/maxbin/megahit_individual_assembly/{sample}/maxbin.txt"),
                outputs.append(f"results/{config["genera"]}/1_metagenome_assembly/6_binning/maxbin/megahit_individual_assembly/{sample}/bins/check.txt"),

                # SemiBin Outputs
                outputs.append(f"results/{config["genera"]}/1_metagenome_assembly/6_binning/semibin2/megahit_individual_assembly/single_binning/{sample}/check.tsv"),

                # DASTool Outputs
                outputs.append(f"results/{config["genera"]}/1_metagenome_assembly/6_binning/DASTool/megahit_individual_assembly/contigs2bin/{sample}.metabat.contigs2bin.tsv"),
                outputs.append(f"results/{config["genera"]}/1_metagenome_assembly/6_binning/DASTool/megahit_individual_assembly/contigs2bin/{sample}.concoct.contigs2bin.tsv"),
                outputs.append(f"results/{config["genera"]}/1_metagenome_assembly/6_binning/DASTool/megahit_individual_assembly/contigs2bin/{sample}.maxbin.contigs2bin.tsv"),
                outputs.append(f"results/{config["genera"]}/1_metagenome_assembly/6_binning/DASTool/megahit_individual_assembly/contigs2bin/{sample}.metabat.contigs2bin.tsv"),
                outputs.append(f"results/{config["genera"]}/1_metagenome_assembly/6_binning/DASTool/megahit_individual_assembly/refined_bins/{sample}/check.tsv"),

                # MetaWRAP Outputs
                #outputs.append(f"results/{config["genera"]}/1_metagenome_assembly/6_binning/MetaWRAP/megahit_individual_assembly/{sample}/bin_refinement/Binning_refiner.stats"),

                # CheckM Outputs
                outputs.append(f"results/{config["genera"]}/1_metagenome_assembly/6_binning/CheckM/megahit_individual_assembly/{sample}/lineage.ms"),
                outputs.append(f"results/{config["genera"]}/1_metagenome_assembly/6_binning/CheckM/megahit_individual_assembly/{sample}/lineage_results.tsv"),
                outputs.append(f"results/{config["genera"]}/1_metagenome_assembly/6_binning/CheckM/megahit_individual_assembly/{sample}/qa_results.tsv"),
                outputs.append(f"results/{config["genera"]}/1_metagenome_assembly/6_binning/CheckM2/megahit_individual_assembly/{sample}/quality_report.tsv"),
                outputs.append(f"results/{config["genera"]}/1_metagenome_assembly/6_binning/GUNC/megahit_individual_assembly/{sample}/check.txt")
    return outputs

rule all:
    input:
        # 1. Read pre-processing pipeline 
        expand("results/{genera}/1_metagenome_assembly/1_pre_processing/adapter_trimming/{sample}/{sample}_val_1.fq.gz", sample=SAMPLES, genera=config["genera"]),
        expand("results/{genera}/1_metagenome_assembly/1_pre_processing/adapter_trimming/{sample}/{sample}_val_2.fq.gz", sample=SAMPLES, genera=config["genera"]),
        expand("results/{genera}/1_metagenome_assembly/1_pre_processing/adapter_trimming/{sample}/{sample}_val_1_fastqc.html", sample=SAMPLES, genera=config["genera"]),
        expand("results/{genera}/1_metagenome_assembly/1_pre_processing/adapter_trimming/{sample}/{sample}_val_2_fastqc.html", sample=SAMPLES, genera=config["genera"]),
        expand("results/{genera}/1_metagenome_assembly/1_pre_processing/aggregate_qc_data/multiqc_report.html", genera=config["genera"]),
        expand("results/{genera}/1_metagenome_assembly/1_pre_processing/mask_beluga_host_genome/beluga_genome_sequence_masked.fa", genera=config["genera"]),
        expand("results/{genera}/1_metagenome_assembly/1_pre_processing/align_to_beluga_host_genome/{sample}/mapped.fq", sample=SAMPLES, genera=config["genera"]),
        expand("results/{genera}/1_metagenome_assembly/1_pre_processing/align_to_beluga_host_genome/{sample}/unmapped.fq", sample=SAMPLES, genera=config["genera"]),
        expand("results/{genera}/1_metagenome_assembly/1_pre_processing/align_to_beluga_host_genome/{sample}/unmapped_R1.fq", sample=SAMPLES, genera=config["genera"]),
        expand("results/{genera}/1_metagenome_assembly/1_pre_processing/align_to_beluga_host_genome/{sample}/unmapped_R2.fq", sample=SAMPLES, genera=config["genera"]),
        expand("results/{genera}/1_metagenome_assembly/1_pre_processing/mask_human_host_genome/human_genome_sequence_masked.fa", genera=config["genera"]),
        expand("results/{genera}/1_metagenome_assembly/1_pre_processing/align_to_human_host_genome/{sample}/mapped.fq", sample=SAMPLES, genera=config["genera"]),
        expand("results/{genera}/1_metagenome_assembly/1_pre_processing/align_to_human_host_genome/{sample}/unmapped.fq", sample=SAMPLES, genera=config["genera"]),
        expand("results/{genera}/1_metagenome_assembly/1_pre_processing/align_to_human_host_genome/{sample}/unmapped_R1.fq", sample=SAMPLES, genera=config["genera"]),
        expand("results/{genera}/1_metagenome_assembly/1_pre_processing/align_to_human_host_genome/{sample}/unmapped_R2.fq", sample=SAMPLES, genera=config["genera"]),
        expand("results/{genera}/1_metagenome_assembly/1_pre_processing/dedup_reads/{sample}/{sample}_host_removed_dedup_R1.fastq", sample=SAMPLES, genera=config["genera"]),
        expand("results/{genera}/1_metagenome_assembly/1_pre_processing/dedup_reads/{sample}/{sample}_host_removed_dedup_R2.fastq", sample=SAMPLES, genera=config["genera"]),

        # 2. Metagenome assembly pipeline
        assembly_outputs(),

        # 3. Deduplicate assembled contigs
        contig_deduplication_outputs(),

        # 4. Align PE reads to assembled contigs
        align_reads_to_individual_assemblies(),

        # 5. Evaluate assemblies
        evaluate_individual_assemblies(),

        # 6. Binning
        binning_individual_assemblies(),

        # 1. Taxonomic Classification
        expand("results/{genera}/2_Taxonomic_Assignment/1_Taxonomic_Classification/Kraken2/{sample}/{sample}.kraken2",sample=SAMPLES, genera=config["genera"]),
        expand("results/{genera}/2_Taxonomic_Assignment/1_Taxonomic_Classification/Kraken2/{sample}/{sample}.kreport2",sample=SAMPLES, genera=config["genera"]),
        expand("results/{genera}/2_Taxonomic_Assignment/1_Taxonomic_Classification/Bracken/{sample}/bracken_level_{level}.txt",sample=SAMPLES, genera=config["genera"], level=["P", "C", "O", "F", "G", "S"]),
        expand("results/{genera}/2_Taxonomic_Assignment/1_Taxonomic_Classification/Bracken/{sample}/bracken_merged.csv", sample=SAMPLES, genera=config["genera"]),
        expand("results/{genera}/2_Taxonomic_Assignment/1_Taxonomic_Classification/Kaiju/{sample}/kaiju.out", sample=SAMPLES, genera=config["genera"]),
        expand("results/{genera}/2_Taxonomic_Assignment/1_Taxonomic_Classification/Kaiju/{sample}/kaiju_summary_{level}.tsv", 
        sample=SAMPLES, 
        genera=config["genera"],
        level = ["phylum","class","order","family","genus","species"])

include: "pipelines/1_Metagenome_Assembly/1_pre_processing.smk"
include: "pipelines/1_Metagenome_Assembly/2_assembly.smk"
include: "pipelines/1_Metagenome_Assembly/3_dedup_contigs.smk"
include: "pipelines/1_Metagenome_Assembly/4_align_reads_to_contigs.smk"
include: "pipelines/1_Metagenome_Assembly/5_evaluate_assemblies.smk"
include: "pipelines/1_Metagenome_Assembly/6_binning.smk"
include: "pipelines/2_Taxonomic_Assignment/1_taxonomic_classification.smk"
#include: pipelines/virome_characterization.smk