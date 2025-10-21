import os
import glob
import pandas as pd

configfile: "config/assembly.yaml"

samples_df = pd.read_csv("tsv/test_beluga_raw_reads.tsv", sep="\t")
SAMPLES = samples_df["sample_id"].tolist()
READS = {row.sample_id: {"r1": row.r1, "r2": row.r2} for row in samples_df.itertuples()}
ASSEMBLERS=config["assembler"]

def assembly_outputs():
    outputs = []
    for assembler in config["assembler"]:
        for sample in SAMPLES:
            if assembler == "spades":
                outputs.append(f"results/{config["genera"]}/testing/2_assembly/spades/individual_metagenome_assembly/{sample}/contigs.fasta")
            elif assembler == "megahit":
                outputs.append(f"results/{config["genera"]}/testing/2_assembly/megahit/individual_metagenome_assembly/{sample}/final.contigs.fa")
    return outputs

def contig_deduplication_outputs():
    outputs=[]
    for assembler in config["assembler"]:
        for sample in SAMPLES:
            if assembler == "spades":
                outputs.append(f"results/{config["genera"]}/testing/3_dedup_contigs/spades/individual_metagenome_assembly/{sample}/{sample}_DEDUP95.fasta")
            elif assembler == "megahit":
                outputs.append(f"results/{config["genera"]}/testing/3_dedup_contigs/megahit/individual_metagenome_assembly/{sample}/{sample}_DEDUP95.fasta")
    return outputs

def align_reads_to_individual_assemblies():
    outputs=[]
    for assembler in config["assembler"]:
        for sample in SAMPLES:
            if assembler == "spades":
                outputs.append(f"results/{config["genera"]}/testing/4_align_reads_to_contigs/spades/contig_index_spades_individual_assemblies/{sample}/{sample}_indexed_contig.1.bt2"),
                outputs.append(f"results/{config["genera"]}/testing/4_align_reads_to_contigs/spades/contig_index_spades_individual_assemblies/{sample}/{sample}_indexed_contig.2.bt2"),
                outputs.append(f"results/{config["genera"]}/testing/4_align_reads_to_contigs/spades/contig_index_spades_individual_assemblies/{sample}/{sample}_indexed_contig.3.bt2"),
                outputs.append(f"results/{config["genera"]}/testing/4_align_reads_to_contigs/spades/contig_index_spades_individual_assemblies/{sample}/{sample}_indexed_contig.4.bt2"),
                outputs.append(f"results/{config["genera"]}/testing/4_align_reads_to_contigs/spades/contig_index_spades_individual_assemblies/{sample}/{sample}_indexed_contig.rev.1.bt2"),
                outputs.append(f"results/{config["genera"]}/testing/4_align_reads_to_contigs/spades/contig_index_spades_individual_assemblies/{sample}/{sample}_indexed_contig.rev.2.bt2"),
                outputs.append(f"results/{config["genera"]}/testing/4_align_reads_to_contigs/spades/contig_read_alignment_individual_assemblies_spades/{sample}_aligned_sorted.bam"),
                outputs.append(f"results/{config["genera"]}/testing/4_align_reads_to_contigs/spades/contig_read_alignment_individual_assemblies_spades/{sample}_aligned_sorted.bam.bai")
            elif assembler == "megahit":
                outputs.append(f"results/{config["genera"]}/testing/4_align_reads_to_contigs/megahit/contig_index_spades_individual_assemblies/{sample}/{sample}_indexed_contig.1.bt2"),
                outputs.append(f"results/{config["genera"]}/testing/4_align_reads_to_contigs/megahit/contig_index_spades_individual_assemblies/{sample}/{sample}_indexed_contig.2.bt2"),
                outputs.append(f"results/{config["genera"]}/testing/4_align_reads_to_contigs/megahit/contig_index_spades_individual_assemblies/{sample}/{sample}_indexed_contig.3.bt2"),
                outputs.append(f"results/{config["genera"]}/testing/4_align_reads_to_contigs/megahit/contig_index_spades_individual_assemblies/{sample}/{sample}_indexed_contig.4.bt2"),
                outputs.append(f"results/{config["genera"]}/testing/4_align_reads_to_contigs/megahit/contig_index_spades_individual_assemblies/{sample}/{sample}_indexed_contig.rev.1.bt2"),
                outputs.append(f"results/{config["genera"]}/testing/4_align_reads_to_contigs/megahit/contig_index_spades_individual_assemblies/{sample}/{sample}_indexed_contig.rev.2.bt2"),
                outputs.append(f"results/{config["genera"]}/testing/4_align_reads_to_contigs/megahit/contig_read_alignment_individual_assemblies_spades/{sample}_aligned_sorted.bam"),
                outputs.append(f"results/{config["genera"]}/testing/4_align_reads_to_contigs/megahit/contig_read_alignment_individual_assemblies_spades/{sample}_aligned_sorted.bam.bai")
    return outputs

rule all:
    input:
        # 1. Read pre-processing pipeline 
        expand("results/{genera}/1_pre_processing/adapter_trimming/{sample}/{sample}_val_1.fq.gz", sample=SAMPLES, genera=config["genera"]),
        expand("results/{genera}/1_pre_processing/adapter_trimming/{sample}/{sample}_val_2.fq.gz", sample=SAMPLES, genera=config["genera"]),
        expand("results/{genera}/1_pre_processing/adapter_trimming/{sample}/{sample}_val_1_fastqc.html", sample=SAMPLES, genera=config["genera"]),
        expand("results/{genera}/1_pre_processing/adapter_trimming/{sample}/{sample}_val_2_fastqc.html", sample=SAMPLES, genera=config["genera"]),
        expand("results/{genera}/1_pre_processing/aggregate_qc_data/multiqc_report.html", genera=config["genera"]),
        expand("results/{genera}/1_pre_processing/mask_beluga_host_genome/beluga_genome_sequence_masked.fa", genera=config["genera"]),
        expand("results/{genera}/1_pre_processing/align_to_beluga_host_genome/{sample}/mapped.fq", sample=SAMPLES, genera=config["genera"]),
        expand("results/{genera}/1_pre_processing/align_to_beluga_host_genome/{sample}/unmapped.fq", sample=SAMPLES, genera=config["genera"]),
        expand("results/{genera}/1_pre_processing/align_to_beluga_host_genome/{sample}/unmapped_R1.fq", sample=SAMPLES, genera=config["genera"]),
        expand("results/{genera}/1_pre_processing/align_to_beluga_host_genome/{sample}/unmapped_R2.fq", sample=SAMPLES, genera=config["genera"]),
        expand("results/{genera}/1_pre_processing/mask_human_host_genome/human_genome_sequence_masked.fa", genera=config["genera"]),
        expand("results/{genera}/1_pre_processing/align_to_human_host_genome/{sample}/mapped.fq", sample=SAMPLES, genera=config["genera"]),
        expand("results/{genera}/1_pre_processing/align_to_human_host_genome/{sample}/unmapped.fq", sample=SAMPLES, genera=config["genera"]),
        expand("results/{genera}/1_pre_processing/align_to_human_host_genome/{sample}/unmapped_R1.fq", sample=SAMPLES, genera=config["genera"]),
        expand("results/{genera}/1_pre_processing/align_to_human_host_genome/{sample}/unmapped_R2.fq", sample=SAMPLES, genera=config["genera"]),
        expand("results/{genera}/1_pre_processing/dedup_reads/{sample}/{sample}_host_removed_dedup_R1.fastq", sample=SAMPLES, genera=config["genera"]),
        expand("results/{genera}/1_pre_processing/dedup_reads/{sample}/{sample}_host_removed_dedup_R2.fastq", sample=SAMPLES, genera=config["genera"]),

        # 2. Metagenome assembly pipeline
        assembly_outputs(),

        # 3. Deduplicate assembled contigs
        contig_deduplication_outputs(),

        # 4. Align PE reads to assembled contigs
        align_reads_to_individual_assemblies(),

        # 5. Evaluate assemblies
       # evaluate_individual_assemblies(),

        # 6. Binning

        # Concoct Outputs
        expand("results/{genera}/6_binning/concoct/SPAdes_individual_assembly/{sample}/contigs_20k.fa", sample=SAMPLES, genera=config["genera"]),
        expand("results/{genera}/6_binning/concoct/SPAdes_individual_assembly/{sample}/contigs_20k.bed", sample=SAMPLES, genera=config["genera"]),
        expand("results/{genera}/6_binning/concoct/SPAdes_individual_assembly/{sample}/coverage_table.tsv", sample=SAMPLES, genera=config["genera"]),
        expand("results/{genera}/6_binning/concoct/SPAdes_individual_assembly/{sample}/bins/clustering.csv",sample=SAMPLES,genera=config["genera"]),
        expand("results/{genera}/6_binning/concoct/SPAdes_individual_assembly/{sample}/bins/clustering_merged.csv", sample=SAMPLES, genera=config["genera"]),
        expand("results/{genera}/6_binning/concoct/SPAdes_individual_assembly/{sample}/fasta_bins/manifest.txt", sample=SAMPLES, genera=config["genera"]),

        # Metabat Outputs
        expand("results/{genera}/6_binning/metabat/SPAdes_individual_assembly/{sample}/METABAT.txt", sample=SAMPLES, genera=config["genera"]),
        expand("results/{genera}/6_binning/metabat/SPAdes_individual_assembly/{sample}/bins/check.txt", sample=SAMPLES, genera=config["genera"]),

        # Maxbin Outputs
        expand("results/{genera}/6_binning/maxbin/SPAdes_individual_assembly/{sample}/maxbin.txt", sample=SAMPLES, genera=config["genera"]),
        expand("results/{genera}/6_binning/maxbin/SPAdes_individual_assembly/{sample}/bins/check.txt", sample=SAMPLES, genera=config["genera"]),

        # Semibin Outputs
        #expand("results/{genera}/6_binning/semibin2/SPAdes_individual_assembly/generate_concatenated_db/concatenated.fa", genera=config["genera"]),
        #expand("results/{genera}/6_binning/semibin2/SPAdes_individual_assembly/align_to_concatenated_db/{sample}/{sample}_indexed_contig.1.bt2", sample=SAMPLES, genera=config["genera"]),
        #expand("results/{genera}/6_binning/semibin2/SPAdes_individual_assembly/align_to_concatenated_db/{sample}/{sample}_indexed_contig.2.bt2", sample=SAMPLES, genera=config["genera"]),
        #expand("results/{genera}/6_binning/semibin2/SPAdes_individual_assembly/align_to_concatenated_db/{sample}/{sample}_indexed_contig.3.bt2", sample=SAMPLES, genera=config["genera"]),
        #expand("results/{genera}/6_binning/semibin2/SPAdes_individual_assembly/align_to_concatenated_db/{sample}/{sample}_indexed_contig.4.bt2", sample=SAMPLES, genera=config["genera"]),
        #expand("results/{genera}/6_binning/semibin2/SPAdes_individual_assembly/align_to_concatenated_db/{sample}/{sample}_indexed_contig.rev.1.bt2", sample=SAMPLES, genera=config["genera"]),
        #expand("results/{genera}/6_binning/semibin2/SPAdes_individual_assembly/align_to_concatenated_db/{sample}/{sample}_indexed_contig.rev.2.bt2", sample=SAMPLES, genera=config["genera"]),
        #expand("results/{genera}/6_binning/semibin2/SPAdes_individual_assembly/align_to_concatenated_db/{sample}_aligned_sorted.bam", sample=SAMPLES, genera=config["genera"]),
        #expand("results/{genera}/6_binning/semibin2/SPAdes_individual_assembly/features_and_model/{sample}/data_split.csv", sample=SAMPLES, genera=config["genera"]),
        #expand("results/{genera}/6_binning/semibin2/SPAdes_individual_assembly/features_and_model/{sample}/data.csv", sample=SAMPLES, genera=config["genera"]),
        #expand("results/{genera}/6_binning/semibin2/SPAdes_individual_assembly/train_model/{sample}/model.pt", sample=SAMPLES, genera=config["genera"]),
        expand("results/{genera}/6_binning/semibin2/SPAdes_individual_assembly/binning/{sample}/check.tsv", sample=SAMPLES, genera=config["genera"]),

        # DASTool Outputs
        expand("results/{genera}/6_binning/DASTool/SPAdes_individual_assembly/contigs2bin/{sample}.metabat.contigs2bin.tsv", sample=SAMPLES, genera=config["genera"]),
        expand("results/{genera}/6_binning/DASTool/SPAdes_individual_assembly/contigs2bin/{sample}.concoct.contigs2bin.tsv", sample=SAMPLES, genera=config["genera"]),
        expand("results/{genera}/6_binning/DASTool/SPAdes_individual_assembly/contigs2bin/{sample}.maxbin.contigs2bin.tsv", sample=SAMPLES, genera=config["genera"]),
        expand("results/{genera}/6_binning/DASTool/SPAdes_individual_assembly/contigs2bin/{sample}.metabat.contigs2bin.tsv", sample=SAMPLES, genera=config["genera"]),
        expand("results/{genera}/6_binning/DASTool/SPAdes_individual_assembly/refined_bins/{sample}/check.tsv", sample=SAMPLES, genera=config["genera"]),

        # MetaWRAP Outputs
        #expand("results/{genera}/6_binning/MetaWRAP/SPAdes_individual_assembly/{sample}/bin_refinement/Binning_refiner.stats", sample=SAMPLES, genera=config["genera"]),

        # Bin Quality Check Outputs
        expand("results/{genera}/6_binning/CheckM/SPAdes_individual_assembly/{sample}/lineage.ms", sample=SAMPLES, genera=config["genera"]),
        expand("results/{genera}/6_binning/CheckM/SPAdes_individual_assembly/{sample}/lineage_results.tsv", sample=SAMPLES, genera=config["genera"]),
        expand("results/{genera}/6_binning/CheckM/SPAdes_individual_assembly/{sample}/qa_results.tsv", sample=SAMPLES, genera=config["genera"]),
        expand("results/{genera}/6_binning/CheckM2/SPAdes_individual_assembly/{sample}/quality_report.tsv", sample=SAMPLES, genera=config["genera"]),
        expand("results/{genera}/6_binning/GUNC/SPAdes_individual_assembly/{sample}/check.txt", sample=SAMPLES, genera=config["genera"])

# Pipelines to call on 
include: "pipelines/1_Metagenome_Assembly_And_Evaluation/1_pre_processing.smk"
include: "pipelines/1_Metagenome_Assembly_And_Evaluation/2_assembly.smk"
include: "pipelines/1_Metagenome_Assembly_And_Evaluation/3_dedup_contigs.smk"
include: "pipelines/1_Metagenome_Assembly_And_Evaluation/4_align_reads_to_contigs.smk"
#include: "pipelines/1_Metagenome_Assembly_And_Evaluation/5_evaluate_assemblies.smk"
#include: "pipelines/1_Metagenome_Assembly_And_Evaluation/6_binning.smk"
#include: pipelines/2_Taxonomic_Assignment_And_AMR_Surveillance/1_taxonomic_classification.smk
#include: pipelines/virome_characterization.smk