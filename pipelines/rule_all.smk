import os
import glob
import pandas as pd

configfile: "config/assembly.yaml"

samples_df = pd.read_csv("tsv/beluga_raw_reads.tsv", sep="\t")
SAMPLES = samples_df["sample_id"].tolist()
READS = {row.sample_id: {"r1": row.r1, "r2": row.r2} for row in samples_df.itertuples()}

rule all:
    input:
        # 1. Read pre-processing pipeline results 
        expand("results/{genera}/1_pre_processing/adapter_trimming/{sample}/{sample}_val_1.fq.gz", sample=SAMPLES, genera=config["genera"]),
        expand("results/{genera}/1_pre_processing/adapter_trimming/{sample}/{sample}_val_2.fq.gz", sample=SAMPLES, genera=config["genera"]),
        expand("results/{genera}/1_pre_processing/adapter_trimming/{sample}/{sample}_val_1_fastqc.html", sample=SAMPLES, genera=config["genera"]),
        expand("results/{genera}/1_pre_processing/adapter_trimming/{sample}/{sample}_val_2_fastqc.html", sample=SAMPLES, genera=config["genera"]),
        expand("results/{genera}/1_pre_processing/aggregate_qc_data/multiqc_report.html", genera=config["genera"]),
        expand("results/{genera}/1_pre_processing/mask_beluga_host_genome/beluga_genome_sequence_masked.fa", genera=config["genera"]),
        expand("results/{genera}/1_pre_processing/align_to_beluga_host_genome/{sample}/mapped.fq", sample=SAMPLES, genera=config["genera"]),
        expand("results/{genera}/1_pre_processing/align_to_beluga_host_genome/{sample}/unmapped.fq", sample=SAMPLES, genera=config["genera"]),
        expand("results/{genera}/1_pre_processing/split_unmapped_reads_beluga/{sample}/unmapped_R1.fq", sample=SAMPLES, genera=config["genera"]),
        expand("results/{genera}/1_pre_processing/split_unmapped_reads_beluga/{sample}/unmapped_R2.fq", sample=SAMPLES, genera=config["genera"]),
        expand("results/{genera}/1_pre_processing/mask_human_host_genome/human_genome_sequence_masked.fa", genera=config["genera"]),
        expand("results/{genera}/1_pre_processing/align_to_human_host_genome/{sample}/mapped.fq", sample=SAMPLES, genera=config["genera"]),
        expand("results/{genera}/1_pre_processing/align_to_human_host_genome/{sample}/unmapped.fq", sample=SAMPLES, genera=config["genera"]),
        expand("results/{genera}/1_pre_processing/split_unmapped_reads_human/{sample}/unmapped_R1.fq", sample=SAMPLES, genera=config["genera"]),
        expand("results/{genera}/1_pre_processing/split_unmapped_reads_human/{sample}/unmapped_R2.fq", sample=SAMPLES, genera=config["genera"]),
        expand("results/{genera}/1_pre_processing/dedup_reads/{sample}/{sample}_host_removed_dedup_R1.fastq", sample=SAMPLES, genera=config["genera"]),
        expand("results/{genera}/1_pre_processing/dedup_reads/{sample}/{sample}_host_removed_dedup_R2.fastq", sample=SAMPLES, genera=config["genera"]),

        # 2. Metagenome assembly pipeline results 
        expand("results/{genera}/2_assembly/SPAdes/individual_metagenome_assembly/{sample}/contigs.fasta", sample=SAMPLES, genera=config["genera"]),
        expand("results/{genera}/2_assembly/megahit/individual_metagenome_assembly/{sample}/final.contigs.fa", sample=SAMPLES, genera=config["genera"]),

        # 3. Deduplicate assembled contigs
        expand("results/{genera}/3_dedup_contigs/SPAdes_single/{sample}/{sample}_DEDUP95.fasta", sample=SAMPLES, genera=config["genera"]),
        expand("results/{genera}/3_dedup_contigs/megahit_single/{sample}/{sample}_DEDUP95.fasta", sample=SAMPLES, genera=config["genera"]),

        # 4. Align PE reads to assembled contigs 
        expand("results/{genera}/4_align_reads_to_contigs/contig_index_spades_individual_assemblies/{sample}/{sample}_indexed_contig.1.bt2", sample=SAMPLES, genera=config["genera"]),
        expand("results/{genera}/4_align_reads_to_contigs/contig_index_spades_individual_assemblies/{sample}/{sample}_indexed_contig.2.bt2", sample=SAMPLES, genera=config["genera"]),
        expand("results/{genera}/4_align_reads_to_contigs/contig_index_spades_individual_assemblies/{sample}/{sample}_indexed_contig.3.bt2", sample=SAMPLES, genera=config["genera"]),
        expand("results/{genera}/4_align_reads_to_contigs/contig_index_spades_individual_assemblies/{sample}/{sample}_indexed_contig.4.bt2", sample=SAMPLES, genera=config["genera"]),
        expand("results/{genera}/4_align_reads_to_contigs/contig_index_spades_individual_assemblies/{sample}/{sample}_indexed_contig.rev.1.bt2", sample=SAMPLES, genera=config["genera"]),
        expand("results/{genera}/4_align_reads_to_contigs/contig_index_spades_individual_assemblies/{sample}/{sample}_indexed_contig.rev.2.bt2", sample=SAMPLES, genera=config["genera"]),
        expand("results/{genera}/4_align_reads_to_contigs/contig_read_alignment_individual_assemblies_spades/{sample}_aligned_sorted.bam", sample=SAMPLES, genera=config["genera"]),
        expand("results/{genera}/4_align_reads_to_contigs/contig_index_megahit_individual_assemblies/{sample}/{sample}_indexed_contig.1.bt2", sample=SAMPLES, genera=config["genera"]),
        expand("results/{genera}/4_align_reads_to_contigs/contig_index_megahit_individual_assemblies/{sample}/{sample}_indexed_contig.2.bt2", sample=SAMPLES, genera=config["genera"]),
        expand("results/{genera}/4_align_reads_to_contigs/contig_index_megahit_individual_assemblies/{sample}/{sample}_indexed_contig.3.bt2", sample=SAMPLES, genera=config["genera"]),
        expand("results/{genera}/4_align_reads_to_contigs/contig_index_megahit_individual_assemblies/{sample}/{sample}_indexed_contig.4.bt2", sample=SAMPLES, genera=config["genera"]),
        expand("results/{genera}/4_align_reads_to_contigs/contig_index_megahit_individual_assemblies/{sample}/{sample}_indexed_contig.rev.1.bt2", sample=SAMPLES, genera=config["genera"]),
        expand("results/{genera}/4_align_reads_to_contigs/contig_index_megahit_individual_assemblies/{sample}/{sample}_indexed_contig.rev.2.bt2", sample=SAMPLES, genera=config["genera"]),
        expand("results/{genera}/4_align_reads_to_contigs/contig_read_alignment_individual_assemblies_megahit/{sample}_aligned_sorted.bam", sample=SAMPLES, genera=config["genera"]),

        # 5. Evaluate assemblies
        expand("results/{genera}/2_assembly/assembly_eval/{sample}/metaspades_assembly_DEDUP95_m1500.fasta", sample=SAMPLES, genera=config["genera"]),
        expand("results/{genera}/2_assembly/assembly_eval/{sample}/assembly_stats.csv", sample=SAMPLES, genera=config["genera"]),
        expand("results/{genera}/2_assembly/assembly_eval/{sample}/report.html", sample=SAMPLES, genera=config["genera"]),

        # 3. Binning pipeline results 
        expand("results/{genera}/3_binning/concoct/{sample}/CONCOCT.*.fa", genera=config["genera"], sample=SAMPLES),
        expand("results/{genera}/3_binning/concoct/{sample}/concoct_output/clustering_merged.csv", genera=config["genera"], sample=SAMPLES),

        expand("results/{genera}/3_binning/metabat/{sample}/METABAT.txt", genera=config["genera"], sample=SAMPLES),
        expand("results/{genera}/3_binning/metabat/{sample}/METABAT.*.fa", genera=config["genera"], sample=SAMPLES),

        expand("results/{genera}/3_binning/maxbin/maxbin.txt", genera=config["genera"], sample=SAMPLES),
        expand("results/{genera}/3_binning/maxbin/{sample}/MAXBIN.*.fa", genera=config["genera"], sample=SAMPLES),
        expand("results/{genera}/3_binning/maxbin/{sample}/MAXBIN.summary", genera=config["genera"], sample=SAMPLES),

        expand("results/{genera}/3_binning/semibin2/generate_concatenated_db/concatenated.fa", genera=config["genera"]),
        expand("results/{genera}/3_binning/semibin2/align_to_concatenated_db/{sample}/{sample}_indexed_contig.1.bt2", genera=config["genera"], sample=SAMPLES),
        expand("results/{genera}/3_binning/semibin2/align_to_concatenated_db/{sample}/{sample}_indexed_contig.2.bt2", genera=config["genera"], sample=SAMPLES),
        expand("results/{genera}/2_binning/semibin2/align_to_concatenated_db/{sample}/{sample}_indexed_contig.3.bt2", genera=config["genera"], sample=SAMPLES),
        expand("results/{genera}/3_binning/semibin2/align_to_concatenated_db/{sample}/{sample}_indexed_contig.4.bt2", genera=config["genera"], sample=SAMPLES),
        expand("results/{genera}/3_binning/semibin2/align_to_concatenated_db/{sample}/{sample}_indexed_contig.rev.1.bt2", genera=config["genera"], sample=SAMPLES),
        expand("results/{genera}/3_binning/semibin2/align_to_concatenated_db/{sample}/{sample}_indexed_contig.rev.2.bt2", genera=config["genera"], sample=SAMPLES),
        expand("results/{genera}/3_binning/semibin2/align_to_concatenated_db/{sample}_aligned_sorted.bam", genera=config["genera"], sample=SAMPLES),
        expand("results/{genera}/3_binning/semibin2/features_and_model/{sample}/data_split.csv", genera=config["genera"], sample=SAMPLES),
        expand("results/{genera}/3_binning/semibin2/features_and_model/{sample}/data.csv", genera=config["genera"], sample=SAMPLES),
        expand("results/{genera}/3_binning/semibin2/train_model/{sample}/model.pt", genera=config["genera"], sample=SAMPLES),
        expand("results/{genera}/3_binning/semibin2/bin/{sample}/bin.*.fa", genera=config["genera"], sample=SAMPLES)

include: "pipelines/1_Metagenome_Assembly_And_Evaluation/1_pre_processing.smk"
include: "pipelines/1_Metagenome_Assembly_And_Evaluation/2_assembly.smk"
include: "pipelines/1_Metagenome_Assembly_And_Evaluation/3_dedup_contigs.smk"
include: "pipelines/1_Metagenome_Assembly_And_Evaluation/4_align_reads_to_contigs.smk"
include: "pipelines/1_Metagenome_Assembly_And_Evaluation/5_binning.smk"
#include: pipelines/2_Taxonomic_Assignment_And_AMR_Surveillance/1_taxonomic_classification.smk
#include: pipelines/virome_characterization.smk