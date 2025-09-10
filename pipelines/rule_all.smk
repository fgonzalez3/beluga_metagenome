import os
import glob
import pandas as pd

configfile: "config/assembly.yaml"

samples_df = pd.read_csv("tsv/beluga_raw_reads.tsv", sep="\t")
SAMPLES = samples_df["sample_id"].tolist()
READS = {row.sample_id: {"r1": row.r1, "r2": row.r2} for row in samples_df.itertuples()}

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
        expand("results/{genera}/1_pre_processing/split_unmapped_reads_beluga/{sample}/unmapped_R1.fq", sample=SAMPLES, genera=config["genera"]),
        expand("results/{genera}/1_pre_processing/split_unmapped_reads_beluga/{sample}/unmapped_R2.fq", sample=SAMPLES, genera=config["genera"]),
        expand("results/{genera}/1_pre_processing/mask_human_host_genome/human_genome_sequence_masked.fa", genera=config["genera"]),
        expand("results/{genera}/1_pre_processing/align_to_human_host_genome/{sample}/mapped.fq", sample=SAMPLES, genera=config["genera"]),
        expand("results/{genera}/1_pre_processing/align_to_human_host_genome/{sample}/unmapped.fq", sample=SAMPLES, genera=config["genera"]),
        expand("results/{genera}/1_pre_processing/split_unmapped_reads_human/{sample}/unmapped_R1.fq", sample=SAMPLES, genera=config["genera"]),
        expand("results/{genera}/1_pre_processing/split_unmapped_reads_human/{sample}/unmapped_R2.fq", sample=SAMPLES, genera=config["genera"]),
        expand("results/{genera}/1_pre_processing/dedup_reads/{sample}/{sample}_host_removed_dedup_R1.fastq", sample=SAMPLES, genera=config["genera"]),
        expand("results/{genera}/1_pre_processing/dedup_reads/{sample}/{sample}_host_removed_dedup_R2.fastq", sample=SAMPLES, genera=config["genera"]),

        # 2. Metagenome assembly pipeline
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
        expand("results/{genera}/5_evaluate_assemblies/filter_individual_assemblies/{sample}/metaspades_assembly_DEDUP95_m1500.fasta", sample=SAMPLES, genera=config["genera"]),
        expand("results/{genera}/5_evaluate_assemblies/filter_individual_assemblies/{sample}/megahit_assembly_DEDUP95_m1500.fasta", sample=SAMPLES, genera=config["genera"]),
        expand("results/{genera}/5_evaluate_assemblies/individual_assembly_eval/assembly_stats.csv", genera=config["genera"]),
        expand("results/{genera}/5_evaluate_assemblies/individual_assembly_eval/report.html", genera=config["genera"]),

        # 6. Binning
        expand("results/{genera}/3_binning/concoct/SPAdes_individual_assembly/{sample}/CONCOCT.*.fa", sample=SAMPLES, genera=config["genera"]),
        expand("results/{genera}/3_binning/concoct/SPAdes_individual_assembly/{sample}/concoct_output/clustering_merged.csv", sample=SAMPLES, genera=config["genera"]),
        expand("results/{genera}/3_binning/concoct/megahit_individual_assembly/{sample}/CONCOCT.*.fa", sample=SAMPLES, genera=config["genera"]),
        expand("results/{genera}/3_binning/concoct/megahit_individual_assembly/{sample}/concoct_output/clustering_merged.csv", sample=SAMPLES, genera=config["genera"]),
        expand("results/{genera}/3_binning/metabat/SPAdes_individual_assembly/{sample}/METABAT.txt", sample=SAMPLES, genera=config["genera"]),
        expand("results/{genera}/3_binning/metabat/SPAdes_individual_assembly/{sample}/METABAT.*.fa", sample=SAMPLES, genera=config["genera"]),
        expand("results/{genera}/3_binning/metabat/megahit_individual_assembly/{sample}/METABAT.txt", sample=SAMPLES, genera=config["genera"]),
        expand("results/{genera}/3_binning/metabat/megahit_individual_assembly/{sample}/METABAT.*.fa", sample=SAMPLES, genera=config["genera"]),
        expand("results/{genera}/3_binning/maxbin/SPAdes_individual_assembly/maxbin.txt", genera=config["genera"]),
        expand("results/{genera}/3_binning/maxbin/megahit_individual_assembly/maxbin.txt", genera=config["genera"]),
        expand("results/{genera}/3_binning/maxbin/SPAdes_individual_assembly/{sample}/MAXBIN.*.fa", sample=SAMPLES, genera=config["genera"]),
        expand("results/{genera}/3_binning/maxbin/SPAdes_individual_assembly/{sample}/MAXBIN.summary", sample=SAMPLES, genera=config["genera"]),
        expand("results/{genera}/3_binning/maxbin/megahit_individual_assembly/{sample}/MAXBIN.*.fa", sample=SAMPLES, genera=config["genera"]),
        expand("results/{genera}/3_binning/maxbin/megahit_individual_assembly/{sample}/MAXBIN.summary", sample=SAMPLES, genera=config["genera"]),
        expand("results/{genera}/3_binning/semibin2/SPAdes_individual_assembly/generate_concatenated_db/concatenated.fa", genera=config["genera"]),
        expand("results/{genera}/3_binning/semibin2/megahit_individual_assembly/generate_concatenated_db/concatenated.fa", genera=config["genera"])
        expand("results/{genera}/3_binning/semibin2/SPAdes_individual_assembly/align_to_concatenated_db/{sample}/{sample}_indexed_contig.1.bt2", sample=SAMPLES, genera=config["genera"]),
        expand("results/{genera}/3_binning/semibin2/SPAdes_individual_assembly/align_to_concatenated_db/{sample}/{sample}_indexed_contig.2.bt2", sample=SAMPLES, genera=config["genera"]),
        expand("results/{genera}/3_binning/semibin2/SPAdes_individual_assembly/align_to_concatenated_db/{sample}/{sample}_indexed_contig.3.bt2", sample=SAMPLES, genera=config["genera"]),
        expand("results/{genera}/3_binning/semibin2/SPAdes_individual_assembly/align_to_concatenated_db/{sample}/{sample}_indexed_contig.4.bt2", sample=SAMPLES, genera=config["genera"]),
        expand("results/{genera}/3_binning/semibin2/SPAdes_individual_assembly/align_to_concatenated_db/{sample}/{sample}_indexed_contig.rev.1.bt2", sample=SAMPLES, genera=config["genera"]),
        expand("results/{genera}/3_binning/semibin2/SPAdes_individual_assembly/align_to_concatenated_db/{sample}/{sample}_indexed_contig.rev.2.bt2", sample=SAMPLES, genera=config["genera"]),
        expand("results/{genera}/3_binning/semibin2/SPAdes_individual_assembly/align_to_concatenated_db/{sample}_aligned_sorted.bam", sample=SAMPLES, genera=config["genera"]),
        expand("results/{genera}/3_binning/semibin2/megahit_individual_assembly/align_to_concatenated_db/{sample}/{sample}_indexed_contig.1.bt2", sample=SAMPLES, genera=config["genera"]),
        expand("results/{genera}/3_binning/semibin2/megahit_individual_assembly/align_to_concatenated_db/{sample}/{sample}_indexed_contig.2.bt2", sample=SAMPLES, genera=config["genera"]),
        expand("results/{genera}/3_binning/semibin2/megahit_individual_assembly/align_to_concatenated_db/{sample}/{sample}_indexed_contig.3.bt2", sample=SAMPLES, genera=config["genera"]),
        expand("results/{genera}/3_binning/semibin2/megahit_individual_assembly/align_to_concatenated_db/{sample}/{sample}_indexed_contig.4.bt2", sample=SAMPLES, genera=config["genera"]),
        expand("results/{genera}/3_binning/semibin2/megahit_individual_assembly/align_to_concatenated_db/{sample}/{sample}_indexed_contig.rev.1.bt2", sample=SAMPLES, genera=config["genera"]),
        expand("results/{genera}/3_binning/semibin2/megahit_individual_assembly/align_to_concatenated_db/{sample}/{sample}_indexed_contig.rev.2.bt2", sample=SAMPLES, genera=config["genera"]),
        expand("results/{genera}/3_binning/semibin2/megahit_individual_assembly/align_to_concatenated_db/{sample}_aligned_sorted.bam", sample=SAMPLES, genera=config["genera"]),
        expand("results/{genera}/3_binning/semibin2/SPAdes_individual_assembly/features_and_model/{sample}/data_split.csv", sample=SAMPLES, genera=config["genera"]),
        expand("results/{genera}/3_binning/semibin2/SPAdes_individual_assembly/features_and_model/{sample}/data.csv", sample=SAMPLES, genera=config["genera"]),
        expand("results/{genera}/3_binning/semibin2/megahit_individual_assembly/features_and_model/{sample}/data_split.csv", sample=SAMPLES, genera=config["genera"]),
        expand("results/{genera}/3_binning/semibin2/megahit_individual_assembly/features_and_model/{sample}/data.csv", sample=SAMPLES, genera=config["genera"]),
        expand("results/{genera}/3_binning/semibin2/SPAdes_individual_assembly/train_model/{sample}/model.pt", sample=SAMPLES, genera=config["genera"]),
        expand("results/{genera}/3_binning/semibin2/megahit_individual_assembly/train_model/{sample}/model.pt", sample=SAMPLES, genera=config["genera"]),
        expand("results/{genera}/3_binning/SPAdes_individual_assembly/semibin2/bin/{sample}/bin.*.fa", sample=SAMPLES, genera=config["genera"]),
        expand("results/{genera}/3_binning/megahit_individual_assembly/semibin2/bin/{sample}/bin.*.fa", sample=SAMPLES, genera=config["genera"]),
        expand("results/{genera}/3_binning/SPAdes_individual_assembly/aggregate_bins/{sample}/metabat_associations.tsv", sample=SAMPLES, genera=config["genera"]),
        expand("results/{genera}/3_binning/SPAdes_individual_assembly/aggregate_bins/{sample}/maxbin_associations.tsv", sample=SAMPLES, genera=config["genera"]),
        expand("results/{genera}/3_binning/SPAdes_individual_assembly/aggregate_bins/{sample}/concoct_associations.tsv", sample=SAMPLES, genera=config["genera"]),
        expand("results/{genera}/3_binning/SPAdes_individual_assembly/aggregate_bins/{sample}/semibin_associations.tsv", sample=SAMPLES, genera=config["genera"]),
        expand("results/{genera}/3_binning/SPAdes_individual_assembly/aggregate_bins/{sample}/DASTOOL.*.fa", sample=SAMPLES, genera=config["genera"]),
        expand("results/{genera}/3_binning/megahit_individual_assembly/aggregate_bins/{sample}/metabat_associations.tsv", sample=SAMPLES, genera=config["genera"]),
        expand("results/{genera}/3_binning/megahit_individual_assembly/aggregate_bins/{sample}/maxbin_associations.tsv", sample=SAMPLES, genera=config["genera"]),
        expand("results/{genera}/3_binning/megahit_individual_assembly/aggregate_bins/{sample}/concoct_associations.tsv", sample=SAMPLES, genera=config["genera"]),
        expand("results/{genera}/3_binning/megahit_individual_assembly/aggregate_bins/{sample}/semibin_associations.tsv", sample=SAMPLES, genera=config["genera"]),
        expand("results/{genera}/3_binning/megahit_individual_assembly/aggregate_bins/{sample}/DASTOOL.*.fa", sample=SAMPLES, genera=config["genera"]),
        expand("results/{genera}/3_binning/SPAdes_individual_assembly/binning_qc/{sample}/quality_report.tsv", sample=SAMPLES, genera=config["genera"]),
        expand("results/{genera}/3_binning/megahit_individual_assembly/binning_qc/{sample}/quality_report.tsv", sample=SAMPLES, genera=config["genera"])


# Pipelines to call on 
include: "pipelines/1_Metagenome_Assembly_And_Evaluation/1_pre_processing.smk"
include: "pipelines/1_Metagenome_Assembly_And_Evaluation/2_assembly.smk"
include: "pipelines/1_Metagenome_Assembly_And_Evaluation/3_dedup_contigs.smk"
include: "pipelines/1_Metagenome_Assembly_And_Evaluation/4_align_reads_to_contigs.smk"
include: "pipelines/1_Metagenome_Assembly_And_Evaluation/5_binning.smk"
#include: pipelines/2_Taxonomic_Assignment_And_AMR_Surveillance/1_taxonomic_classification.smk
#include: pipelines/virome_characterization.smk