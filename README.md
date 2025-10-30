# Beluga Whale Metagenome

Created by Freddy L. Gonzalez (Turner Lab)

This repo contains data, code, results, and figures related to the analysis of captive beluga whale microbiomes at Mystic Aquarium, Mystic, CT. 

# Metadata

Metadata for all sequenced specimens is available in the 'metadata' folder.

# Pipelines

There are currently two pipelines with various sub-workflows. To run this pipeline, simply run the rule_all.smk, which will create data for the following processes:

- 1_Metagenome_Assembly
    - 1_pre_processing.smk - Includes adapter trimming and low quality read removal, filtering of host reads, and deduplication of reads
    - 2_assembly.smk - Currently includes contig assembly using both metaSPAdes and megahit
    - 3_dedup_contigs.smk - De-duplication of contigs assembled in Step 2
    - 4_align_reads_to_contigs.smk - Alignment of de-duplicated short reads to contig assemblies
    - 5_evaluate_assemblies.smk - Assessment of contig assemblies
    - 6_binning.smk - Binning, bin aggregation, and quality checks of bins
- 2_Taxonomic_Assignment
    - 1_taxonomic_classification.smk - Read- and MAG-based taxonomic classification
    - 2_phylogenomics.smk - Placement of MAGs within phylogenetic trees
    - 3_functional_profiling.smk - Gene calling and functional profiling of MAGs
    - 4_AMR_surveillance.smk - Surveying for chromosomal and mobile AMR genes
    - 5_pangenomics.smk - Depicting microbial diversity via pangenomic models

# Results

Relevant results necessary to create Figures. 

# Figures

Main and supplementary figures that will make it to the final manuscript.

# Scripts 

Contains some helper scripts for running the pipeline. 