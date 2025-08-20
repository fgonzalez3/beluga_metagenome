import os
import glob
import pandas as pd

rule all:
    input:
        "results/{genera}/2_binning/concoct/{sample}/concoct_output/clustering_merged.csv",
        "results/{genera}/2_binning/metabat/{sample}/CONCOCT.*.fa",
        "results/{genera}/2_binning/metabat/{sample}/concoct_clustering_gt1000.csv",
        "results/{genera}/2_binning/metabat/METABAT.txt",
        "results/{genera}/2_binning/metabat/{sample}/METABAT.*.fa",
        "results/{genera}/2_binning/maxbin/maxbin.txt",
        "results/{genera}/2_binning/metabat/{sample}/MAXBIN.*.fa",
        "results/{genera}/2_binning/maxbin/{sample}/MAXBIN.summary",
        "results/{genera}/2_binning/aggregate_bins/reporter_files/metabat_associations.tsv",
        "results/{genera}/2_binning/aggregate_bins/reporter_files/maxbin_associations.tsv",
        "results/{genera}/2_binning/aggregate_bins/reporter_files/concoct_associations.tsv",
        "results/{genera}/2_binning/aggregate_bins/reporter_files/DASTOOL.*.fa",
        "results/{genera}/magpurify/{sample}/clean_bins.fa",
        "results/{genera}/binning_qc/{sample}/quality_report.tsv",
        ""

rule concoct:
    """
    Group assembled contigs into bins that represent individual genomes or closely related organisms using Concoct

    *I am supplying the original unfiltered contigs as input for each of my binning tools here, because I am selecting 
    for contigs at minimum of 1.5kb anyways. This should match the eval length in the last rule of assembly.smk.*
    """
    input:
        og_contigs = "results/{genera}/1_assembly/dedup_contigs/{sample}/{sample}_DEDUP95.fasta",
        bams = "results/{genera}/1_assembly/contig_read_alignment/{sample}_aligned_sorted.bam"
    output:
        bins = "results/{genera}/2_binning/concoct/{sample}/concoct_output/clustering_merged.csv",
        "results/{genera}/2_binning/metabat/{sample}/CONCOCT.*.fa"
    params:
        genera=config["genera"],
        fasta_bins_dir = "results/{genera}/2_binning/concoct/{sample}/concoct_output/fasta_bins",
        threads = 8,
        outdir = "results/{genera}/2_binning/concoct/{sample}/concoct_output/CONCOCT",
        contig_len = 20000,
        min_len = 1500
    log:
        stdout = "logs/{genera}/2_binning/concoct/{sample}/concoct.out",
        stderr = "logs/{genera}/2_binning/concoct/{sample}/concoct.err"
    shell:
        """
        module unload minicondasq
        source activate /home/flg9/.conda/envs/concoct_env

        # 1. Shred contigs greater than 20kb into smaller subsets
        cut_up_fasta.py \
        {input.og_contigs} -c {params.contig_len} -o 0 --merge_last -b contigs_20k.bed > contigs_20k.fa \
        1>> {log.stdout} 2>> {log.stderr}

        # 2. Pull alignment data
        concoct_coverage_table.py \
        contigs_20k.bed {input.bams} > coverage_table.tsv \
        1>> {log.stdout} 2>> {log.stderr}
 
        # 3. Bin
        concoct \
        --composition_file contigs_20k.fa --coverage_file coverage_table.tsv \
        -b {params.outdir} -t {params.threads} -l {params.min_len} -d \
        1>> {log.stdout} 2>> {log.stderr}

        # 4. Merge subcontig clustering into original contig clustering
        merge_cutup_clustering.py \
        {params.outdir}/clustering_gt20000.csv > {params.outdir}/clustering_merged.csv \
        1>> {log.stdout} 2>> {log.stderr}

        # 5. Extract bins as individual FASTAs
        mkdir -p \
        {params.fasta_bins_dir} 

        extract_Fasta_bins.py \
        {input.og_contigs} {params.outdir}/clustering_merged.csv --output_path {params.fasta_bins_dir} \
        1>> {log.stdout} 2>> {log.stderr}
        """

rule metabat2:
    """
    Group assembled contigs into bins that represent individual genomes or closely related organisms usin Metabat
    """
    input:
        contigs = "results/{genera}/1_assembly/dedup_contigs/{sample}/{sample}_DEDUP95.fasta",
        bams = "results/{genera}/1_assembly/contig_read_alignment/{sample}_aligned_sorted.bam"
    output:
        depth_file = "results/{genera}/2_binning/metabat/METABAT.txt",
        bins = "results/{genera}/2_binning/metabat/{sample}/METABAT.*.fa"
    params:
        genera=config["genera"],
        threads = 8,
        min_size = 1500,
        outdir = "results/{genera}/2_binning/metabat/{sample}/METABAT"
    log:
        stdout = "logs/{genera}/2_binning/metabat/{sample}/metabat.out",
        stderr = "logs/{genera}/2_binning/metabat/{sample}/metabat.err"
    shell:
        """
        module unload miniconda
        module load docker/6.0.1

        # 1. Set the work directory to ensure docker starts where it should
        WORKDIR = $(pwd)

        # 2. Generation of depth files for binning using previously sorted bams
        docker run --rm \
        --workdir $WORKDIR \
        --volume $WORKDIR:$WORKDIR \
        jgi_summarize_bam_contig_depths \
        --outputDepth {params.depth_file} {input.bams} \
        1>> {log.stdout} 2>> {log.stderr}

        # 3. Bin using previously generated contigs and depth file
        docker run --rm \
        --workdir $WORKDIR \
        --volume $WORKDIR:$WORKDIR 
        metabat/metabat:latest runMetaBat.sh \
        -t {params.threads} -m {params.min_size} \
        -i {input.contigs} \
        -a {output.depth_file} \
        -o {params.outdir} \
        --verbose --debug \
        1>> {log.stdout} 2>> {log.stderr}
        """

rule maxbin2:
    """
    Group assembled contigs into bins that represent individual genomes or closely related organisms usin MaxBin2
    """
    input:
        contigs = "results/{genera}/1_assembly/dedup_contigs/{sample}/{sample}_DEDUP95.fasta",
        r1 = "results/{genera}/1_assembly/dedup_reads/{sample}/{sample}_host_removed_dedup_R1.fastq",
        r2 = "results/{genera}/1_assembly/dedup_reads/{sample}/{sample}_host_removed_dedup_R2.fastq",
        metabat_depth_file = "results/{genera}/2_binning/metabat/metabat.txt"
    output:
        maxbin_depth_file = "results/{genera}/2_binning/maxbin/maxbin.txt",
        "results/{genera}/2_binning/metabat/{sample}/MAXBIN.*.fa",
        "results/{genera}/2_binning/maxbin/{sample}/MAXBIN.summary"
    params:
        threads=4,
        genera = config["genera"],
        contig_len = 1500,
        outdir = "results/{params.genera}/2_binning/maxbin/{sample}/MAXBIN"
    log:
        stdout = "logs/{genera}/2_binning/maxbin/{sample}/maxbin.out",
        stderr = "logs/{genera}/2_binning/maxbin/{sample}/maxbin.err"
    shell:
        """
        module unload miniconda
        module load MaxBin/2.2.7-gompi-2020b 

        # 1. Reformat Metabat depth file as to use as input for Maxbin
        cut \
        -f1,4,6,8,10 {input.metabat_epth_file} > {output.maxbin_depth_file}

        # 2. Bin
        run_MaxBin.pl \
        -thread {params.threads} --min_contig_length {params.contig_len} \
        -contig {input.contigs} -reads {input.r1} -reads2 {input.r2} 
        --abund {out.maxbin_depth_file} \
        -out {params.outdir} \
        1> {log.stdout} 2> {log.stderr}
        """

rule DASTool:
    """
    Calculate an optimized set of bins from multiple binning tools that were previously used
    """
    input:
        maxbin_contigs = "results/{genera}/2_binning/metabat/{sample}/MAXBIN.*.fa",
        concoct_csv = "results/{genera}/2_binning/concoct/{sample}/concoct_output/clustering_merged.csv",
        metabat_contigs = "results/{genera}/2_binning/metabat/{sample}/METABAT.*.fa",
        contigs = "results/{genera}/1_assembly/dedup_contigs/{sample}/{sample}_DEDUP95.fasta"
    output:
        metabat_summ = "results/{genera}/2_binning/aggregate_bins/reporter_files/metabat_associations.tsv",
        maxbin_summ = "results/{genera}/2_binning/aggregate_bins/reporter_files/maxbin_associations.tsv",
        concoct_summ = "results/{genera}/2_binning/aggregate_bins/reporter_files/concoct_associations.tsv"
        "results/{genera}/2_binning/aggregate_bins/reporter_files/DASTOOL.*.fa"
    params:
        genera=config["genera"],
        threads = 8,
        engine = "diamond", 
        dastool_out = "results/{genera}/2_binning/aggregate_bins/{sample}/aggregate_bin_{sample}"
    log:
        stdout = "logs/{genera}/2_binning/aggregate_bins/{sample}/das_tool.out",
        stderr = "logs/{genera}/2_binning/aggregate_bins/{sample}/das_tool.err"
    shell:
        """
        module unload miniconda 
        source activate XXXXXXXXX

        # 1. Generate reporter file of bins and bin quality as input for DasTool

        # First, Maxbin
        Fasta_to_Contigs2Bin.sh \
        -i {input.maxbin_contigs} -e fasta > {output.maxbin_summ}

        # Second, Concoct
        perl \
        -pe "s/,/\tconcoct./g;" {input.concoct_csv} > {output.concoct_summ}

        # Last, Metabat
        Fasta_to_Contigs2Bin.sh \
        -i {input.metabat_contigs} -e fasta > {output.metabat_summ}

        # 2. Bin de-replication with Das_Tool

        mkdir -p \
        {params.dastool_out}

        DAS_Tool \
        -i {output.metabat_summ},{output.maxbin_summ},{output.concoct_summ} \ 
        -l Metabat,Maxbin,Concoct -o {params.dastool_out} -c {input.contigs} \
        --write_bin_evals --write_bins --threads {params.threads} --debug \
        --write_unbinned --write_bins --search_engine {params.engine} \
        1> {log.stdout} 2> {log.stderr}
        """

rule magpurify2:
    """
    Refine bins further by removing contaminants with MagPurify2
    """
    input:
        maxbin_bins=lambda wildcards: sorted(glob.glob(f"results/{wildcards.genera}/binning/{wildcards.sample}/MAXBIN.*.fasta")),
        metabat_bins=lambda wildcards: sorted(glob.glob(f"results/{genera}/2_binning/metabat/{sample}/METABAT.*.fa")),
        concoct_bins=lambda wildcards: sorted(glob.glob(f"results/{genera}/2_binning/metabat/{sample}/CONCOCT.*.fa")),
        dastool_bins=lambda wildcards: sorted(glob.glob(f"results/{genera}/2_binning/aggregate_bins/reporter_files/DASTOOL.*.fa")),
        magpurify_bins=lambda wildcards: sorted(glob.glob(f"results/{genera}/magurify/{sample}/clean_bins.fa")),
        metawrap_bins= "",
        bams = "results/{genera}/1_assembly/contig_read_alignment/{sample}_aligned_sorted.bam"

    output:
        "results/{genera}/magpurify/{sample}/clean_bins.fa"
    params:
        genera=config["genera"],
        threads=4,
        outdir = "results/{genera}/magpurify/{sample}"
    log:
        stdout = "logs/{genera}/magpurify/{sample}/magpurify.out",
        stderr = "logs/{genera}/magpurify/{sample}/magpurify.err"
    shell:
        """
        module unload miniconda 
        source activate XXXXXX

        # 1. Run the composition module
        # This module identifies putative contaminants using tetranucleotide frequencies
        
        magpurify2 composition \
        {input.maxbin_bins} {input.concoct_bins} {input.metabat_bins} \
        {input.dastool_bins} {input.magpurify_bins} {input.metawrap_bins} \
        -t {params.threads} {params.outdir} \
        1>> {log.stdout} 2>> {log.stderr}

        # 2. Run the coverage module
        # This module assesses contig relatedness via differential contig coverage

        magpurify2 coverage \
        {input.maxbin_bins} {input.concoct_bins} {input.metabat_bins} \
        {input.dastool_bins} {input.magpurify_bins} {input.metawrap_bins} \
        --bam_files {input.bams} {params.outdir} \
        1>> {log.stdout} 2>> {log.stderr}

        # 3. Run the taxonomy module
        # This module assess contig contamination via taxonomical analysis

        magpurify2 taxonomy \
        {input.maxbin_bins} {input.concoct_bins} {input.metabat_bins} \
        {input.dastool_bins} {input.magpurify_bins} {input.metawrap_bins} \
        {params.outdir} magpurify2DB \
        1>> {log.stdout} 2>> {log.stderr}

        # 4. Run the filter module
        # This module filters out contaminants

        magpurify2 filter \
        {input.maxbin_bins} {input.concoct_bins} {input.metabat_bins} \
        {input.dastool_bins} {input.magpurify_bins} {input.metawrap_bins} \
        {params.outdir} filtered_genomes \
        1>> {log.stdout} 2>> {log.stderr}
        """

rule refine_bins_MetaWRAP:
    """
    Bin with binning module from MetaWrap to compare with the above methods
    """
    input:
        assembly = "results/{genera}/1_assembly/assembly_eval/{sample}/metaspades_assembly_DEDUP95_m1500.fasta",
        fwd = lambda wildcards: READS[wildcards.sample]["r1"],
        rev = lambda wildcards: READS[wildcards.sample]["r2"]
    output:
        bin1 = "results/{genera}/metawrap/binning/metabat2_bins",
        bin2 = "results/{genera}/metawrap/binning/maxbin2_bins",
        bin3 = "results/{genera}/metawrap/binning/concoct_bins"
    params:
        threads = 10,
        binning_outdir = "results/{genera}/metawrap/binning",
        refinement_outdir = "results/{genera}/metawrap/bin_refinement",
        blobology_outdir = "results/{genera}/metawrap/blobology,
        abund_outdir = "results/{genera}/metawrap/abund",
        reassembly_outdir = "results/{genera}/metawrap/reassembly",
        classify_outdir = "results/{genera}/metawrap/classify",
        annot_outdir = "results/{genera}/metawrap/annotate"
    log:
        stdout = "logs/{genera}/metawrap/{sample}/metawrap.out",
        stderr = "logs/{genera}/metawrap/{sample}/metawrap.err"
    shell:
        """
        module unload miniconda 
        source activate XXXXXX

        # 1. Run binning module
        metaWRAP binning \
        -a {input.assembly} -t {params.threads} -o {params.binning_outdir} \
        --metabat2 --maxbin2 --concoct {input.fwd} {input.rev} --run-checkm \
        1>> {log.stdout} 2>> {log.stderr}

        # 2. Run bin refinement module
        metaWRAP bin_refinement \
        -o {params.refinement_outdir} -t {params.threads} \
        -A {output.bin1} -B {output.bin2} -C {output.bin3} \
        1>> {log.stdout} 2>> {log.stderr}

        # 3. Visualize communities in our newly refined bins
        metaWRAP blobology \
        -a {input.assembly} -t {params.threads} -o {params.blobology_outdir} \ 
        {input.fwd} {input.rev} \
        1>> {log.stdout} 2>> {log.stderr}

        # 4. Plot abundance of bins
        metaWRAP quant_bins \
        -b {params.refinement_outdir}/metawrap_bins -o {params.abund_outdir} \
        -a {input.assembly} {input.fwd} {input.rev} \
        1>> {log.stdout} 2>> {log.stderr}

        # 5. Re-assemble MAGs using newly populated bins
        metaWRAP reassemble_bins \
        -o {params.reassembly_outdir} -1 {input.fwd} -2 {input.rev} \
        -t {params.threads} -b {params.refinement_outdir}/metawrap_bins \
        1>> {log.stdout} 2>> {log.stderr}

        # 6. Determine the taxonomy of each refined bin
        metaWRAP classify_bins \
        -b {params.reassembly_outdir}/reassembled_bins \
        -o {params.classify_outdir} -t {params.threads} \
        1>> {log.stdout} 2>> {log.stderr}

        # 7. Functionally annotate refined bins
        metaWRAP annotate_bins \
        -o {params.annot_outdir} -t {params.threads} \
        -b {params.reassembly_outdir}/reassembled_bins \
        1>> {log.stdout} 2>> {log.stderr}
        """ 

rule bin_quality_check:
    """
    Check bin quality with CheckM2 from our custom binning methods & compare with those of MetaWrap
    """
    input:
        # Compare between individual and refined bin datasets
        maxbin_bins=lambda wildcards: sorted(glob.glob(f"results/{wildcards.genera}/binning/{wildcards.sample}/MAXBIN.*.fasta")),
        metabat_bins=lambda wildcards: sorted(glob.glob(f"results/{genera}/2_binning/metabat/{sample}/METABAT.*.fa")),
        concoct_bins=lambda wildcards: sorted(glob.glob(f"results/{genera}/2_binning/metabat/{sample}/CONCOCT.*.fa")),
        dastool_bins=lambda wildcards: sorted(glob.glob(f"results/{genera}/2_binning/aggregate_bins/reporter_files/DASTOOL.*.fa")),
        magpurify_bins=lambda wildcards: sorted(glob.glob(f"results/{genera}/magurify/{sample}/clean_bins.fa")),
        metawrap_bins= "results/{genera}/metawrap/reassembly/reassembled_bins"
    output:
        "results/{genera}/binning_qc/{sample}/quality_report.tsv"
    params:
        threads=4,
        outdir="results/{genera}/binning_qc/{sample}",
    log:
        stdout = "logs/{genera}/binning_qc/{sample}/checkm2.out",
        stderr = "logs/{genera}/binning_qc/{sample}/checkm2.err"
    shell:
        """
        module unload miniconda
        source activate /home/flg9/.conda/envs/checkm2
        export CHECKM2DB="/vast/palmer/pi/turner/data/db/CheckM2/CheckM2_database"

        mkdir -p \
        {params.output_dir}

        # Assess bins using CheckM2 ML models
        checkm2 predict \
        --threads {params.threads} \
        --input {input.maxbin_bins} {input.concoct_bins} {input.metabat_bins} \
        {input.dastool_bins} {input.magpurify_bins} {input.metawrap_bins} \
        --output_directory {params.outdir} \
        1>> {log.stdout} 2>> {log.stderr}
        """

rule align_reads_to_bins:
    """
    Filter out any reads that align to our pre-assembled and refined bins
    """
    input:
    output:
    params:
    log:
    shell:
        """
        """

rule assemble_bin_reads:
    """
    Re-assemble reads that aligned to our pre-assembled and refined bins
    """
    input:
    output:
    params:
    log:
    shell:
        """
        """

rule bin_reassemblies:
    """
    Bin and refine re-assemblies
    """
    input:
    output:
    params:
    log:
    shell:
        """
        """

rule GUNC:
    """
    Visualize binning accuracy
    """
    input:
    output:
    params:
    log:
    shell:
        """
        """