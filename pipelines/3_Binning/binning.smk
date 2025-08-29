import os
import glob
import pandas as pd

rule concoct: #done
    """
    Group assembled contigs into bins that represent individual genomes or closely related organisms using Concoct

    *I am supplying the original unfiltered contigs as input for each of my binning tools here, because I am selecting 
    for contigs at minimum of 1.5kb anyways. This should match the eval length in the last rule of assembly.smk.*
    """
    input:
        contigs = "results/{genera}/1_assembly/dedup_contigs/{sample}/{sample}_DEDUP95.fasta",
        bams = "results/{genera}/1_assembly/contig_read_alignment/{sample}_aligned_sorted.bam"
    output:
        bins = "results/{genera}/2_binning/concoct/{sample}/CONCOCT.*.fa",
        csv = "results/{genera}/2_binning/concoct/{sample}/concoct_output/clustering_merged.csv"
    params:
        outdir = "results/{genera}/2_binning/concoct/{sample}",
        basename = "results/{genera}/2_binning/concoct/{sample}/CONCOCT",
        threads = 4,
        contig_len = 20000,
        min_len = 1500,
        img = "svg"
    log:
        stdout = "logs/{genera}/2_binning/concoct/{sample}/concoct.out",
        stderr = "logs/{genera}/2_binning/concoct/{sample}/concoct.err"
    shell:
        """
        module unload minicondasq
        source activate /home/flg9/.conda/envs/concoct_env

        # 1. Shred contigs into non-overlapping parts of equal length
        cut_up_fasta.py \
        {input.contigs} -c {params.contig_len} --merge_last -o 0 \
        -b {params.outdir}/contigs_20k.bed > {params.outdir}/contigs_20k.fa \
        1>> {log.stdout} 2>> {log.stderr}

        # 2. Generate input coverage table for CONCOCT using previously curated BED file and BAMs
        concoct_coverage_table.py \
        {params.outdir}/contigs_20k.bed {input.bams} > {params.outdir}/coverage_table.tsv \
        1>> {log.stdout} 2>> {log.stderr}
 
        # 3. Bin
        concoct \
        --composition_file {params.outdir}/contigs_20k.fa \
        --coverage_file {params.outdir}/coverage_table.tsv \
        -b {params.basename} -t {params.threads} -l {params.min_len} -d \
        1>> {log.stdout} 2>> {log.stderr}

        # 4. Merge subcontig clustering into original contig clustering
        merge_cutup_clustering.py \
        {params.outdir}/clustering_gt20000.csv > {output.csv} \
        1>> {log.stdout} 2>> {log.stderr}

        # 5. Extract bins as individual FASTAs
        extract_Fasta_bins.py \
        {input.contigs} {output.csv} \
        --output_path {params.outdir} \
        1>> {log.stdout} 2>> {log.stderr}

        # 6. Generate distance matrix between bins
        python dnadiff_dist_matrix.py \
        {params.outdir} {output.bins} \
        --plot_image_extension {params.img} \
        1>> {log.stdout} 2>> {log.stderr}
        """

rule metabat2: #done
    """
    Group assembled contigs into bins that represent individual genomes or closely related organisms using Metabat
    """
    input:
        contigs = "results/{genera}/1_assembly/dedup_contigs/{sample}/{sample}_DEDUP95.fasta",
        bams = "results/{genera}/1_assembly/contig_read_alignment/{sample}_aligned_sorted.bam"
    output:
        depth_file = "results/{genera}/2_binning/metabat/{sample}/METABAT.txt",
        bins = "results/{genera}/2_binning/metabat/{sample}/METABAT.*.fa"
    params:
        genera=config["genera"],
        threads = 4,
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
        --outputDepth {output.depth_file} {input.bams} \
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

rule maxbin2_depth: # done
    """
    Get depth file for MaxBin using previously generated read alignment to contigs
    """
    input:
        bams = expand("results/{genera}/1_assembly/contig_read_alignment/{sample}_aligned_sorted.bam", genera=config["genera"], sample=SAMPLES)
    output:
        depth_file = "results/{genera}/2_binning/maxbin/maxbin.txt"
    log:
        stdout = "logs/{genera}/2_binning/maxbin/{sample}/maxbin.out",
        stderr = "logs/{genera}/2_binning/maxbin/{sample}/maxbin.err"
    shell:
        """
        module unload miniconda
        module load MaxBin/2.2.7-gompi-2020b 

        jgi_summarize_bam_contig_depths \
        --outputDepth {output.depth_file} {input.bams} \
        1>> {log.stdout} 2>> {log.stderr}
        """

rule maxbin2: # done
    """
    Bin contigs using MaxBin and previously generated depth file
    """
    input:
        contigs = "results/{genera}/1_assembly/dedup_contigs/{sample}/{sample}_DEDUP95.fasta",
        r1 = "results/{genera}/1_assembly/dedup_reads/{sample}/{sample}_host_removed_dedup_R1.fastq",
        r2 = "results/{genera}/1_assembly/dedup_reads/{sample}/{sample}_host_removed_dedup_R2.fastq",
        maxbin_depth_file = "results/{genera}/2_binning/maxbin/maxbin.txt"
    output:
        bins = "results/{genera}/2_binning/maxbin/{sample}/MAXBIN.*.fa",
        summ = "results/{genera}/2_binning/maxbin/{sample}/MAXBIN.summary"
    params:
        threads=4,
        contig_len = 1500,
        outdir = "results/{genera}/2_binning/maxbin/{sample}/MAXBIN"
    log:
        stdout = "logs/{genera}/2_binning/maxbin/{sample}/maxbin.out",
        stderr = "logs/{genera}/2_binning/maxbin/{sample}/maxbin.err"
    shell:
        """
        module unload miniconda
        module load MaxBin/2.2.7-gompi-2020b 

        run_MaxBin.pl \
        -thread {params.threads} --min_contig_length {params.contig_len} \
        -contig {input.contigs} -reads {input.r1} -reads2 {input.r2} \
        --abund {input.maxbin_depth_file} \
        -out {params.outdir} \
        1>> {log.stdout} 2>> {log.stderr}
        """

rule semibin2_features_and_model:
    """
    Generate sequence features and train model for SemiBin2
    """
    input:
        contigs = "results/{genera}/1_assembly/dedup_contigs/{sample}/{sample}_DEDUP95.fasta",
        bams = "results/{genera}/1_assembly/contig_read_alignment/{sample}_aligned_sorted.bam"
    output:
        split = "results/{genera}/2_binning/semibin2/data_split.csv",
        csv = "results/{genera}/2_binning/semibin2/data.csv",
        model = "results/{genera}/2_binning/semibin2/model.pt"
    params:
        outdir = "results/{genera}/2_binning/semibin2",
        threads = 4
    log:
        stdout = "logs/{genera}/2_binning/semibin2/features_model.out",
        stderr = "logs/{genera}/2_binning/semibin2/features_model.err"
    shell:
        """
        module unload miniconda
        source activate /vast/palmer/pi/turner/flg9/conda_envs/semibin

        # 1. Generate sequence features
        SemiBin2 generate_sequence_features_single \
        -i {input.contigs} \
        -b {input.bams} \
        -o {params.outdir} \
        -t {params.threads} \
        1>> {log.stdout} 2>> {log.stderr}

        # 2. Train model
        SemiBin2 train_self \
        --data {output.csv} \
        --data-split {output.split} \
        -o {params.outdir} \
        -t {params.threads} \
        1>> {log.stdout} 2>> {log.stderr}
        """

rule semibin2_bin:
    """
    Bin contigs using SemiBin's co-assembly binning model.
    This modelco-assembles samples as if the pool of samples were a single sample.
    Bins are then constructed from this pool of co-assembled contigs. 
    This model is most appropriate for samples that are very similar and can be expected to contain
    overlapping sets of organisms (i.e. time series from the same habitat).
    Since our samples are from three co-habiting individuals, this seemed most appropriate.
    """
    input:
        contigs = "results/{genera}/1_assembly/dedup_contigs/{sample}/{sample}_DEDUP95.fasta",
        csv = "results/{genera}/2_binning/semibin2/data.csv",
        model = "results/{genera}/2_binning/semibin2/model.pt"
    output:
        bins = "results/{genera}/2_binning/semibin2/{sample}/bin.*.fa"
    params:
        outdir = "results/{genera}/2_binning/semibin2/{sample}",
        seq_type = "short_reads",
        GTDB_path = "/vast/palmer/pi/turner/data/db/gtdbtk-2.4.1",
        minlen = "1500",
        threads = 4
    log:
        stdout = "logs/{genera}/2_binning/semibin2/{sample}/semibin2.out",
        stderr = "logs/{genera}/2_binning/semibin2/{sample}/semibin2.err"
    shell:
        """
        module unload miniconda
        source activate /vast/palmer/pi/turner/flg9/conda_envs/semibin

        # Bin
        SemiBin2 bin_short \
        -i {input.contigs} \
        --model {input.model} \
        --data {input.csv} \
        -o {params.outdir} \
        -t {params.threads} \
        --sequencing-type={params.seq_type} \
        --reference-db-data-dir={params.GTDB_path} \
        --min-len={params.minlen} \
        1>> {log.stdout} 2>> {log.stderr}
        """