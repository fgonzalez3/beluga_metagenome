# This pipeline goes through the following steps - 
    # 1. Binning of contigs with Concoct
    # 2. Binning of contigs with Metabat2
    # 3. Collection of depth file required for MaxBin2 from BAM files produced in prior pipeline
    # 4. Binning of contigs with MaxBin2
    # 5. Model procurement and binning with Semibin2
    # 6. Bin correction and convergence with DasTool
    # 7. Bin quality check with CheckM2

rule concoct_bins_spades:
    """
    Group assembled contigs into bins that represent individual genomes or closely related organisms using Concoct (via Docker)
    """
    input:
        contigs = "results/{genera}/3_dedup_contigs/SPAdes/individual_metagenome_assembly/{sample}/{sample}_DEDUP95.fasta",
        bams = "results/{genera}/4_align_reads_to_contigs/contig_read_alignment_individual_assemblies_spades/{sample}_aligned_sorted.bam",
        bai "results/{genera}/4_align_reads_to_contigs/contig_read_alignment_individual_assemblies_spades/{sample}_aligned_sorted.bam.bai"
    output:
        bins = "results/{genera}/6_binning/concoct/SPAdes_individual_assembly/{sample}/CONCOCT.*.fa",
        csv = "results/{genera}/6_binning/concoct/SPAdes_individual_assembly/{sample}/concoct_output/clustering_merged.csv"
    params:
        outdir = "results/{genera}/6_binning/concoct/SPAdes_individual_assembly/{sample}",
        basename = "results/{genera}/6_binning/concoct/SPAdes_individual_assembly/{sample}/CONCOCT",
        threads = 4,
        contig_len = 20000,
        min_len = 1500,
        img = "svg"
    log:
        stdout = "logs/{genera}/6_binning/concoct/SPAdes_individual_assembly/{sample}/concoct.out",
        stderr = "logs/{genera}/6_binning/concoct/SPAdes_individual_assembly/{sample}/concoct.err"
    shell:
        """
        module unload miniconda

        # 1. Shred contigs into non-overlapping parts of equal length
        apptainer exec containers/concoct-1.1.0.sif \
        cut_up_fasta.py {input.contigs} -c {params.contig_len} --merge_last -o 0 \
        -b {params.outdir}/contigs_20k.bed > {params.outdir}/contigs_20k.fa \
        1>> {log.stdout} 2>> {log.stderr}

        # 2. Generate input coverage table for CONCOCT
        apptainer exec containers/concoct-1.1.0.sif \
        concoct_coverage_table.py {params.outdir}/contigs_20k.bed {input.bams} > {params.outdir}/coverage_table.tsv \
        1>> {log.stdout} 2>> {log.stderr}

        # 3. Bin
        apptainer exec containers/concoct-1.1.0.sif \
        concoct --composition_file {params.outdir}/contigs_20k.fa \
        --coverage_file {params.outdir}/coverage_table.tsv \
        -b {params.basename} -t {params.threads} -l {params.min_len} -d \
        1>> {log.stdout} 2>> {log.stderr}

        # 4. Merge subcontig clustering into original contig clustering
        apptainer exec containers/concoct-1.1.0.sif \
        merge_cutup_clustering.py {params.outdir}/clustering_gt20000.csv > {output.csv} \
        1>> {log.stdout} 2>> {log.stderr}

        # 5. Extract bins as individual FASTAs
        apptainer exec containers/concoct-1.1.0.sif \
        extract_fasta_bins.py {input.contigs} {output.csv} --output_path {params.outdir} \
        1>> {log.stdout} 2>> {log.stderr}

        # 6. Generate distance matrix between bins
        apptainer exec containers/concoct-1.1.0.sif \
        python dnadiff_dist_matrix.py {params.outdir} {output.bins} --plot_image_extension {params.img} \
        1>> {log.stdout} 2>> {log.stderr}
        """

rule concoct_bins_megahit: #test
    """
    Group assembled contigs into bins that represent individual genomes or closely related organisms using Concoct
    """
    input:
        contigs = "results/{genera}/3_dedup_contigs/megahit/individual_metagenome_assembly/{sample}/{sample}_DEDUP95.fasta",
        bams = "results/{genera}/4_align_reads_to_contigs/contig_read_alignment_individual_assemblies_megahit/{sample}_aligned_sorted.bam",
        bai "results/{genera}/4_align_reads_to_contigs/contig_read_alignment_individual_assemblies_megahit/{sample}_aligned_sorted.bam.bai"
    output:
        bins = "results/{genera}/6_binning/concoct/megahit_individual_assembly/{sample}/CONCOCT.*.fa",
        csv = "results/{genera}/6_binning/concoct/megahit_individual_assembly/{sample}/concoct_output/clustering_merged.csv"
    params:
        outdir = "results/{genera}/6_binning/concoct/megahit_individual_assembly/{sample}",
        basename = "results/{genera}/6_binning/concoct/megahit_individual_assembly/{sample}/CONCOCT",
        threads = 4,
        contig_len = 20000,
        min_len = 1500,
        img = "svg"
    log:
        stdout = "logs/{genera}/6_binning/concoct/megahit_individual_assembly/{sample}/concoct.err",
        stderr = "logs/{genera}/6_binning/concoct/megahit_individual_assembly/{sample}/concoct.err"
    shell:
        """
        module unload minicondasq

        # 1. Shred contigs into non-overlapping parts of equal length
        apptainer exec containers/concoct-1.1.0.sif \
        cut_up_fasta.py \
        {input.contigs} -c {params.contig_len} --merge_last -o 0 \
        -b {params.outdir}/contigs_20k.bed > {params.outdir}/contigs_20k.fa \
        1>> {log.stdout} 2>> {log.stderr}

        # 2. Generate input coverage table for CONCOCT using previously curated BED file and BAMs
        apptainer exec containers/concoct-1.1.0.sif \
        concoct_coverage_table.py \
        {params.outdir}/contigs_20k.bed {input.bams} > {params.outdir}/coverage_table.tsv \
        1>> {log.stdout} 2>> {log.stderr}
 
        # 3. Bin
        apptainer exec containers/concoct-1.1.0.sif \
        concoct \
        --composition_file {params.outdir}/contigs_20k.fa \
        --coverage_file {params.outdir}/coverage_table.tsv \
        -b {params.basename} -t {params.threads} -l {params.min_len} -d \
        1>> {log.stdout} 2>> {log.stderr}

        # 4. Merge subcontig clustering into original contig clustering
        apptainer exec containers/concoct-1.1.0.sif \
        merge_cutup_clustering.py \
        {params.outdir}/clustering_gt20000.csv > {output.csv} \
        1>> {log.stdout} 2>> {log.stderr}

        # 5. Extract bins as individual FASTAs
        apptainer exec containers/concoct-1.1.0.sif \
        extract_Fasta_bins.py \
        {input.contigs} {output.csv} \
        --output_path {params.outdir} \
        1>> {log.stdout} 2>> {log.stderr}

        # 6. Generate distance matrix between bins
        apptainer exec containers/concoct-1.1.0.sif \
        python dnadiff_dist_matrix.py \
        {params.outdir} {output.bins} \
        --plot_image_extension {params.img} \
        1>> {log.stdout} 2>> {log.stderr}
        """

rule metabat2_bin_spades:
    """
    Group assembled contigs into bins that represent individual genomes or closely related organisms using Metabat2 (via Docker)
    """
    input:
        contigs = "results/{genera}/3_dedup_contigs/SPAdes/individual_metagenome_assembly/{sample}/{sample}_DEDUP95.fasta",
        bams = "results/{genera}/4_align_reads_to_contigs/contig_read_alignment_individual_assemblies_spades/{sample}_aligned_sorted.bam"
    output:
        depth_file = "results/{genera}/6_binning/metabat/SPAdes_individual_assembly/{sample}/METABAT.txt",
        bins = "results/{genera}/6_binning/metabat/SPAdes_individual_assembly/{sample}/METABAT.*.fa"
    params:
        genera=config["genera"],
        threads = 4,
        min_size = 1500,
        outdir = "results/{genera}/6_binning/metabat/SPAdes_individual_assembly/{sample}/METABAT"
    log:
        stdout = "logs/{genera}/6_binning/metabat/SPAdes_individual_assembly/{sample}/metabat.out",
        stderr = "logs/{genera}/6_binning/metabat/SPAdes_individual_assembly/{sample}/metabat.err"
    shell:
        """
        module unload miniconda

        # Generate depth file
        apptainer exec containers/metabat2-2.15.sif \
        jgi_summarize_bam_contig_depths \
        --outputDepth {output.depth_file} {input.bams} \
        1>> {log.stdout} 2>> {log.stderr}

        # Run MetaBAT2
        apptainer exec containers/metabat2-2.15.sif \
        runMetaBat.sh \
        -t {params.threads} -m {params.min_size} \
        -i {input.contigs} \
        -a {output.depth_file} \
        -o {params.outdir} \
        --verbose --debug \
        1>> {log.stdout} 2>> {log.stderr}
        """

rule metabat2_bin_megahit: #test
    """
    Group assembled contigs into bins that represent individual genomes or closely related organisms using Metabat
    """
    input:
        contigs = "results/{genera}/3_dedup_contigs/megahit/individual_metagenome_assembly/{sample}/{sample}_DEDUP95.fasta",
        bams = "results/{genera}/4_align_reads_to_contigs/contig_read_alignment_individual_assemblies_megahit/{sample}_aligned_sorted.bam"
    output:
        depth_file = "results/{genera}/6_binning/metabat/megahit_individual_assembly/{sample}/METABAT.txt",
        bins = "results/{genera}/6_binning/metabat/megahit_individual_assembly/{sample}/METABAT.*.fa"
    params:
        genera=config["genera"],
        threads = 4,
        min_size = 1500,
        outdir = "results/{genera}/6_binning/metabat/megahit_individual_assembly/{sample}/METABAT"
    log:
        stdout = "logs/{genera}/6_binning/metabat/megahit_individual_assembly/{sample}/metabat.out",
        stderr = "logs/{genera}/6_binning/metabat/megahit_individual_assembly/{sample}/metabat.err"
    shell:
        """
        module unload miniconda

        # Generate depth file
        apptainer exec containers/metabat2-2.15.sif \
        jgi_summarize_bam_contig_depths \
        --outputDepth {output.depth_file} {input.bams} \
        1>> {log.stdout} 2>> {log.stderr}

        # Run MetaBAT2
        apptainer exec containers/metabat2-2.15.sif \
        runMetaBat.sh \
        -t {params.threads} -m {params.min_size} \
        -i {input.contigs} \
        -a {output.depth_file} \
        -o {params.outdir} \
        --verbose --debug \
        1>> {log.stdout} 2>> {log.stderr}
        """

rule maxbin2_depth_spades: # test
    """
    Get depth file for MaxBin using previously generated read alignment to contigs
    """
    input:
        metabat_txt = "results/{genera}/6_binning/metabat/SPAdes_individual_assembly/{sample}/METABAT.txt"
    output:
        depth_file = "results/{genera}/6_binning/maxbin/SPAdes_individual_assembly/maxbin.txt"
    log:
        stdout = "logs/{genera}/6_binning/maxbin/SPAdes_individual_assembly/maxbin.out",
        stderr = "logs/{genera}/6_binning/maxbin/SPAdes_individual_assembly/maxbin.err"
    shell:
        """
        cut \
        -f1,4,6,8,10 {input.metabat_txt} > {output.depth_file} \
        1>> {log.stdout} 2>> {log.stderr}
        """

rule maxbin2_depth_megahit: # test
    """
    Get depth file for MaxBin using previously generated read alignment to contigs
    """
    input:
        metabat_txt = "results/{genera}/6_binning/metabat/megahit_individual_assembly/{sample}/METABAT.txt"
    output:
        depth_file = "results/{genera}/6_binning/maxbin/megahit_individual_assembly/maxbin.txt"
    log:
        stdout = "logs/{genera}/6_binning/maxbin/megahit_individual_assembly/maxbin.out",
        stderr = "logs/{genera}/6_binning/maxbin/megahit_individual_assembly/maxbin.err"
    shell:
        """
        cut \
        -f1,4,6,8,10 {input.metabat_txt} > {output.depth_file} \
        1>> {log.stdout} 2>> {log.stderr}
        """

rule maxbin2_bin_spades: # test
    """
    Bin contigs using MaxBin and previously generated depth file
    """
    input:
        contigs = "results/{genera}/3_dedup_contigs/SPAdes/individual_metagenome_assembly/{sample}/{sample}_DEDUP95.fasta",
        r1 = "results/{genera}/1_pre_processing/dedup_reads/{sample}/{sample}_host_removed_dedup_R1.fastq",
        r2 = "results/{genera}/1_pre_processing/dedup_reads/{sample}/{sample}_host_removed_dedup_R2.fastq",
        maxbin_depth_file = "results/{genera}/6_binning/maxbin/SPAdes_individual_assembly/maxbin.txt"
    output:
        bins = "results/{genera}/6_binning/maxbin/SPAdes_individual_assembly/{sample}/MAXBIN.*.fa",
        summ = "results/{genera}/6_binning/maxbin/SPAdes_individual_assembly/{sample}/MAXBIN.summary"
    params:
        threads=4,
        contig_len = 1500,
        outdir = "results/{genera}/6_binning/maxbin/SPAdes_individual_assembly/{sample}/MAXBIN"
    log:
        stdout = "logs/{genera}/6_binning/maxbin/SPAdes_individual_assembly/{sample}/maxbin.out",
        stderr = "logs/{genera}/6_binning/maxbin/SPAdes_individual_assembly/{sample}/maxbin.err"
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

rule maxbin2_bin_megahit: # test
    """
    Bin contigs using MaxBin and previously generated depth file
    """
    input:
        contigs = "results/{genera}/3_dedup_contigs/megahit/individual_metagenome_assembly/{sample}/{sample}_DEDUP95.fasta",
        r1 = "results/{genera}/1_pre_processing/dedup_reads/{sample}/{sample}_host_removed_dedup_R1.fastq",
        r2 = "results/{genera}/1_pre_processing/dedup_reads/{sample}/{sample}_host_removed_dedup_R2.fastq",
        maxbin_depth_file = "results/{genera}/6_binning/maxbin/megahit_individual_assembly/maxbin.txt"
    output:
        bins = "results/{genera}/6_binning/maxbin/megahit_individual_assembly/{sample}/MAXBIN.*.fa",
        summ = "results/{genera}/6_binning/maxbin/megahit_individual_assembly/{sample}/MAXBIN.summary"
    params:
        threads=4,
        contig_len = 1500,
        outdir = "results/{genera}/6_binning/maxbin/megahit_individual_assembly/{sample}/MAXBIN"
    log:
        stdout = "logs/{genera}/6_binning/maxbin/megahit_individual_assembly/{sample}/maxbin.out",
        stderr = "logs/{genera}/6_binning/maxbin/megahit_individual_assembly/{sample}/maxbin.err"
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

##### Semibin - binning multi-sample workflow for individual binning assemblies #####

rule semibin2_generate_concatenated_db_spades: # test
    """
    Generate concatenated FASTA file necessary for SembiBin's multi-sample binning pipeline
    """
    input:
        contigs = expand("results/{genera}/3_dedup_contigs/SPAdes/individual_metagenome_assembly/{sample}/{sample}_DEDUP95.fasta", genera=config["genera"], sample=SAMPLES)
    output:
        "results/{genera}/6_binning/semibin2/SPAdes_individual_assembly/generate_concatenated_db/concatenated.fa"
    params:
        outdir = "results/{genera}/6_binning/semibin2/SPAdes_individual_assembly/generate_concatenated_db",
        threads = 4
    log:
        stdout = "logs/{genera}/6_binning/semibin2/SPAdes_individual_assembly/generate_concatenated_db/concatenate_fa.out",
        stderr = "logs/{genera}/6_binning/semibin2/SPAdes_individual_assembly/generate_concatenated_db/concatenate_fa.err"
    shell:
        """
        module unload miniconda
        source activate /vast/palmer/pi/turner/flg9/conda_envs/semibin

        SemiBin2 concatenate_fasta \
        --input-fasta {input.contigs} \
        --output {params.outdir} --compression=none
        """

rule semibin2_generate_concatenated_db_megahit: # test
    """
    Generate concatenated FASTA file necessary for SembiBin's multi-sample binning pipeline
    """
    input:
        contigs = expand("results/{genera}/3_dedup_contigs/megahit/individual_metagenome_assembly/{sample}/{sample}_DEDUP95.fasta", genera=config["genera"], sample=SAMPLES)
    output:
        "results/{genera}/6_binning/semibin2/megahit_individual_assembly/generate_concatenated_db/concatenated.fa"
    params:
        outdir = "results/{genera}/6_binning/semibin2/megahit_individual_assembly/generate_concatenated_db",
        threads = 4
    log:
        stdout = "logs/{genera}/6_binning/semibin2/megahit_individual_assembly/generate_concatenated_db/concatenate_fa.out",
        stderr = "logs/{genera}/6_binning/semibin2/megahit_individual_assembly/generate_concatenated_db/concatenate_fa.err"
    shell:
        """
        module unload miniconda
        source activate /vast/palmer/pi/turner/flg9/conda_envs/semibin

        SemiBin2 concatenate_fasta \
        --input-fasta {input.contigs} \
        --output {params.outdir} --compression=none
        """

rule sembin2_align_to_concatenated_db_spades: # test
    """
    Align reads from each sample to our concatenated FASTA db, necessary for SemiBin pipeline
    """
    input:
        contigs = "results/{genera}/3_dedup_contigs/SPAdes/individual_metagenome_assembly/{sample}/{sample}_DEDUP95.fasta",
        r1 = "results/{genera}/1_pre_processing/dedup_reads/{sample}/{sample}_host_removed_dedup_R1.fastq",
        r2 = "results/{genera}/1_pre_processing/dedup_reads/{sample}/{sample}_host_removed_dedup_R2.fastq"
    output:
        "results/{genera}/6_binning/semibin2/SPAdes_individual_assembly/align_to_concatenated_db/{sample}/{sample}_indexed_contig.1.bt2",
        "results/{genera}/6_binning/semibin2/SPAdes_individual_assembly/align_to_concatenated_db/{sample}/{sample}_indexed_contig.2.bt2",
        "results/{genera}/6_binning/semibin2/SPAdes_individual_assembly/align_to_concatenated_db/{sample}/{sample}_indexed_contig.3.bt2",
        "results/{genera}/6_binning/semibin2/SPAdes_individual_assembly/align_to_concatenated_db/{sample}/{sample}_indexed_contig.4.bt2",
        "results/{genera}/6_binning/semibin2/SPAdes_individual_assembly/align_to_concatenated_db/{sample}/{sample}_indexed_contig.rev.1.bt2",
        "results/{genera}/6_binning/semibin2/SPAdes_individual_assembly/align_to_concatenated_db/{sample}/{sample}_indexed_contig.rev.2.bt2",
        "results/{genera}/6_binning/semibin2/SPAdes_individual_assembly/align_to_concatenated_db/{sample}_aligned_sorted.bam"
    params:
        genera=config["genera"],
        outdir = "results/{genera}/6_binning/semibin2/SPAdes_individual_assembly/align_to_concatenated_db/{sample}"
    log:
        stdout = "logs/{genera}/6_binning/semibin2/SPAdes_individual_assembly/align_to_concatenated_db/{sample}_aln.out",
        stderr = "logs/{genera}/6_binning/semibin2/SPAdes_individual_assembly/align_to_concatenated_db/{sample}_aln.err"
    shell:
        """
        module unload miniconda 
        module load Bowtie2/2.5.1-GCC-12.2.0
        module load SAMtools/1.21-GCC-12.2.0

        # 1. Build db index
        bowtie2-build \
        -f {input.contigs} {params.outdir}/{wildcards.sample}_indexed_contig \
        1>> {log.stdout} 2>> {log.stderr}

        # 2. Align reads back to SemiBin db index
        bowtie2 \
        -x {params.outdir}/{wildcards.sample}_indexed_contig -1 {input.r1} -2 {input.r2} | samtools view -b -F 4 -F 2048 | samtools sort -o {output[6]} \
        1>> {log.stdout} 2>> {log.stderr}
        """

rule sembin2_align_to_concatenated_db_megahit: # test
    """
    Align reads from each sample to our concatenated FASTA db, necessary for SemiBin pipeline
    """
    input:
        contigs = "results/{genera}/3_dedup_contigs/megahit/individual_metagenome_assembly/{sample}/{sample}_DEDUP95.fasta",
        r1 = "results/{genera}/1_pre_processing/dedup_reads/{sample}/{sample}_host_removed_dedup_R1.fastq",
        r2 = "results/{genera}/1_pre_processing/dedup_reads/{sample}/{sample}_host_removed_dedup_R2.fastq"
    output:
        "results/{genera}/6_binning/semibin2/megahit_individual_assembly/align_to_concatenated_db/{sample}/{sample}_indexed_contig.1.bt2",
        "results/{genera}/6_binning/semibin2/megahit_individual_assembly/align_to_concatenated_db/{sample}/{sample}_indexed_contig.2.bt2",
        "results/{genera}/6_binning/semibin2/megahit_individual_assembly/align_to_concatenated_db/{sample}/{sample}_indexed_contig.3.bt2",
        "results/{genera}/6_binning/semibin2/megahit_individual_assembly/align_to_concatenated_db/{sample}/{sample}_indexed_contig.4.bt2",
        "results/{genera}/6_binning/semibin2/megahit_individual_assembly/align_to_concatenated_db/{sample}/{sample}_indexed_contig.rev.1.bt2",
        "results/{genera}/6_binning/semibin2/megahit_individual_assembly/align_to_concatenated_db/{sample}/{sample}_indexed_contig.rev.2.bt2",
        "results/{genera}/6_binning/semibin2/megahit_individual_assembly/align_to_concatenated_db/{sample}_aligned_sorted.bam"
    params:
        genera=config["genera"],
        outdir = "results/{genera}/6_binning/semibin2/megahit_individual_assembly/align_to_concatenated_db/{sample}"
    log:
        stdout = "logs/{genera}/6_binning/semibin2/megahit_individual_assembly/align_to_concatenated_db/{sample}_aln.out",
        stderr = "logs/{genera}/6_binning/semibin2/megahit_individual_assembly/align_to_concatenated_db/{sample}_aln.err"
    shell:
        """
        module unload miniconda 
        module load Bowtie2/2.5.1-GCC-12.2.0
        module load SAMtools/1.21-GCC-12.2.0

        # 1. Build db index
        bowtie2-build \
        -f {input.contigs} {params.outdir}/{wildcards.sample}_indexed_contig \
        1>> {log.stdout} 2>> {log.stderr}

        # 2. Align reads back to SemiBin db index
        bowtie2 \
        -x {params.outdir}/{wildcards.sample}_indexed_contig -1 {input.r1} -2 {input.r2} | samtools view -b -F 4 -F 2048 | samtools sort -o {output[6]} \
        1>> {log.stdout} 2>> {log.stderr}
        """

rule semibin2_features_and_model_spades: # test
    """
    Generate sequence features and train model for SemiBin2
    """
    input:
        cat_fa = "results/{genera}/6_binning/semibin2/SPAdes_individual_assembly/generate_concatenated_db/concatenated.fa",
        bams = expand("results/{genera}/6_binning/semibin2/SPAdes_individual_assembly/align_to_concatenated_db/{sample}_aligned_sorted.bam", genera=config["genera"], sample=SAMPLES)
    output:
        split = "results/{genera}/6_binning/semibin2/SPAdes_individual_assembly/features_and_model/{sample}/data_split.csv",
        csv = "results/{genera}/6_binning/semibin2/SPAdes_individual_assembly/features_and_model/{sample}/data.csv"
    params:
        outdir = "results/{genera}/6_binning/semibin2/SPAdes_individual_assembly/features_and_model/{sample}",
        threads = 4
    log:
        stdout = "logs/{genera}/6_binning/semibin2/SPAdes_individual_assembly/features_and_model/{sample}/sequence_features.out",
        stderr = "logs/{genera}/6_binning/semibin2/SPAdes_individual_assembly/features_and_model/{sample}/sequence_features.err"
    shell:
        """
        module unload miniconda
        source activate /vast/palmer/pi/turner/flg9/conda_envs/semibin

        # Generate sequence features data.csv & data_split.csv files
        SemiBin2 generate_sequence_features_single \
        -i {input.cat_fa} \
        -b {input.bams} \
        -o {params.outdir} \
        -t {params.threads} \
        1>> {log.stdout} 2>> {log.stderr}
        """

rule semibin2_features_and_model_megahit: # test
    """
    Generate sequence features and train model for SemiBin2
    """
    input:
        cat_fa = "results/{genera}/6_binning/semibin2/megahit_individual_assembly/generate_concatenated_db/concatenated.fa",
        bams = expand("results/{genera}/6_binning/semibin2/megahit_individual_assembly/align_to_concatenated_db/{sample}_aligned_sorted.bam", genera=config["genera"], sample=SAMPLES)
    output:
        split = "results/{genera}/6_binning/semibin2/megahit_individual_assembly/features_and_model/{sample}/data_split.csv",
        csv = "results/{genera}/6_binning/semibin2/megahit_individual_assembly/features_and_model/{sample}/data.csv"
    params:
        outdir = "results/{genera}/6_binning/semibin2/megahit_individual_assembly/features_and_model/{sample}",
        threads = 4
    log:
        stdout = "logs/{genera}/6_binning/semibin2/megahit_individual_assembly/features_and_model/{sample}/sequence_features.out",
        stderr = "logs/{genera}/6_binning/semibin2/megahit_individual_assembly/features_and_model/{sample}/sequence_features.err"
    shell:
        """
        module unload miniconda
        source activate /vast/palmer/pi/turner/flg9/conda_envs/semibin

        # Generate sequence features data.csv & data_split.csv files
        SemiBin2 generate_sequence_features_single \
        -i {input.cat_fa} \
        -b {input.bams} \
        -o {params.outdir} \
        -t {params.threads} \
        1>> {log.stdout} 2>> {log.stderr}
        """

rule semibin2_train_model_spades: # test
    """
    Train ML model on previously curated SemiBin2 feature data
    """
    input:
        split = "results/{genera}/6_binning/semibin2/SPAdes_individual_assembly/features_and_model/{sample}/data_split.csv",
        csv = "results/{genera}/6_binning/semibin2/SPAdes_individual_assembly/features_and_model/{sample}/data.csv"
    output:
        model = "results/{genera}/6_binning/semibin2/SPAdes_individual_assembly/train_model/{sample}/model.pt"
    params:
        outdir = "results/{genera}/6_binning/semibin2/SPAdes_individual_assembly/train_model/{sample}",
        threads = 4
    log:
        stdout = "logs/{genera}/6_binning/semibin2/SPAdes_individual_assembly/train_model/{sample}/ML_train.out",
        stderr = "logs/{genera}/6_binning/semibin2/SPAdes_individual_assembly/train_model/{sample}/ML_train.err"
    shell:
        """
        module unload miniconda
        source activate /vast/palmer/pi/turner/flg9/conda_envs/semibin

        # Train model
        SemiBin2 train_self \
        --data {input.csv} \
        --data-split {input.split} \
        -o {params.outdir} \
        -t {params.threads} \
        1>> {log.stdout} 2>> {log.stderr}
        """

rule semibin2_train_model_megahit: # test
    """
    Train ML model on previously curated SemiBin2 feature data
    """
    input:
        split = "results/{genera}/6_binning/semibin2/megahit_individual_assembly/features_and_model/{sample}/data_split.csv",
        csv = "results/{genera}/6_binning/semibin2/megahit_individual_assembly/features_and_model/{sample}/data.csv"
    output:
        model = "results/{genera}/6_binning/semibin2/megahit_individual_assembly/train_model/{sample}/model.pt"
    params:
        outdir = "results/{genera}/6_binning/semibin2/megahit_individual_assembly/train_model/{sample}",
        threads = 4
    log:
        stdout = "logs/{genera}/6_binning/semibin2/megahit_individual_assembly/train_model/{sample}/ML_train.out",
        stderr = "logs/{genera}/6_binning/semibin2/megahit_individual_assembly/train_model/{sample}/ML_train.err"
    shell:
        """
        module unload miniconda
        source activate /vast/palmer/pi/turner/flg9/conda_envs/semibin

        # Train model
        SemiBin2 train_self \
        --data {input.csv} \
        --data-split {input.split} \
        -o {params.outdir} \
        -t {params.threads} \
        1>> {log.stdout} 2>> {log.stderr}
        """

rule semibin2_bin_spades: # test
    """
    Bin contigs using SemiBin's multi-sample binning model for individual binning of samples.
    This method often returns the most bins and is most optimized for complex samples.    
    """
    input:
        contigs = "results/{genera}/3_dedup_contigs/SPAdes/individual_metagenome_assembly/{sample}/{sample}_DEDUP95.fasta",
        csv = "results/{genera}/6_binning/semibin2/SPAdes_individual_assembly/features_and_model/{sample}/data.csv",
        model = "results/{genera}/6_binning/semibin2/SPAdes_individual_assembly/train_model/{sample}/model.pt"
    output:
        bins = "results/{genera}/6_binning/SPAdes_individual_assembly/semibin2/bin/{sample}/bin.*.fa"
    params:
        outdir = "results/{genera}/6_binning/SPAdes_individual_assembly/semibin2/bin/{sample}",
        seq_type = "short_reads",
        GTDB_path = "/vast/palmer/pi/turner/data/db/gtdbtk-2.4.1",
        minlen = "1500",
        threads = 4
    log:
        stdout = "logs/{genera}/6_binning/SPAdes_individual_assembly/semibin2/bin/{sample}/bin.out",
        stderr = "logs/{genera}/6_binning/SPAdes_individual_assembly/semibin2/bin/{sample}/bin.err"
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

rule semibin2_bin_megahit: # test
    """
    Bin contigs using SemiBin's multi-sample binning model for individual binning of samples.
    This method often returns the most bins and is most optimized for complex samples.    
    """
    input:
        contigs = "results/{genera}/3_dedup_contigs/megahit/individual_metagenome_assembly/{sample}/{sample}_DEDUP95.fasta",
        csv = "results/{genera}/6_binning/semibin2/megahit_individual_assembly/features_and_model/{sample}/data.csv",
        model = "results/{genera}/6_binning/semibin2/megahit_individual_assembly/train_model/{sample}/model.pt"
    output:
        bins = "results/{genera}/6_binning/megahit_individual_assembly/semibin2/bin/{sample}/bin.*.fa"
    params:
        outdir = "results/{genera}/6_binning/megahit_individual_assembly/semibin2/bin/{sample}",
        seq_type = "short_reads",
        GTDB_path = "/vast/palmer/pi/turner/data/db/gtdbtk-2.4.1",
        minlen = "1500",
        threads = 4
    log:
        stdout = "logs/{genera}/6_binning/megahit_individual_assembly/semibin2/bin/{sample}/bin.out",
        stderr = "logs/{genera}/6_binning/megahit_individual_assembly/semibin2/bin/{sample}/bin.err"
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

rule DASTool_spades: # test
    """
    Calculate an optimized set of bins from multiple binning tools that were previously used
    """
    input:
        maxbin_contigs = "results/{genera}/6_binning/maxbin/SPAdes_individual_assembly/{sample}/MAXBIN.*.fa",
        concoct_csv = "results/{genera}/6_binning/concoct/SPAdes_individual_assembly/{sample}/concoct_output/clustering_merged.csv",
        metabat_contigs = "results/{genera}/6_binning/metabat/SPAdes_individual_assembly/{sample}/METABAT.*.fa",
        semibin_contigs = "results/{genera}/6_binning/SPAdes_individual_assembly/semibin2/bin/{sample}/bin.*.fa",
        contigs = "results/{genera}/3_dedup_contigs/SPAdes/individual_metagenome_assembly/{sample}/{sample}_DEDUP95.fasta"
    output:
        metabat_summ = "results/{genera}/6_binning/SPAdes_individual_assembly/aggregate_bins/{sample}/metabat_associations.tsv",
        maxbin_summ = "results/{genera}/6_binning/SPAdes_individual_assembly/aggregate_bins/{sample}/maxbin_associations.tsv",
        concoct_summ = "results/{genera}/6_binning/SPAdes_individual_assembly/aggregate_bins/{sample}/concoct_associations.tsv",
        semibin_summ = "results/{genera}/6_binning/SPAdes_individual_assembly/aggregate_bins/{sample}/semibin_associations.tsv",
        bins = "results/{genera}/6_binning/SPAdes_individual_assembly/aggregate_bins/{sample}/DASTOOL.*.fa"
    params:
        genera=config["genera"],
        threads = 4,
        engine = "diamond",
        names = "Metabat,Maxbin,Concoct,SemiBin",
        dastool_outdir = "results/{genera}/6_binning/SPAdes_individual_assembly/aggregate_bins/{sample}"
    log:
        stdout = "logs/{genera}/6_binning/SPAdes_individual_assembly/aggregate_bins/{sample}/das_tool.out",
        stderr = "logs/{genera}/6_binning/SPAdes_individual_assembly/aggregate_bins/{sample}/das_tool.err"
    shell:
        """
        module unload miniconda 
        source activate /vast/palmer/pi/turner/flg9/conda_envs/das_tool

        # 1. Generate reporter file of bins and bin quality as input for DasTool

        # First, Maxbin
        Fasta_to_Contigs2Bin.sh \
        -i {input.maxbin_contigs} -e fasta > {output.maxbin_summ} \
        1>> {log.stdout} 2>> {log.stderr}

        # Second, Concoct
        perl \
        -pe "s/,/\tconcoct./g;" {input.concoct_csv} > {output.concoct_summ} \
        1>> {log.stdout} 2>> {log.stderr}

        # Third, Metabat
        Fasta_to_Contigs2Bin.sh \
        -i {input.metabat_contigs} -e fasta > {output.metabat_summ} \
        1>> {log.stdout} 2>> {log.stderr}

        # Fourth, SemiBin2
        Fasta_to_Contigs2Bin.sh \
        -i {input.semibin_contigs} -e fasta > {output.semibin_summ} \
        1>> {log.stdout} 2>> {log.stderr}

        # 2. Bin de-replication with Das_Tool
        mkdir -p \
        {params.dastool_outdir}

        DAS_Tool \
        -i {output.metabat_summ},{output.maxbin_summ},{output.concoct_summ},{output.semibin_summ} \
        -l {params.names} -o {params.dastool_outdir} -c {input.contigs} \
        --write_bin_evals --write_bins --threads {params.threads} --debug \
        --write_unbinned --write_bins --search_engine {params.engine} \
        1>> {log.stdout} 2>> {log.stderr}
        """

rule DASTool_megahit: # test
    """
    Calculate an optimized set of bins from multiple binning tools that were previously used
    """
    input:
        maxbin_contigs = "results/{genera}/6_binning/maxbin/megahit_individual_assembly/{sample}/MAXBIN.*.fa",
        concoct_csv = "results/{genera}/6_binning/concoct/megahit_individual_assembly/{sample}/concoct_output/clustering_merged.csv",
        metabat_contigs = "results/{genera}/6_binning/metabat/megahit_individual_assembly/{sample}/METABAT.*.fa",
        semibin_contigs = "results/{genera}/6_binning/megahit_individual_assembly/semibin2/bin/{sample}/bin.*.fa",
        contigs = "results/{genera}/3_dedup_contigs/megahit/individual_metagenome_assembly/{sample}/{sample}_DEDUP95.fasta"
    output:
        metabat_summ = "results/{genera}/6_binning/megahit_individual_assembly/aggregate_bins/{sample}/metabat_associations.tsv",
        maxbin_summ = "results/{genera}/6_binning/megahit_individual_assembly/aggregate_bins/{sample}/maxbin_associations.tsv",
        concoct_summ = "results/{genera}/6_binning/megahit_individual_assembly/aggregate_bins/{sample}/concoct_associations.tsv",
        semibin_summ = "results/{genera}/6_binning/megahit_individual_assembly/aggregate_bins/{sample}/semibin_associations.tsv",
        bins = "results/{genera}/6_binning/megahit_individual_assembly/aggregate_bins/{sample}/DASTOOL.*.fa"
    params:
        genera=config["genera"],
        threads = 4,
        engine = "diamond",
        names = "Metabat,Maxbin,Concoct,SemiBin",
        dastool_outdir = "results/{genera}/6_binning/megahit_individual_assembly/aggregate_bins/{sample}"
    log:
        stdout = "logs/{genera}/6_binning/megahit_individual_assembly/aggregate_bins/{sample}/das_tool.out",
        stderr = "logs/{genera}/6_binning/megahit_individual_assembly/aggregate_bins/{sample}/das_tool.err"
    shell:
        """
        module unload miniconda 
        source activate /vast/palmer/pi/turner/flg9/conda_envs/das_tool

        # 1. Generate reporter file of bins and bin quality as input for DasTool

        # First, Maxbin
        Fasta_to_Contigs2Bin.sh \
        -i {input.maxbin_contigs} -e fasta > {output.maxbin_summ} \
        1>> {log.stdout} 2>> {log.stderr}

        # Second, Concoct
        perl \
        -pe "s/,/\tconcoct./g;" {input.concoct_csv} > {output.concoct_summ} \
        1>> {log.stdout} 2>> {log.stderr}

        # Third, Metabat
        Fasta_to_Contigs2Bin.sh \
        -i {input.metabat_contigs} -e fasta > {output.metabat_summ} \
        1>> {log.stdout} 2>> {log.stderr}

        # Fourth, SemiBin2
        Fasta_to_Contigs2Bin.sh \
        -i {input.semibin_contigs} -e fasta > {output.semibin_summ} \
        1>> {log.stdout} 2>> {log.stderr}

        # 2. Bin de-replication with Das_Tool
        mkdir -p \
        {params.dastool_outdir}

        DAS_Tool \
        -i {output.metabat_summ},{output.maxbin_summ},{output.concoct_summ},{output.semibin_summ} \
        -l {params.names} -o {params.dastool_outdir} -c {input.contigs} \
        --write_bin_evals --write_bins --threads {params.threads} --debug \
        --write_unbinned --write_bins --search_engine {params.engine} \
        1>> {log.stdout} 2>> {log.stderr}
        """

rule bin_quality_check_spades:
    input:
        bins=lambda wildcards: (
            sorted(glob.glob(f"results/{wildcards.genera}/6_binning/maxbin/SPAdes_individual_assembly/{wildcards.sample}/MAXBIN.*.fa")) +
            sorted(glob.glob(f"results/{wildcards.genera}/6_binning/metabat/SPAdes_individual_assembly/{wildcards.sample}/METABAT.*.fa")) +
            sorted(glob.glob(f"results/{wildcards.genera}/6_binning/concoct/SPAdes_individual_assembly/{wildcards.sample}/CONCOCT.*.fa")) +
            sorted(glob.glob(f"results/{wildcards.genera}/6_binning/SPAdes_individual_assembly/semibin2/bin/{wildcards.sample}/bin.*.fa")) +
            sorted(glob.glob(f"results/{wildcards.genera}/6_binning/SPAdes_individual_assembly/aggregate_bins/{wildcards.sample}/DASTOOL.*.fa"))
        )
    output:
        "results/{genera}/6_binning/SPAdes_individual_assembly/binning_qc/{sample}/quality_report.tsv"
    params:
        threads=4,
        genera=config["genera"],
        outdir="results/{genera}/6_binning/SPAdes_individual_assembly/binning_qc/{sample}"
    log:
        stdout = "logs/{genera}/6_binning/SPAdes_individual_assembly/binning_qc/{sample}/checkm2.out",
        stderr = "logs/{genera}/6_binning/SPAdes_individual_assembly/binning_qc/{sample}/checkm2.err"
    shell:
        """
        module unload miniconda
        source activate /home/flg9/.conda/envs/checkm2
        export CHECKM2DB="/vast/palmer/pi/turner/data/db/CheckM2/CheckM2_database"

        mkdir -p {params.outdir}

        checkm2 predict \
        --threads {params.threads} \
        --input {input.bins} \
        --output_directory {params.outdir} \
        1>> {log.stdout} 2>> {log.stderr}
        """

rule bin_quality_check_megahit:
    input:
        bins=lambda wildcards: (
            sorted(glob.glob(f"results/{wildcards.genera}/6_binning/maxbin/megahit_individual_assembly/{wildcards.sample}/MAXBIN.*.fa")) +
            sorted(glob.glob(f"results/{wildcards.genera}/6_binning/metabat/megahit_individual_assembly/{wildcards.sample}/METABAT.*.fa")) +
            sorted(glob.glob(f"results/{wildcards.genera}/6_binning/concoct/megahit_individual_assembly/{wildcards.sample}/CONCOCT.*.fa")) +
            sorted(glob.glob(f"results/{wildcards.genera}/6_binning/megahit_individual_assembly/semibin2/bin/{wildcards.sample}/bin.*.fa")) +
            sorted(glob.glob(f"results/{wildcards.genera}/6_binning/megahit_individual_assembly/aggregate_bins/{wildcards.sample}/DASTOOL.*.fa"))
        )
    output:
        "results/{genera}/6_binning/megahit_individual_assembly/binning_qc/{sample}/quality_report.tsv"
    params:
        threads=4,
        genera=config["genera"],
        outdir="results/{genera}/6_binning/megahit_individual_assembly/binning_qc/{sample}"
    log:
        stdout = "logs/{genera}/6_binning/megahit_individual_assembly/binning_qc/{sample}/checkm2.out",
        stderr = "logs/{genera}/6_binning/megahit_individual_assembly/binning_qc/{sample}/checkm2.err"
    shell:
        """
        module unload miniconda
        source activate /home/flg9/.conda/envs/checkm2
        export CHECKM2DB="/vast/palmer/pi/turner/data/db/CheckM2/CheckM2_database"

        mkdir -p {params.outdir}

        checkm2 predict \
        --threads {params.threads} \
        --input {input.bins} \
        --output_directory {params.outdir} \
        1>> {log.stdout} 2>> {log.stderr}
        """