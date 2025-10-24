    
checkpoint concoct_bins_megahit: #test
    """
    Group assembled contigs into bins that represent individual genomes or closely related organisms using Concoct
    """
    input:
        contigs = "results/{genera}/3_dedup_contigs/megahit/individual_metagenome_assembly/{sample}/{sample}_DEDUP95.fasta",
        bams = "results/{genera}/4_align_reads_to_contigs/contig_read_alignment_individual_assemblies_megahit/{sample}_aligned_sorted.bam",
        bai = "results/{genera}/4_align_reads_to_contigs/contig_read_alignment_individual_assemblies_megahit/{sample}_aligned_sorted.bam.bai"
    output:
        contigs_fa = "results/{genera}/6_binning/concoct/megahit_individual_assembly/{sample}/contigs_20k.fa",
        contigs_bed = "results/{genera}/6_binning/concoct/megahit_individual_assembly/{sample}/contigs_20k.bed",
        coverage_tsv = "results/{genera}/6_binning/concoct/megahit_individual_assembly/{sample}/coverage_table.tsv",
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
        stdout = "logs/{genera}/6_binning/concoct/megahit_individual_assembly/{sample}/concoct.out",
        stderr = "logs/{genera}/6_binning/concoct/megahit_individual_assembly/{sample}/concoct.err"
    shell:
        """
        module unload miniconda
        mkdir -p {params.outdir}

        # 1. Shred contigs into non-overlapping parts of equal length
        apptainer exec containers/concoct-1.1.0.sif \
        cut_up_fasta.py {input.contigs} -c {params.contig_len} --merge_last \
        -o {output.contigs_fa} -b {output.contigs_bed} \
        1>> {log.stdout} 2>> {log.stderr}

        # 2. Generate input coverage table for CONCOCT
        apptainer exec containers/concoct-1.1.0.sif \
        concoct_coverage_table.py {output.contigs_bed} {input.bams} > {output.coverage_tsv} \
        1>> {log.stdout} 2>> {log.stderr}
 
        # 3. Bin
        apptainer exec containers/concoct-1.1.0.sif \
        concoct --composition_file {output.contigs_fa} \
        --coverage_file {output.coverage_tsv} \
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

rule metabat2_depth_spades:
    input:
        bams = "results/{genera}/4_align_reads_to_contigs/contig_read_alignment_individual_assemblies_spades/{sample}_aligned_sorted.bam"
    output:
        depth_file = "results/{genera}/6_binning/metabat/SPAdes_individual_assembly/{sample}/METABAT.txt"
    params:
        dir = "results/{genera}/6_binning/metabat/SPAdes_individual_assembly/{sample}"
    log:
        stdout = "logs/{genera}/6_binning/metabat/SPAdes_individual_assembly/{sample}/depth.out",
        stderr = "logs/{genera}/6_binning/metabat/SPAdes_individual_assembly/{sample}/depth.err"
    shell:
        """
        mkdir -p {params.dir}
        echo "BAM file: {input.bams}"
        ls -lh {input.bams}
        echo "Output depth file: {output.depth_file}"

        apptainer exec --bind $PWD containers/metabat2-2.15.sif \
        jgi_summarize_bam_contig_depths \
        --outputDepth {output.depth_file} {input.bams} \
        1>> {log.stdout} 2>> {log.stderr}
        """

checkpoint metabat2_bin_spades:
    input:
        contigs = "results/{genera}/3_dedup_contigs/SPAdes/individual_metagenome_assembly/{sample}/{sample}_DEDUP95.fasta",
        depth_file = "results/{genera}/6_binning/metabat/SPAdes_individual_assembly/{sample}/METABAT.txt"
    output:
        directory("results/{genera}/6_binning/metabat/SPAdes_individual_assembly/{sample}")
    params:
        threads = 4,
        min_size = 1500,
        outdir = "results/{genera}/6_binning/metabat/SPAdes_individual_assembly/{sample}/METABAT"
    log:
        stdout = "logs/{genera}/6_binning/metabat/SPAdes_individual_assembly/{sample}/metabat.out",
        stderr = "logs/{genera}/6_binning/metabat/SPAdes_individual_assembly/{sample}/metabat.err"
    shell:
        """
        apptainer exec containers/metabat2-2.15.sif \
        metabat2 \
        -t {params.threads} -m {params.min_size} \
        -i {input.contigs} \
        -a {input.depth_file} \
        -o {params.outdir} \
        --verbose --debug \
        1>> {log.stdout} 2>> {log.stderr}
        """

def get_bins_spades(wildcards):
    ckpt = checkpoints.metabat2_bin_spades.get(genera=wildcards.genera, sample=wildcards.sample)
    # This resolves the checkpoint output directory
    outdir = ckpt.output[0]
    # Glob all .fa files inside
    bins = glob_wildcards(os.path.join(outdir, "METABAT.{i}.fa")).i
    return expand(os.path.join(outdir, "METABAT.{i}.fa"), i=bins)

rule metabat2_depth_megahit:
    input:
        bams = "results/{genera}/4_align_reads_to_contigs/contig_read_alignment_individual_assemblies_megahit/{sample}_aligned_sorted.bam"
    output:
        depth_file = "results/{genera}/6_binning/metabat/megahit_individual_assembly/{sample}/METABAT.txt"
    params:
        dir = "results/{genera}/6_binning/metabat/megahit_individual_assembly/{sample}"
    log:
        stdout = "logs/{genera}/6_binning/metabat/megahit_individual_assembly/{sample}/depth.out",
        stderr = "logs/{genera}/6_binning/metabat/megahit_individual_assembly/{sample}/depth.err"
    shell:
        """
        mkdir -p {params.dir}
        echo "BAM file: {input.bams}"
        ls -lh {input.bams}
        echo "Output depth file: {output.depth_file}"

        apptainer exec --bind $PWD containers/metabat2-2.15.sif \
        jgi_summarize_bam_contig_depths \
        --outputDepth {output.depth_file} {input.bams} \
        1>> {log.stdout} 2>> {log.stderr}
        """

checkpoint metabat2_bin_megahit:
    input:
        contigs = "results/{genera}/3_dedup_contigs/megahit/individual_metagenome_assembly/{sample}/{sample}_DEDUP95.fasta",
        depth_file = "results/{genera}/6_binning/metabat/megahit_individual_assembly/{sample}/METABAT.txt"
    output:
        directory("results/{genera}/6_binning/metabat/megahit_individual_assembly/{sample}")
    params:
        threads = 4,
        min_size = 1500,
        outdir = "results/{genera}/6_binning/metabat/megahit_individual_assembly/{sample}/METABAT"
    log:
        stdout = "logs/{genera}/6_binning/metabat/megahit_individual_assembly/{sample}/metabat.out",
        stderr = "logs/{genera}/6_binning/metabat/megahit_individual_assembly/{sample}/metabat.err"
    shell:
        """
        apptainer exec containers/metabat2-2.15.sif \
        metabat2 \
        -t {params.threads} -m {params.min_size} \
        -i {input.contigs} \
        -a {input.depth_file} \
        -o {params.outdir} \
        --verbose --debug \
        1>> {log.stdout} 2>> {log.stderr}
        """

def get_bins_megahit(wildcards):
    ckpt = checkpoints.metabat2_bin_megahit.get(genera=wildcards.genera, sample=wildcards.sample)
    # This resolves the checkpoint output directory
    outdir = ckpt.output[0]
    # Glob all .fa files inside
    bins = glob_wildcards(os.path.join(outdir, "METABAT.{i}.fa")).i
    return expand(os.path.join(outdir, "METABAT.{i}.fa"), i=bins)

rule maxbin2_depth_spades: # test
    """
    Get depth file for MaxBin using previously generated read alignment to contigs
    """
    input:
        metabat_txt = "results/{genera}/6_binning/metabat/SPAdes_individual_assembly/{sample}/METABAT.txt"
    output:
        depth_file = "results/{genera}/6_binning/maxbin/SPAdes_individual_assembly/{sample}/maxbin.txt"
    log:
        stdout = "logs/{genera}/6_binning/maxbin/SPAdes_individual_assembly/{sample}/maxbin.out",
        stderr = "logs/{genera}/6_binning/maxbin/SPAdes_individual_assembly/{sample}/maxbin.err"
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
        depth_file = "results/{genera}/6_binning/maxbin/megahit_individual_assembly/{sample}/maxbin.txt"
    log:
        stdout = "logs/{genera}/6_binning/maxbin/megahit_individual_assembly/{sample}/maxbin.out",
        stderr = "logs/{genera}/6_binning/maxbin/megahit_individual_assembly/{sample}/maxbin.err"
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
        maxbin_depth_file = "results/{genera}/6_binning/maxbin/SPAdes_individual_assembly/{sample}/maxbin.txt"
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
        maxbin_depth_file = "results/{genera}/6_binning/maxbin/megahit_individual_assembly/{sample}/maxbin.txt"
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




# working semibin process

rule semibin2_generate_concatenated_db_spades:
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

rule sembin2_align_to_concatenated_db_spades:
    """
    Align reads from each sample to our concatenated FASTA db, necessary for SemiBin pipeline
    """
    input:
        db = "results/{genera}/6_binning/semibin2/SPAdes_individual_assembly/generate_concatenated_db/concatenated.fa",
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
        -f {input.db} {params.outdir}/{wildcards.sample}_indexed_contig \
        1>> {log.stdout} 2>> {log.stderr}

        # 2. Align reads back to SemiBin db index
        bowtie2 \
        -x {params.outdir}/{wildcards.sample}_indexed_contig -1 {input.r1} -2 {input.r2} | samtools view -b -F 4 -F 2048 | samtools sort -o {output[6]} \
        1>> {log.stdout} 2>> {log.stderr}
        """

rule easy_multi_bin:
    """
    Run Semibin2 on easy multi binning mode
    """
    input:
        fa = "results/{genera}/6_binning/semibin2/SPAdes_individual_assembly/generate_concatenated_db/concatenated.fa",
        bams = expand("results/{genera}/6_binning/semibin2/SPAdes_individual_assembly/align_to_concatenated_db/{sample}_aligned_sorted.bam", genera=config["genera"], sample=SAMPLES)
    output:
        outdir = directory("results/{genera}/6_binning/semibin2/SPAdes_individual_assembly/binning/output_bins"),
        check = "results/{genera}/6_binning/semibin2/SPAdes_individual_assembly/binning/check.tsv"
    params:
        base = "results/{genera}/6_binning/semibin2/SPAdes_individual_assembly/binning/",
        outdir = "results/{genera}/6_binning/semibin2/SPAdes_individual_assembly/binning/output_bins",
        compress = "none",
        minlen = 1500,
        threads = 1
    log:
        stdout = "logs/{genera}/6_binning/semibin2/SPAdes_individual_assembly/binning/multi_binning.out",
        stderr = "logs/{genera}/6_binning/semibin2/SPAdes_individual_assembly/binning/multi_binning.err"
    shell:
        """
        module unload miniconda 
        source activate /vast/palmer/pi/turner/flg9/conda_envs/semibin

        SemiBin2 multi_easy_bin \
        -i {input.fa} \
        -b {input.bams} \
        --self-supervised \
        --min-len {params.minlen} \
        --threads {params.threads} \
        -o {params.base} \
        --compression {params.compress} \
        1>> {log.stdout} 2>> {log.stderr}

        ls {params.outdir}/*.fa > {output.check}
        """

rule semibin2_features_and_model_spades:
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
        outdir = "results/{genera}/6_binning/semibin2/SPAdes_individual_assembly/features_and_model",
        outdir2 =  "results/{genera}/6_binning/semibin2/SPAdes_individual_assembly/features_and_model/{sample}",
        threads = 4
    log:
        stdout = "logs/{genera}/6_binning/semibin2/SPAdes_individual_assembly/features_and_model/{sample}/sequence_features.out",
        stderr = "logs/{genera}/6_binning/semibin2/SPAdes_individual_assembly/features_and_model/{sample}/sequence_features.err"
    shell:
        """
        module unload miniconda
        source activate /vast/palmer/pi/turner/flg9/conda_envs/semibin

        # Generate sequence features data.csv & data_split.csv files
        SemiBin2 generate_sequence_features_multi \
        -i {input.cat_fa} \
        -b {input.bams} \
        -o {params.outdir} \
        -t {params.threads} \
        1>> {log.stdout} 2>> {log.stderr}

        mv {params.outdir}/samples/{wildcards.sample}_DEDUP95/data.csv {params.outdir2}
        mv {params.outdir}/samples/{wildcards.sample}_DEDUP95/data_split.csv {params.outdir2}
        """

rule semibin2_train_model_spades:
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

checkpoint semibin2_bin_spades:
    """
    Bin contigs using SemiBin's multi-sample binning model for individual binning of samples.
    This method often returns the most complete bins and is most optimized for complex samples.    
    """
    input:
        contigs = "results/{genera}/3_dedup_contigs/SPAdes/individual_metagenome_assembly/{sample}/{sample}_DEDUP95.fasta",
        csv = "results/{genera}/6_binning/semibin2/SPAdes_individual_assembly/features_and_model/{sample}/data.csv",
        model = "results/{genera}/6_binning/semibin2/SPAdes_individual_assembly/train_model/{sample}/model.pt"
    output:
        outdir = directory("results/{genera}/6_binning/semibin2/SPAdes_individual_assembly/binning/{sample}/output_bins"),
        check = "results/{genera}/6_binning/semibin2/SPAdes_individual_assembly/binning/{sample}/contig_bins.tsv"
    params:
        base = "results/{genera}/6_binning/semibin2/SPAdes_individual_assembly/binning/{sample}",
        outdir = "results/{genera}/6_binning/semibin2/SPAdes_individual_assembly/binning/{sample}/output_bins",
        minlen = 1500,
        threads = 4
    log:
        stdout = "logs/{genera}/6_binning/semibin2/SPAdes_individual_assembly/binning/{sample}/bin.out",
        stderr = "logs/{genera}/6_binning/semibin2/SPAdes_individual_assembly/binning/{sample}/bin.err"
    shell:
        """
        module unload miniconda
        source activate /vast/palmer/pi/turner/flg9/conda_envs/semibin

        # Bin
        SemiBin2 bin_short \
        -i {input.contigs} \
        --model {input.model} \
        --data {input.csv} \
        -o {params.base} \
        -t {params.threads} \
        --min-len={params.minlen} \
        1>> {log.stdout} 2>> {log.stderr}

        ls {params.outdir}/*.fa.gz > {output.check}
        """

def get_semibin2_bin_spades(wc):
    ckpt = checkpoints.semibin2_bin_spades.get(genera=wc.genera, sample=wc.sample)
    import os, glob
    bins_dir = ckpt.output.outdir
    return sorted(glob.glob(os.path.join(bins_dir, "*.fa.gz")))


