import os, glob

rule concoct_cut_spades:
    """
    Split contigs into sizes appropriate for binning
    """
    input:
        contigs = "results/{genera}/3_dedup_contigs/SPAdes/individual_metagenome_assembly/{sample}/{sample}_DEDUP95.fasta"
    output:
        contigs_bed  = "results/{genera}/6_binning/concoct/SPAdes_individual_assembly/{sample}/contigs_20k.bed",
        contigs_fa   = "results/{genera}/6_binning/concoct/SPAdes_individual_assembly/{sample}/contigs_20k.fa"
    params:
        outdir = "results/{genera}/6_binning/concoct/SPAdes_individual_assembly/{sample}",
        contig_len = 20000,
        overlap = 0
    log:
        stdout = "logs/{genera}/6_binning/concoct/SPAdes_individual_assembly/{sample}/cut.out",
        stderr = "logs/{genera}/6_binning/concoct/SPAdes_individual_assembly/{sample}/cut.err"
    shell:
        """
        module unload miniconda
        mkdir -p {params.outdir}

        apptainer exec containers/concoct-1.1.0.sif \
        cut_up_fasta.py {input.contigs} -c {params.contig_len} -o {params.overlap} --merge_last \
        -b {output.contigs_bed} \
        2>> {log.stderr} | tee {output.contigs_fa} >> {log.stdout}
        """

rule concoct_coverage_spades:
    """
    Generate table with coverage depth information per sample and subcontig
    """
    input:
        contigs_bed = "results/{genera}/6_binning/concoct/SPAdes_individual_assembly/{sample}/contigs_20k.bed",
        bams = "results/{genera}/4_align_reads_to_contigs/contig_read_alignment_individual_assemblies_spades/{sample}_aligned_sorted.bam"
    output:
        coverage_tsv = "results/{genera}/6_binning/concoct/SPAdes_individual_assembly/{sample}/coverage_table.tsv"
    log:
        stdout = "logs/{genera}/6_binning/concoct/SPAdes_individual_assembly/{sample}/coverage.out",
        stderr = "logs/{genera}/6_binning/concoct/SPAdes_individual_assembly/{sample}/coverage.err"
    shell:
        """
        module unload miniconda

        apptainer exec containers/concoct-1.1.0.sif \
        concoct_coverage_table.py {input.contigs_bed} {input.bams} \
        2>> {log.stderr} | tee {output.coverage_tsv} >> {log.stdout}
        """

rule concoct_bin_spades:
    """
    Bin using CONCOCT (checkpoint because number of bin FASTAs is dynamic)
    """
    input:
        contigs_fa   = "results/{genera}/6_binning/concoct/SPAdes_individual_assembly/{sample}/contigs_20k.fa",
        coverage_tsv = "results/{genera}/6_binning/concoct/SPAdes_individual_assembly/{sample}/coverage_table.tsv"
    output:
        clustering = "results/{genera}/6_binning/concoct/SPAdes_individual_assembly/{sample}/bins/clustering.csv"
    params:
        threads = 4,
        min_len = 1500,
        outdir  = "results/{genera}/6_binning/concoct/SPAdes_individual_assembly/{sample}/bins"
    log:
        stdout = "logs/{genera}/6_binning/concoct/SPAdes_individual_assembly/{sample}/binning.out",
        stderr = "logs/{genera}/6_binning/concoct/SPAdes_individual_assembly/{sample}/binning.err"
    shell:
        """
        module unload miniconda
        mkdir -p {params.outdir}

        apptainer exec containers/concoct-1.1.0.sif \
        concoct --composition_file {input.contigs_fa} \
        --coverage_file {input.coverage_tsv} \
        -b {params.outdir}/ -t {params.threads} -l {params.min_len} \
        1>> {log.stdout} 2>> {log.stderr}

        mv {params.outdir}/clustering_gt{params.min_len}.csv {output.clustering}
        """

rule concoct_merge_spades:
    """
    Merge subcontig clustering into original contig clustering
    """
    input:
        clustering = "results/{genera}/6_binning/concoct/SPAdes_individual_assembly/{sample}/bins/clustering.csv"
    output:
        csv = "results/{genera}/6_binning/concoct/SPAdes_individual_assembly/{sample}/bins/clustering_merged.csv"
    params:
    log:
        stdout = "logs/{genera}/6_binning/concoct/SPAdes_individual_assembly/{sample}/merge.out",
        stderr = "logs/{genera}/6_binning/concoct/SPAdes_individual_assembly/{sample}/merge.err"
    shell:
        """
        module unload miniconda

        apptainer exec containers/concoct-1.1.0.sif \
        merge_cutup_clustering.py {input.clustering} \
        2>> {log.stderr} | tee {output.csv} >> {log.stdout} 
        """

checkpoint concoct_extract_bins_spades:
    """
    Extract bins as individual FASTAs
    """
    input:
        contigs = "results/{genera}/3_dedup_contigs/SPAdes/individual_metagenome_assembly/{sample}/{sample}_DEDUP95.fasta",
        csv = "results/{genera}/6_binning/concoct/SPAdes_individual_assembly/{sample}/bins/clustering_merged.csv"
    output:
        binsdir = directory("results/{genera}/6_binning/concoct/SPAdes_individual_assembly/{sample}/fasta_bins"),
        manifest  = "results/{genera}/6_binning/concoct/SPAdes_individual_assembly/{sample}/fasta_bins/manifest.txt"
    params:
        outdir = "results/{genera}/6_binning/concoct/SPAdes_individual_assembly/{sample}/fasta_bins",
    log:
        stdout = "logs/{genera}/6_binning/concoct/SPAdes_individual_assembly/{sample}/extract_fa.out",
        stderr = "logs/{genera}/6_binning/concoct/SPAdes_individual_assembly/{sample}/extract_fa.err"
    shell:
        """
        module unload miniconda 
        mkdir -p {params.outdir}

        apptainer exec containers/concoct-1.1.0.sif \
        extract_fasta_bins.py {input.contigs} {input.csv} --output_path {params.outdir} \
        1>> {log.stdout} 2>> {log.stderr}

        ls {params.outdir}/*.fa > {output.manifest}
        """

def get_concoct_bins_spades(wc):
    ckpt = checkpoints.concoct_extract_bins_spades.get(genera=wc.genera, sample=wc.sample)
    import os, glob
    bins_dir = ckpt.output.binsdir
    return sorted(glob.glob(os.path.join(bins_dir, "*.fa")))

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
        directory("results/{genera}/6_binning/metabat/SPAdes_individual_assembly/{sample}/bins"),
        check = "results/{genera}/6_binning/metabat/SPAdes_individual_assembly/{sample}/bins/check.txt"
    params:
        threads = 4,
        min_size = 1500,
        basename = "results/{genera}/6_binning/metabat/SPAdes_individual_assembly/{sample}/bins/METABAT",
        outdir="results/{genera}/6_binning/metabat/SPAdes_individual_assembly/{sample}/bins"
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
        -o {params.basename} \
        --verbose --debug \
        1>> {log.stdout} 2>> {log.stderr}

        ls {params.outdir}/*.fa > {output.check}
        """

def get_metabat_bins_spades(wc):
    ckpt = checkpoints.metabat2_bin_spades.get(genera=wc.genera, sample=wc.sample)
    import os, glob
    bins_dir = ckpt.output.outdir
    return sorted(glob.glob(os.path.join(bins_dir, "*.fa")))

rule maxbin2_depth_spades: # test
    """
    Get depth file for MaxBin using previously generated read alignment to contigs
    """
    input:
        metabat_txt = "results/{genera}/6_binning/metabat/SPAdes_individual_assembly/{sample}/METABAT.txt"
    output:
        depth_file = "results/{genera}/6_binning/maxbin/SPAdes_individual_assembly/{sample}/maxbin.txt"
    log:
        stdout = "logs/{genera}/6_binning/maxbin/SPAdes_individual_assembly/{sample}/maxbin_depth.out",
        stderr = "logs/{genera}/6_binning/maxbin/SPAdes_individual_assembly/{sample}/maxbin_depth.err"
    shell:
        """
        cut \
        -f1,4,6,8,10 {input.metabat_txt} \
        2>> {log.stderr} | tee {output.depth_file} 1>> {log.stdout}
        """

checkpoint maxbin2_bin_spades: # test
    """
    Bin contigs using MaxBin and previously generated depth file
    """
    input:
        contigs = "results/{genera}/3_dedup_contigs/SPAdes/individual_metagenome_assembly/{sample}/{sample}_DEDUP95.fasta",
        r1 = "results/{genera}/1_pre_processing/dedup_reads/{sample}/{sample}_host_removed_dedup_R1.fastq",
        r2 = "results/{genera}/1_pre_processing/dedup_reads/{sample}/{sample}_host_removed_dedup_R2.fastq",
        maxbin_depth_file = "results/{genera}/6_binning/maxbin/SPAdes_individual_assembly/{sample}/maxbin.txt"
    output:
        directory("results/{genera}/6_binning/maxbin/SPAdes_individual_assembly/{sample}/bins"),
        check = "results/{genera}/6_binning/maxbin/SPAdes_individual_assembly/{sample}/bins/check.txt"
    params:
        threads=4,
        contig_len = 1500,
        basename = "results/{genera}/6_binning/maxbin/SPAdes_individual_assembly/{sample}/bins/MAXBIN",
        outdir = "results/{genera}/6_binning/maxbin/SPAdes_individual_assembly/{sample}/bins"
    log:
        stdout = "logs/{genera}/6_binning/maxbin/SPAdes_individual_assembly/{sample}/maxbin.out",
        stderr = "logs/{genera}/6_binning/maxbin/SPAdes_individual_assembly/{sample}/maxbin.err"
    shell:
        """
        module unload miniconda
        module load MaxBin/2.2.7-gompi-2020b 

        mkdir -p {params.outdir}

        run_MaxBin.pl \
        -thread {params.threads} -min_contig_length {params.contig_len} \
        -contig {input.contigs} -reads {input.r1} -reads2 {input.r2} \
        -abund {input.maxbin_depth_file} \
        -out {params.basename} \
        1>> {log.stdout} 2>> {log.stderr}

        ls {params.outdir}/*.fasta > {output.check}
        """

def get_maxbin_bins_spades(wc):
    ckpt = checkpoints.maxbin2_bin_spades.get(genera=wc.genera, sample=wc.sample)
    import os, glob
    bins_dir = ckpt.output.outdir
    return sorted(glob.glob(os.path.join(bins_dir, "*.fasta")))

checkpoint easy_single_bin:
    """
    Run Semibin2 on easy single binning mode
    """
    input:
        fa = "results/{genera}/3_dedup_contigs/SPAdes/individual_metagenome_assembly/{sample}/{sample}_DEDUP95.fasta",
        bams = "results/{genera}/4_align_reads_to_contigs/contig_read_alignment_individual_assemblies_spades/{sample}_aligned_sorted.bam"
    output:
        outdir = directory("results/{genera}/6_binning/semibin2/SPAdes_individual_assembly/binning/{sample}/output_bins"),
        check = "results/{genera}/6_binning/semibin2/SPAdes_individual_assembly/binning/{sample}/check.tsv"
    params:
        base = "results/{genera}/6_binning/semibin2/SPAdes_individual_assembly/binning/{sample}",
        outdir = "results/{genera}/6_binning/semibin2/SPAdes_individual_assembly/binning/{sample}/output_bins",
        compress = "none",
        minlen = 1500,
        threads = 1
    log:
        stdout = "logs/{genera}/6_binning/semibin2/SPAdes_individual_assembly/binning/{sample}/single_binning.out",
        stderr = "logs/{genera}/6_binning/semibin2/SPAdes_individual_assembly/binning/{sample}/single_binning.err"
    shell:
        """
        module unload miniconda 
        source activate /vast/palmer/pi/turner/flg9/conda_envs/semibin

        SemiBin2 single_easy_bin \
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

def get_semibin2_bin_spades(wc):
    ckpt = checkpoints.easy_single_bin.get(genera=wc.genera, sample=wc.sample)
    import os, glob
    bins_dir = ckpt.output.outdir
    return sorted(glob.glob(os.path.join(bins_dir, "*.fa")))