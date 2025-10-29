rule concoct_cut:
    """
    Split contigs into sizes appropriate for binning
    """
    input:
        contigs = "results/{genera}/1_metagenome_assembly/3_dedup_contigs/{assembler}/individual_metagenome_assembly/{sample}/{sample}_DEDUP95.fasta"
    output:
        contigs_bed  = "results/{genera}/1_metagenome_assembly/6_binning/concoct/{assembler}_individual_assembly/{sample}/contigs_20k.bed",
        contigs_fa   = "results/{genera}/1_metagenome_assembly/6_binning/concoct/{assembler}_individual_assembly/{sample}/contigs_20k.fa"
    params:
        outdir = "results/{genera}/1_metagenome_assembly/6_binning/concoct/{assembler}_individual_assembly/{sample}",
        contig_len = 20000,
        overlap = 0
    log:
        stdout = "logs/{genera}/1_metagenome_assembly/6_binning/concoct/{assembler}_individual_assembly/{sample}/cut.out",
        stderr = "logs/{genera}/1_metagenome_assembly/6_binning/concoct/{assembler}_individual_assembly/{sample}/cut.err"
    shell:
        """
        module unload miniconda
        mkdir -p {params.outdir}

        apptainer exec containers/concoct-1.1.0.sif \
        cut_up_fasta.py {input.contigs} \
        -c {params.contig_len} \
        -o {params.overlap} \
        --merge_last \
        -b {output.contigs_bed} \
        2>> {log.stderr} | tee {output.contigs_fa} >> {log.stdout}
        """

rule concoct_coverage:
    """
    Generate table with coverage depth information per sample and subcontig
    """
    input:
        contigs_bed = "results/{genera}/1_metagenome_assembly/6_binning/concoct/{assembler}_individual_assembly/{sample}/contigs_20k.bed",
        bams = "results/{genera}/1_metagenome_assembly/4_align_reads_to_contigs/{assembler}/contig_read_alignment_individual_assemblies/{sample}_aligned_sorted.bam"
    output:
        coverage_tsv = "results/{genera}/1_metagenome_assembly/6_binning/concoct/{assembler}_individual_assembly/{sample}/coverage_table.tsv"
    log:
        stdout = "logs/{genera}/1_metagenome_assembly/6_binning/concoct/{assembler}_individual_assembly/{sample}/coverage.out",
        stderr = "logs/{genera}/1_metagenome_assembly/6_binning/concoct/{assembler}_individual_assembly/{sample}/coverage.err"
    shell:
        """
        module unload miniconda

        apptainer exec containers/concoct-1.1.0.sif \
        concoct_coverage_table.py {input.contigs_bed} {input.bams} \
        2>> {log.stderr} | tee {output.coverage_tsv} >> {log.stdout}
        """

rule concoct_bin:
    """
    Bin using CONCOCT
    """
    input:
        contigs_fa   = "results/{genera}/1_metagenome_assembly/6_binning/concoct/{assembler}_individual_assembly/{sample}/contigs_20k.fa",
        coverage_tsv = "results/{genera}/1_metagenome_assembly/6_binning/concoct/{assembler}_individual_assembly/{sample}/coverage_table.tsv"
    output:
        clustering = "results/{genera}/1_metagenome_assembly/6_binning/concoct/{assembler}_individual_assembly/{sample}/bins/clustering.csv"
    params:
        threads = 4,
        min_len = 1500,
        outdir  = "results/{genera}/1_metagenome_assembly/6_binning/concoct/{assembler}_individual_assembly/{sample}/bins"
    log:
        stdout = "logs/{genera}/1_metagenome_assembly/6_binning/concoct/{assembler}_individual_assembly/{sample}/binning.out",
        stderr = "logs/{genera}/1_metagenome_assembly/6_binning/concoct/{assembler}_individual_assembly/{sample}/binning.err"
    shell:
        """
        module unload miniconda
        mkdir -p {params.outdir}

        apptainer exec containers/concoct-1.1.0.sif \
        concoct --composition_file {input.contigs_fa} \
        --coverage_file {input.coverage_tsv} \
        -b {params.outdir}/ \
        -t {params.threads} \
        -l {params.min_len} \
        1>> {log.stdout} 2>> {log.stderr}

        mv {params.outdir}/clustering_gt{params.min_len}.csv {output.clustering}
        """

rule concoct_merge:
    """
    Merge subcontig clustering into original contig clustering
    """
    input:
        clustering = "results/{genera}/1_metagenome_assembly/6_binning/concoct/{assembler}_individual_assembly/{sample}/bins/clustering.csv"
    output:
        csv = "results/{genera}/1_metagenome_assembly/6_binning/concoct/{assembler}_individual_assembly/{sample}/bins/clustering_merged.csv"
    params:
    log:
        stdout = "logs/{genera}/1_metagenome_assembly/6_binning/concoct/{assembler}_individual_assembly/{sample}/merge.out",
        stderr = "logs/{genera}/1_metagenome_assembly/6_binning/concoct/{assembler}_individual_assembly/{sample}/merge.err"
    shell:
        """
        module unload miniconda

        apptainer exec containers/concoct-1.1.0.sif \
        merge_cutup_clustering.py {input.clustering} \
        2>> {log.stderr} | tee {output.csv} >> {log.stdout} 
        """

checkpoint concoct_extract_bins:
    """
    Extract bins as individual FASTAs
    """
    input:
        contigs = "results/{genera}/1_metagenome_assembly/3_dedup_contigs/{assembler}/individual_metagenome_assembly/{sample}/{sample}_DEDUP95.fasta",
        csv = "results/{genera}/1_metagenome_assembly/6_binning/concoct/{assembler}_individual_assembly/{sample}/bins/clustering_merged.csv"
    output:
        binsdir = directory("results/{genera}/1_metagenome_assembly/6_binning/concoct/{assembler}_individual_assembly/{sample}/fasta_bins"),
        manifest  = "results/{genera}/1_metagenome_assembly/6_binning/concoct/{assembler}_individual_assembly/{sample}/fasta_bins/manifest.txt"
    params:
        outdir = "results/{genera}/1_metagenome_assembly/6_binning/concoct/{assembler}_individual_assembly/{sample}/fasta_bins",
    log:
        stdout = "logs/{genera}/1_metagenome_assembly/6_binning/concoct/{assembler}_individual_assembly/{sample}/extract_fa.out",
        stderr = "logs/{genera}/1_metagenome_assembly/6_binning/concoct/{assembler}_individual_assembly/{sample}/extract_fa.err"
    shell:
        """
        module unload miniconda 
        mkdir -p {params.outdir}

        apptainer exec containers/concoct-1.1.0.sif \
        extract_fasta_bins.py {input.contigs} {input.csv} \
        --output_path {params.outdir} \
        1>> {log.stdout} 2>> {log.stderr}

        ls {params.outdir}/*.fa > {output.manifest}
        """

def get_concoct_bins(wc):
    ckpt = checkpoints.concoct_extract_bins.get(genera=wc.genera, sample=wc.sample)
    import os, glob
    bins_dir = ckpt.output.binsdir
    return sorted(glob.glob(os.path.join(bins_dir, "*.fa")))

rule metabat2_depth:
    input:
        bams = "results/{genera}/1_metagenome_assembly/4_align_reads_to_contigs/{assembler}/contig_read_alignment_individual_assemblies/{sample}_aligned_sorted.bam"
    output:
        depth_file = "results/{genera}/1_metagenome_assembly/6_binning/metabat/{assembler}_individual_assembly/{sample}/METABAT.txt"
    params:
        dir = "results/{genera}/1_metagenome_assembly/6_binning/metabat/{assembler}_individual_assembly/{sample}"
    log:
        stdout = "logs/{genera}/1_metagenome_assembly/6_binning/metabat/{assembler}_individual_assembly/{sample}/depth.out",
        stderr = "logs/{genera}/1_metagenome_assembly/6_binning/metabat/{assembler}_individual_assembly/{sample}/depth.err"
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

checkpoint metabat2_bin:
    input:
        contigs = "results/{genera}/1_metagenome_assembly/3_dedup_contigs/{assembler}/individual_metagenome_assembly/{sample}/{sample}_DEDUP95.fasta",
        depth_file = "results/{genera}/1_metagenome_assembly/6_binning/metabat/{assembler}_individual_assembly/{sample}/METABAT.txt"
    output:
        directory("results/{genera}/1_metagenome_assembly/6_binning/metabat/{assembler}_individual_assembly/{sample}/bins"),
        check = "results/{genera}/1_metagenome_assembly/6_binning/metabat/{assembler}_individual_assembly/{sample}/bins/check.txt"
    params:
        threads = 4,
        min_size = 1500,
        basename = "results/{genera}/1_metagenome_assembly/6_binning/metabat/{assembler}_individual_assembly/{sample}/bins/METABAT",
        outdir="results/{genera}/1_metagenome_assembly/6_binning/metabat/{assembler}_individual_assembly/{sample}/bins"
    log:
        stdout = "logs/{genera}/1_metagenome_assembly/6_binning/metabat/{assembler}_individual_assembly/{sample}/metabat.out",
        stderr = "logs/{genera}/1_metagenome_assembly/6_binning/metabat/{assembler}_individual_assembly/{sample}/metabat.err"
    shell:
        """
        apptainer exec containers/metabat2-2.15.sif \
        metabat2 \
        -t {params.threads} \
        -m {params.min_size} \
        -i {input.contigs} \
        -a {input.depth_file} \
        -o {params.basename} \
        --verbose --debug \
        1>> {log.stdout} 2>> {log.stderr}

        ls {params.outdir}/*.fa > {output.check}
        """

def get_metabat_bins_spades(wc):
    ckpt = checkpoints.metabat2_bin.get(genera=wc.genera, sample=wc.sample)
    import os, glob
    bins_dir = ckpt.output.outdir
    return sorted(glob.glob(os.path.join(bins_dir, "*.fa")))

rule maxbin2_depth:
    """
    Get depth file for MaxBin using previously generated read alignment to contigs
    """
    input:
        metabat_txt = "results/{genera}/1_metagenome_assembly/6_binning/metabat/{assembler}_individual_assembly/{sample}/METABAT.txt"
    output:
        depth_file = "results/{genera}/1_metagenome_assembly/6_binning/maxbin/{assembler}_individual_assembly/{sample}/maxbin.txt"
    log:
        stdout = "logs/{genera}/1_metagenome_assembly/6_binning/maxbin/{assembler}_individual_assembly/{sample}/maxbin_depth.out",
        stderr = "logs/{genera}/1_metagenome_assembly/6_binning/maxbin/{assembler}_individual_assembly/{sample}/maxbin_depth.err"
    shell:
        """
        cut \
        -f1,4,6,8,10 {input.metabat_txt} \
        2>> {log.stderr} | tee {output.depth_file} 1>> {log.stdout}
        """

checkpoint maxbin2_bin:
    """
    Bin contigs using MaxBin and previously generated depth file
    """
    input:
        contigs = "results/{genera}/1_metagenome_assembly/3_dedup_contigs/{assembler}/individual_metagenome_assembly/{sample}/{sample}_DEDUP95.fasta",
        r1 = "results/{genera}/1_metagenome_assembly/1_pre_processing/dedup_reads/{sample}/{sample}_host_removed_dedup_R1.fastq",
        r2 = "results/{genera}/1_metagenome_assembly/1_pre_processing/dedup_reads/{sample}/{sample}_host_removed_dedup_R2.fastq",
        maxbin_depth_file = "results/{genera}/1_metagenome_assembly/6_binning/maxbin/{assembler}_individual_assembly/{sample}/maxbin.txt"
    output:
        directory("results/{genera}/1_metagenome_assembly/6_binning/maxbin/{assembler}_individual_assembly/{sample}/bins"),
        check = "results/{genera}/1_metagenome_assembly/6_binning/maxbin/{assembler}_individual_assembly/{sample}/bins/check.txt"
    params:
        threads=4,
        contig_len = 1500,
        basename = "results/{genera}/1_metagenome_assembly/6_binning/maxbin/{assembler}_individual_assembly/{sample}/bins/MAXBIN",
        outdir = "results/{genera}/1_metagenome_assembly/6_binning/maxbin/{assembler}_individual_assembly/{sample}/bins"
    log:
        stdout = "logs/{genera}/1_metagenome_assembly/6_binning/maxbin/{assembler}_individual_assembly/{sample}/maxbin.out",
        stderr = "logs/{genera}/1_metagenome_assembly/6_binning/maxbin/{assembler}_individual_assembly/{sample}/maxbin.err"
    shell:
        """
        module unload miniconda
        module load MaxBin/2.2.7-gompi-2020b 

        mkdir -p {params.outdir}

        run_MaxBin.pl \
        -thread {params.threads} \
        -min_contig_length {params.contig_len} \
        -contig {input.contigs} \
        -reads {input.r1} -reads2 {input.r2} \
        -abund {input.maxbin_depth_file} \
        -out {params.basename} \
        1>> {log.stdout} 2>> {log.stderr}

        ls {params.outdir}/*.fasta > {output.check}
        """

def get_maxbin_bins(wc):
    ckpt = checkpoints.maxbin2_bin_spades.get(genera=wc.genera, sample=wc.sample)
    import os, glob
    bins_dir = ckpt.output.outdir
    return sorted(glob.glob(os.path.join(bins_dir, "*.fasta")))

checkpoint easy_single_bin:
    """
    Run Semibin2 on easy single binning mode
    """
    input:
        fa = "results/{genera}/1_metagenome_assembly/3_dedup_contigs/{assembler}/individual_metagenome_assembly/{sample}/{sample}_DEDUP95.fasta",
        bams = "results/{genera}/1_metagenome_assembly/4_align_reads_to_contigs/{assembler}/contig_read_alignment_individual_assemblies/{sample}_aligned_sorted.bam"
    output:
        outdir = directory("results/{genera}/1_metagenome_assembly/6_binning/semibin2/{assembler}_individual_assembly/single_binning/{sample}/output_bins"),
        check = "results/{genera}/1_metagenome_assembly/6_binning/semibin2/{assembler}_individual_assembly/single_binning/{sample}/check.tsv"
    params:
        base = "results/{genera}/1_metagenome_assembly/6_binning/semibin2/{assembler}_individual_assembly/single_binning/{sample}",
        outdir = "results/{genera}/1_metagenome_assembly/6_binning/semibin2/{assembler}_individual_assembly/single_binning/{sample}/output_bins",
        compress = "none",
        minlen = 1500,
        threads = 1
    log:
        stdout = "logs/{genera}/1_metagenome_assembly/6_binning/semibin2/{assembler}_individual_assembly/single_binning/{sample}/single_binning.out",
        stderr = "logs/{genera}/1_metagenome_assembly/6_binning/semibin2/{assembler}_individual_assembly/single_binning/{sample}/single_binning.err"
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

def get_semibin2_bin(wc):
    ckpt = checkpoints.easy_single_bin.get(genera=wc.genera, sample=wc.sample)
    import os, glob
    bins_dir = ckpt.output.outdir
    return sorted(glob.glob(os.path.join(bins_dir, "*.fa")))

rule prep_DASTool_input:
    """
    Create tab separated files of contig IDs and bin IDs required for DASTool input
    """
    input:
        semibin_bins = "results/{genera}/1_metagenome_assembly/6_binning/semibin2/{assembler}_individual_assembly/single_binning/{sample}/output_bins",
        concoct_bins = "results/{genera}/1_metagenome_assembly/6_binning/concoct/{assembler}_individual_assembly/{sample}/fasta_bins",
        maxbin_bins = "results/{genera}/1_metagenome_assembly/6_binning/maxbin/{assembler}_individual_assembly/{sample}/bins",
        metabat_bins = "results/{genera}/1_metagenome_assembly/6_binning/metabat/{assembler}_individual_assembly/{sample}/bins"
    output:
        semibin_prep = "results/{genera}/1_metagenome_assembly/6_binning/DASTool/{assembler}_individual_assembly/contigs2bin/{sample}.semibin.contigs2bin.tsv",
        concoct_prep = "results/{genera}/1_metagenome_assembly/6_binning/DASTool/{assembler}_individual_assembly/contigs2bin/{sample}.concoct.contigs2bin.tsv",
        maxbin_prep = "results/{genera}/1_metagenome_assembly/6_binning/DASTool/{assembler}_individual_assembly/contigs2bin/{sample}.maxbin.contigs2bin.tsv",
        metabat_prep = "results/{genera}/1_metagenome_assembly/6_binning/DASTool/{assembler}_individual_assembly/contigs2bin/{sample}.metabat.contigs2bin.tsv"
    params:
        ext1 = "fa",
        ext2 = "fasta"
    shell:
        """
        # SemiBin
        scripts/Fasta_to_Contig2Bin.sh \
        -i {input.semibin_bins} -e {params.ext1} > {output.semibin_prep}

        # Concoct
        scripts/Fasta_to_Contig2Bin.sh \
        -i {input.concoct_bins} -e {params.ext1} > {output.concoct_prep}

        # MaxBin
        scripts/Fasta_to_Contig2Bin.sh \
        -i {input.maxbin_bins} -e {params.ext2} > {output.maxbin_prep}

        # Metabat
        scripts/Fasta_to_Contig2Bin.sh \
        -i {input.metabat_bins} -e {params.ext1} > {output.metabat_prep}
        """

checkpoint DASTool_bin_refinement:
    """
    Refine bins using DASTool
    """
    input:
        contigs = "results/{genera}/1_metagenome_assembly/3_dedup_contigs/{assembler}/individual_metagenome_assembly/{sample}/{sample}_DEDUP95.fasta",
        input1 = "results/{genera}/1_metagenome_assembly/6_binning/DASTool/{assembler}_individual_assembly/contigs2bin/{sample}.semibin.contigs2bin.tsv",
        input2 = "results/{genera}/1_metagenome_assembly/6_binning/DASTool/{assembler}_individual_assembly/contigs2bin/{sample}.concoct.contigs2bin.tsv",
        input3 = "results/{genera}/1_metagenome_assembly/6_binning/DASTool/{assembler}_individual_assembly/contigs2bin/{sample}.maxbin.contigs2bin.tsv",
        input4 = "results/{genera}/1_metagenome_assembly/6_binning/DASTool/{assembler}_individual_assembly/contigs2bin/{sample}.metabat.contigs2bin.tsv"
    output:
        outdir = directory("results/{genera}/1_metagenome_assembly/6_binning/DASTool/{assembler}_individual_assembly/refined_bins/{sample}/_DASTool_bins"),
        check = "results/{genera}/1_metagenome_assembly/6_binning/DASTool/{assembler}_individual_assembly/refined_bins/{sample}/check.tsv"
    params:
        names = "SemiBin,Concoct,MaxBin,MetaBat",
        outdir = "results/{genera}/1_metagenome_assembly/6_binning/DASTool/{assembler}_individual_assembly/refined_bins/{sample}",
        bins_outdir = "results/{genera}/1_metagenome_assembly/6_binning/DASTool/{assembler}_individual_assembly/refined_bins/{sample}/_DASTool_bins",
        basename = "results/{genera}/1_metagenome_assembly/6_binning/DASTool/{assembler}_individual_assembly/refined_bins/{sample}/",
        engine = "diamond",
        threads = 2
    log:
        stdout = "logs/{genera}/1_metagenome_assembly/6_binning/DASTool/{assembler}_individual_assembly/refined_bins/{sample}/das_tool.out",
        stderr = "logs/{genera}/1_metagenome_assembly/6_binning/DASTool/{assembler}_individual_assembly/refined_bins/{sample}/das_tool.err"
    shell:
        """
        module unload miniconda
        source activate /vast/palmer/pi/turner/flg9/conda_envs/das_tool

        mkdir -p \
        {params.outdir}

        DAS_Tool \
        -i {input.input1},{input.input2},{input.input3},{input.input4} \
        -l {params.names} \
        -c {input.contigs} \
        -o {params.basename} \
        --search_engine {params.engine} \
        -t {params.threads} \
        --write_bin_evals \
        --write_bins \
        --write_unbinned \
        --debug \
        1>> {log.stdout} 2>> {log.stderr}

        ls {params.bins_outdir}/*.fa > {output.check}
        """

def get_DASTool_bins(wc):
    ckpt = checkpoints.DASTool_bin_refinement.get(genera=wc.genera, sample=wc.sample)
    import os, glob
    bins_dir = ckpt.output.outdir
    return sorted(glob.glob(os.path.join(bins_dir, "*.fa")))

rule CheckM:
    """
    Check bin quality with CheckM
    """
    input:
        semibin_bins = "results/{genera}/1_metagenome_assembly/6_binning/semibin2/{assembler}_individual_assembly/single_binning/{sample}/output_bins",
        concoct_bins = "results/{genera}/1_metagenome_assembly/6_binning/concoct/{assembler}_individual_assembly/{sample}/fasta_bins",
        maxbin_bins = "results/{genera}/1_metagenome_assembly/6_binning/maxbin/{assembler}_individual_assembly/{sample}/bins",
        metabat_bins = "results/{genera}/1_metagenome_assembly/6_binning/metabat/{assembler}_individual_assembly/{sample}/bins",
        dastool_bins = "results/{genera}/1_metagenome_assembly/6_binning/DASTool/{assembler}_individual_assembly/refined_bins/{sample}/_DASTool_bins"
    output:
        lineage = "results/{genera}/1_metagenome_assembly/6_binning/CheckM/{assembler}_individual_assembly/{sample}/lineage.ms",
        out1 = "results/{genera}/1_metagenome_assembly/6_binning/CheckM/{assembler}_individual_assembly/{sample}/lineage_results.tsv",
        out2 = "results/{genera}/1_metagenome_assembly/6_binning/CheckM/{assembler}_individual_assembly/{sample}/qa_results.tsv"
    params:
        threads = 1,
        ext = "fa",
        format = 2,
        outdir="results/{genera}/1_metagenome_assembly/6_binning/CheckM/{assembler}_individual_assembly/{sample}"
    log:
        stdout = "logs/{genera}/1_metagenome_assembly/6_binning/CheckM/{assembler}_individual_assembly/{sample}/CheckM.out",
        stderr = "logs/{genera}/1_metagenome_assembly/6_binning/CheckM/{assembler}_individual_assembly/{sample}/CheckM.err"
    shell:
        """
        module unload miniconda
        source activate /home/flg9/.conda/envs/CheckM
        checkm data setRoot /vast/palmer/pi/turner/data/db/CheckM/checkm_data_2015_01_16

        mkdir -p \
        {params.outdir}

        # 1. Assess bins using lineage-specific gene marker sets
        checkm lineage_wf \
        {input.dastool_bins} {params.outdir} \
        -t {params.threads} \
        -x {params.ext} \
        --tab_table \
        -f {output.out1} \
        1>> {log.stdout} 2>> {log.stderr}

        # 2. Produce a detailed report regarding completeness and contamination of our bins
        checkm qa \
        {output.lineage} {params.outdir} \
        -o {params.format} \
        -t {params.threads} \
        -f {output.out2} \
        --tab_table \
        1>> {log.stdout} 2>> {log.stderr}
        """

rule CheckM2:
    """
    Check the quality of our bins with CheckM2
    """
    input:
        semibin_bins = "results/{genera}/1_metagenome_assembly/6_binning/semibin2/{assembler}_individual_assembly/single_binning/{sample}/output_bins",
        concoct_bins = "results/{genera}/1_metagenome_assembly/6_binning/concoct/{assembler}_individual_assembly/{sample}/fasta_bins",
        maxbin_bins = "results/{genera}/1_metagenome_assembly/6_binning/maxbin/{assembler}_individual_assembly/{sample}/bins",
        metabat_bins = "results/{genera}/1_metagenome_assembly/6_binning/metabat/{assembler}_individual_assembly/{sample}/bins",
        dastool_bins = "results/{genera}/1_metagenome_assembly/6_binning/DASTool/{assembler}_individual_assembly/refined_bins/{sample}/_DASTool_bins"
    output:
        "results/{genera}/1_metagenome_assembly/6_binning/CheckM2/{assembler}_individual_assembly/{sample}/quality_report.tsv"
    params:
        outdir = "results/{genera}/1_metagenome_assembly/6_binning/CheckM2/{assembler}_individual_assembly/{sample}/",
        ext = ".fa",
        threads = 2
    log:
        stdout = "logs/{genera}/1_metagenome_assembly/6_binning/CheckM2/{assembler}_individual_assembly/{sample}/checkm2.out",
        stderr = "logs/{genera}/1_metagenome_assembly/6_binning/CheckM2/{assembler}_individual_assembly/{sample}/checkm2.err"
    shell:
        """
        module unload miniconda
        source activate /home/flg9/.conda/envs/checkm2
        export CHECKM2DB=/vast/palmer/pi/turner/data/db/CheckM2/CheckM2_database/uniref100.KO.1.dmnd

        mkdir -p \
        {params.outdir}

        checkm2 predict \
        --input {input.dastool_bins} \
        --output-directory {params.outdir}\
        --threads {params.threads} \
        -x {params.ext} \
        --allmodels \
        --force \
        --debug \
        1>> {log.stdout} 2>> {log.stderr}
        """

rule GUNC_run: 
    """
    Visualize chimerism and contamination of MAGs
    """
    input:
        dastool_bins = "results/{genera}/1_metagenome_assembly/6_binning/DASTool/{assembler}_individual_assembly/refined_bins/{sample}/_DASTool_bins"
    output:
        "results/{genera}/1_metagenome_assembly/6_binning/GUNC/{assembler}_individual_assembly/{sample}/check.txt"
    params:
        db = "/vast/palmer/pi/turner/flg9/TLab/Projects/Beluga_Metagenome_Assembly_Testing/db/gunc_db_progenomes2.1.dmnd",
        tmp = "results/{genera}/1_metagenome_assembly/6_binning/GUNC/{assembler}_individual_assembly/{sample}/tmp",
        outdir = "results/{genera}/1_metagenome_assembly/6_binning/GUNC/{assembler}_individual_assembly/{sample}/output",
        suffix = ".fa",
        threads = 1
    log:
        stdout = "logs/{genera}/1_metagenome_assembly/6_binning/GUNC/{assembler}_individual_assembly/{sample}/gunc.out",
        stderr = "logs/{genera}/1_metagenome_assembly/6_binning/GUNC/{assembler}_individual_assembly/{sample}/gunc.err"
    shell:
        """
        mkdir -p {params.outdir} {params.tmp}

        apptainer exec containers/gunc-1.0.6.sif \
        gunc run \
        --input_dir {input.dastool_bins} \
        --db_file {params.db} \
        --file_suffix {params.suffix} \
        --threads {params.threads} \
        --temp_dir {params.tmp} \
        --out_dir {params.outdir} \
        --detailed_output \
        --contig_taxonomy_output \
        1>> {log.stdout} 2>> {log.stderr}

        touch {output}
        """