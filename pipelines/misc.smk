rule deepurify: # conda install and db install
    """
    Refine converged bins further by removing contaminants with DeepPurify
    """
    input:
        dastool_bins=lambda wildcards: sorted(glob.glob(f"results/{genera}/2_binning/aggregate_bins/reporter_files/DASTOOL.*.fa"))
    output:
        "results/{genera}/deepurify/{sample}/MetaInfo.tsv",
        "results/{genera}/deepurify/{sample}/clean_bins.fa"
    params:
        threads=4,
        db = path/to/db,
        gpus = 4,
        outdir = "results/{genera}/deepurify/{sample}",
        tmp = "results/{genera}/deepurify/tmp",
        suffix = "results/{genera}/deepurify/{sample}/DEEPURIFY"
    log:
        stdout = "logs/{genera}/deepurify/{sample}/deepurify.out",
        stderr = "logs/{genera}/deepurify/{sample}/deepurify.err"
    shell:
        """
        module unload miniconda 
        source activate XXXXXX

        # Run Deepurify clean mode
        deepurify clean \
        -i {input.dastool_bins} -o {params.outdir} \
        -db {params.db} --gpu_num {params.gpus} --bin_suffix {params.suffix} \
        --each_gpu_threads {params.threads} \
        --temp_output_folder {params.tmp} \
        1>> {log.stdout} 2>> {log.stderr}
        """

rule MetaWRAP: #download db
    """
    Refine bins with MetaWrap, comparing with DasTool output
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
        source activate /vast/palmer/pi/turner/flg9/conda_envs/metawrap-env
        export PATH=/vast/palmer/pi/turner/flg9/metaWRAP/bin:$PATH

        # 1. Run bin refinement module
        metaWRAP bin_refinement \
        -o {params.refinement_outdir} -t {params.threads} \
        -A {output.bin1} -B {output.bin2} -C {output.bin3} \
        1>> {log.stdout} 2>> {log.stderr}

        # 2. Visualize communities in our newly refined bins
        metaWRAP blobology \
        -a {input.assembly} -t {params.threads} -o {params.blobology_outdir} \ 
        {input.fwd} {input.rev} \
        1>> {log.stdout} 2>> {log.stderr}

        # 3. Plot abundance of bins
        metaWRAP quant_bins \
        -b {params.refinement_outdir}/metawrap_bins -o {params.abund_outdir} \
        -a {input.assembly} {input.fwd} {input.rev} \
        1>> {log.stdout} 2>> {log.stderr}

        # 4. Re-assemble MAGs using newly populated bins
        metaWRAP reassemble_bins \
        -o {params.reassembly_outdir} -1 {input.fwd} -2 {input.rev} \
        -t {params.threads} -b {params.refinement_outdir}/metawrap_bins \
        1>> {log.stdout} 2>> {log.stderr}

        # 5. Determine the taxonomy of each refined bin
        metaWRAP classify_bins \
        -b {params.reassembly_outdir}/reassembled_bins \
        -o {params.classify_outdir} -t {params.threads} \
        1>> {log.stdout} 2>> {log.stderr}

        # 6. Functionally annotate refined bins
        metaWRAP annotate_bins \
        -o {params.annot_outdir} -t {params.threads} \
        -b {params.reassembly_outdir}/reassembled_bins \
        1>> {log.stdout} 2>> {log.stderr}
        """ 

rule GUNC: # work in progress
    """
    Visualize binning accuracy
    """
    input:
        deepurify_bins=lambda wildcards: sorted(glob.glob(f"results/{genera}/deepurify/{sample}/clean_bins.fa"))
    output:
        "results/{genera}/vis_bins/all_levels.tsv",
        diamond = "results/{genera}/vis_bins/{sample}/diamond_output"
    params:
        db_path = "",
        threads = 4,
        tmp = "results/{genera}/vis_bins/tmp",
        run_outdir = "results/{genera}/vis_bins",
        plot_outdir = "results/{genera}/vis_bins/{sample}"
    log:
        stdout = "logs/{genera}/vis_bins/{sample}/gunc.out",
        stderr = "logs/{genera}/vis_bins/{sample}/gunc.err"
    shell:
        """
        module unload miniconda
        source activate XXXXXX

        # First, run on our custom binning output

        # 1. Run chimerism analysis
        gunc run \
        -i {input} -r {params.db_path} \
        --threads {params.threads} --temp_dir {params.tmp} \
        --sensitive --out_dir {params.run_outdir} --detailed_output \
        --contig_taxonomy_output \
        1>> {log.stdout} 2>> {log.stderr}

        # 2. Produce plot
        gunc run \
        --diamond_file {output.diamond} \
        --outdir {params.plot_outdir}

        # Second, run on wrapper output to compare with our custom set 

        # 1. 

        # 2. 
        """