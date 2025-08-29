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




##### integrate into binning test pipeline 







rule DASTool: #done
    """
    Calculate an optimized set of bins from multiple binning tools that were previously used
    """
    input:
        maxbin_contigs = "results/{genera}/2_binning/metabat/{sample}/MAXBIN.*.fa",
        concoct_csv = "results/{genera}/2_binning/concoct/{sample}/concoct_output/clustering_merged.csv",
        metabat_contigs = "results/{genera}/2_binning/metabat/{sample}/METABAT.*.fa",
        semibin_contigs = "results/{genera}/2_binning/semibin2/output_recluster_bins/semibin2.*.fa",
        contigs = "results/{genera}/1_assembly/dedup_contigs/{sample}/{sample}_DEDUP95.fasta"
    output:
        metabat_summ = "results/{genera}/2_binning/aggregate_bins/reporter_files/metabat_associations.tsv",
        maxbin_summ = "results/{genera}/2_binning/aggregate_bins/reporter_files/maxbin_associations.tsv",
        concoct_summ = "results/{genera}/2_binning/aggregate_bins/reporter_files/concoct_associations.tsv",
        semibin_summ = "results/{genera}/2_binning/aggregate_bins/reporter_files/semibin_associations.tsv",
        "results/{genera}/2_binning/aggregate_bins/reporter_files/DASTOOL.*.fa"
    params:
        genera=config["genera"],
        threads = 4,
        engine = "diamond",
        names = {Metabat,Maxbin,Concoct,SemiBin},
        dastool_outdir = "results/{genera}/2_binning/aggregate_bins/{sample}/aggregate_bin_{sample}"
    log:
        stdout = "logs/{genera}/2_binning/aggregate_bins/{sample}/das_tool.out",
        stderr = "logs/{genera}/2_binning/aggregate_bins/{sample}/das_tool.err"
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
        -i {input.sembin_contigs} -e fasta > {output.semibin_summ} \
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

def get_bins(wildcards):
    return {
        "maxbin_bins": sorted(glob.glob(f"results/{wildcards.genera}/2_binning/maxbin/{wildcards.sample}/MAXBIN.*.fa")),
        "metabat_bins": sorted(glob.glob(f"results/{wildcards.genera}/2_binning/metabat/{wildcards.sample}/METABAT.*.fa")),
        "concoct_bins": sorted(glob.glob(f"results/{wildcards.genera}/2_binning/concoct/{wildcards.sample}/CONCOCT.*.fa")),
        "semibin_bins": sorted(glob.glob(f"results/{wildcards.genera}/2_binning/semibin2/{wildcards.sample}/semibin2.*.fa")),
    }

rule bin_quality_check:
    input:
        get_bins
    output:
        "results/{genera}/2_binning/binning_qc/{sample}/quality_report.tsv"
    params:
        threads=4,
        outdir="results/{genera}/2_binning/binning_qc/{sample}"
    log:
        stdout = "logs/{genera}/2_binning/binning_qc/{sample}/checkm2.out",
        stderr = "logs/{genera}/2_binning/binning_qc/{sample}/checkm2.err"
    shell:
        """
        module unload miniconda
        source activate /home/flg9/.conda/envs/checkm2
        export CHECKM2DB="/vast/palmer/pi/turner/data/db/CheckM2/CheckM2_database"

        mkdir -p {params.outdir}

        checkm2 predict \
        --threads {params.threads} \
        --input {input.maxbin_bins} {input.concoct_bins} {input.metabat_bins} {input.semibin_bins} \
        --output_directory {params.outdir} \
        1>> {log.stdout} 2>> {log.stderr}
        """



        
        
        expand("results/{genera}/2_binning/aggregate_bins/reporter_files/metabat_associations.tsv", genera=config["genera"], sample=SAMPLES),
        expand("results/{genera}/2_binning/aggregate_bins/reporter_files/maxbin_associations.tsv", genera=config["genera"], sample=SAMPLES),
        expand("results/{genera}/2_binning/aggregate_bins/reporter_files/concoct_associations.tsv", genera=config["genera"], sample=SAMPLES),
        expand("results/{genera}/2_binning/aggregate_bins/reporter_files/DASTOOL.*.fa", genera=config["genera"], sample=SAMPLES),
        expand("results/{genera}/2_binning/binning_qc/{sample}/quality_report.tsv",genera=config["genera"],sample=SAMPLES)