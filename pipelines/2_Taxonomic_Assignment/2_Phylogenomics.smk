# This workflow is designed to phylogenetically place our MAGs across our time series

rule phylogenomics:
    """
    Place binned MAGs on reference-based phylogeny w/ nucleotide sequences
    """
    input:
        dir = "results/{genera}/1_metagenome_assembly/6_binning/DASTool/{assembler}_individual_assembly/refined_bins/{sample}/_DASTool_bins"
    output:
        "results/{genera}/2_Taxonomic_Assignment/2_phylogenomics/PhyloPhlAn/{assembler}/{sample}/RAxML_bestTree._DASTool_bins_refined.tre"
    params:
        cpus = 1,
        diversity = "high",
        ext = ".fa",
        db = "phylophlan", # this db uses 400 marker genes # amphora db is the other default option with 136 marker genes
        configfile = "/vast/palmer/pi/turner/flg9/conda_envs/PhyloPhlAn/configs/supermatrix_aa.cfg",
        outdir = "results/{genera}/2_Taxonomic_Assignment/2_phylogenomics/PhyloPhlAn/{assembler}/{sample}"
    log:
        stdout = "logs/{genera}/2_Taxonomic_Assignment/2_phylogenomics/PhyloPhlAn/{assembler}/{sample}/phylo.out",
        stderr = "logs/{genera}/2_Taxonomic_Assignment/2_phylogenomics/PhyloPhlAn/{assembler}/{sample}/phylo.err"
    shell: 
        """
        module unload miniconda
        source activate /vast/palmer/pi/turner/flg9/conda_envs/PhyloPhlAn

        phylophlan \
        --input {input.dir} \
        --database {params.db} \
        --diversity {params.diversity} \
        --config_file {params.configfile} \
        --genome_extension {params.ext} \
        --output {params.outdir} \
        --nproc {params.cpus} \
        --verbose \
        1>> {log.stdout} 2>> {log.stderr}
        """