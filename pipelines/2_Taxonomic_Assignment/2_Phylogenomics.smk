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
        cpus = 4,
        diversity = "high",
        db = "phylophlan", # this db uses 400 marker genes # amphora db is the other default option with 136 marker genes
        configfile = "/vast/palmer/pi/turner/flg9/conda_envs/PhyloPhlAn/configs/supermatrix_nt.cfg",
        outdir = "results/{genera}/2_Taxonomic_Assignment/2_phylogenomics/PhyloPhlAn/{assembler}/{sample}"
    log:
        stdout = "logs/{genera}/2_Taxonomic_Assignment/2_phylogenomics/PhyloPhlAn/{assembler}/{sample}/phylo.out",
        stderr = "logs/{genera}/2_Taxonomic_Assignment/2_phylogenomics/PhyloPhlAn/{assembler}/{sample}/phylo.err"
    shell: 
        """
        module unload miniconda
        source activate /vast/palmer/pi/turner/flg9/conda_envs/PhyloPhlAn

        phylophlan \
        -i {input.dir} \
        -d {params.db} \
        --diversity {params.diversity} \
        -f {params.configfile} \
        -output_folder {params.outdir} \
        --nproc {params.cpus}
        1>> {log.stdout} 2>> {log.stderr}
        """