rule filter_individual_assemblies:
    """
    Filter contigs that are a minimum of 1.5kb in length.
    This is the standard length for binning, but it can be iterative
    and we can change this length if needed.
    """
    input:
        contigs = "results/{genera}/testing/3_dedup_contigs/{assembler}/individual_metagenome_assembly/{sample}/{sample}_DEDUP95.fasta"
    output:
        "results/{genera}/testing/5_evaluate_assemblies/{assembler}/filter_individual_assemblies/{sample}/assembly_DEDUP95_m1500.fasta"
    params:
        len = 1500,
        threads = 4,
        outdir = "results/{genera}/testing/5_evaluate_assemblies/{assembler}/filter_individual_assemblies/{sample}"
    log:
        stdout = "logs/{genera}/testing/5_evaluate_assemblies/{assembler}/filter_individual_assemblies/{sample}/filter_assemblies.out",
        stderr = "logs/{genera}/testing/5_evaluate_assemblies/{assembler}/filter_individual_assemblies/{sample}/filter_assemblies.err"
    shell:
        """
        module unload miniconda
        module load SeqKit/2.8.1

        # 1. Create outdir
        mkdir -p {params.outdir}

        # 2. Filter out contigs <1.5kb
        set -x
        echo "Running seqkit..." 1>> {log.stdout} 2>> {log.stderr}

        seqkit seq \
        {input.contigs} --min-len {params.len} --threads {params.threads} -o {output} 
        """

rule evaluate_individual_assemblies:
    input:
        expand(
            "results/{genera}/testing/5_evaluate_assemblies/{assembler}/filter_individual_assemblies/{sample}/assembly_DEDUP95_m1500.fasta",
            sample=SAMPLES,
            genera=config["genera"],
            assembler=config["assembler"]
        )
    output:
        "results/{genera}/testing/5_evaluate_assemblies/QUAST/report.html"
    params:
        outdir = "results/{genera}/testing/5_evaluate_assemblies/QUAST",
        labels = lambda wildcards, input: ",".join([
            f"{path.split('/')[-4]}.{path.split('/')[-2]}.individual.m1500"
            for path in input
        ]),
        threads = 4
    log:
        stdout = "logs/{genera}/5_evaluate_assemblies/QUAST/assembly_eval.out",
        stderr = "logs/{genera}/5_evaluate_assemblies/QUAST/assembly_eval.err"
    shell:
        """
        module unload miniconda
        module load BBMap/38.90-GCCcore-10.2.0
        source activate /home/flg9/.conda/envs/quast

        mkdir -p {params.outdir}

        echo "Running stats.sh..." 1>> {log.stdout} 2>> {log.stderr}

        statswrapper.sh in={input}

        echo "Running metaquast.py..." 1>> {log.stdout} 2>> {log.stderr}

        metaquast.py \
        -t {params.threads} \
        --labels {params.labels} \
        --output-dir {params.outdir} {input} \
        1>> {log.stdout} 2>> {log.stderr}
        """