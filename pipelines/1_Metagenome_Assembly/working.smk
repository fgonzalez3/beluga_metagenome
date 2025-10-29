



# Compare multi-sample binning to single sample binning with SemiBin

# Multi-Sample Binning

rule semibin_generate_concatenated_db:
    """
    Generate concatenated FASTA file necessary for SembiBin's multi-sample binning pipeline
    """
    input:
        fa = lambda wildcards: expand(
        "results/{genera}/testing/3_dedup_contigs/{assembler}/individual_metagenome_assembly/{sample}/{sample}_DEDUP95.fasta",
        genera=config["genera"],
        sample=SAMPLES,
        assembler=[wildcards.assembler])
    output:
        "results/{genera}/testing/6_binning/semibin2/{assembler}_individual_assembly/generate_concatenated_db/concatenated.fa"
    params:
        outdir = "results/{genera}/testing/6_binning/semibin2/{assembler}_individual_assembly/generate_concatenated_db",
        compress = "none",
        threads = 1
    log:
        stdout = "logs/{genera}/testing/6_binning/semibin2/{assembler}_individual_assembly/generate_concatenated_db/concatenate_fa.out",
        stderr = "logs/{genera}/testing/6_binning/semibin2/{assembler}_individual_assembly/generate_concatenated_db/concatenate_fa.err"
    shell:
        """
        module unload miniconda
        source activate /vast/palmer/pi/turner/flg9/conda_envs/semibin

        SemiBin2 concatenate_fasta \
        --input-fasta {input.fa} \
        --output {params.outdir} --compression={params.compress} \
        1>> {log.stdout} 2>> {log.stderr}
        """

rule semibin_align_to_concatenated_db:
    """
    Align reads from each sample to our concatenated FASTA db, necessary for SemiBin pipeline
    """
    input:
        db = "results/{genera}/testing/6_binning/semibin2/{assembler}_individual_assembly/generate_concatenated_db/concatenated.fa",
        r1 = "results/{genera}/1_pre_processing/dedup_reads/{sample}/{sample}_host_removed_dedup_R1.fastq",
        r2 = "results/{genera}/1_pre_processing/dedup_reads/{sample}/{sample}_host_removed_dedup_R2.fastq"
    output:
        "results/{genera}/testing/6_binning/semibin2/{assembler}_individual_assembly/align_to_concatenated_db/{sample}/{sample}_indexed_contig.1.bt2",
        "results/{genera}/testing/6_binning/semibin2/{assembler}_individual_assembly/align_to_concatenated_db/{sample}/{sample}_indexed_contig.2.bt2",
        "results/{genera}/testing/6_binning/semibin2/{assembler}_individual_assembly/align_to_concatenated_db/{sample}/{sample}_indexed_contig.3.bt2",
        "results/{genera}/testing/6_binning/semibin2/{assembler}_individual_assembly/align_to_concatenated_db/{sample}/{sample}_indexed_contig.4.bt2",
        "results/{genera}/testing/6_binning/semibin2/{assembler}_individual_assembly/align_to_concatenated_db/{sample}/{sample}_indexed_contig.rev.1.bt2",
        "results/{genera}/testing/6_binning/semibin2/{assembler}_individual_assembly/align_to_concatenated_db/{sample}/{sample}_indexed_contig.rev.2.bt2",
        "results/{genera}/testing/6_binning/semibin2/{assembler}_individual_assembly/align_to_concatenated_db/{sample}_aligned_sorted.bam"
    params:
        genera=config["genera"],
        outdir = "results/{genera}/testing/6_binning/semibin2/{assembler}_individual_assembly/align_to_concatenated_db/{sample}"
    log:
        stdout = "logs/{genera}/testing/6_binning/semibin2/{assembler}_individual_assembly/align_to_concatenated_db/{sample}_aln.out",
        stderr = "logs/{genera}/testing/6_binning/semibin2/{assembler}_individual_assembly/align_to_concatenated_db/{sample}_aln.err"
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

checkpoint semibin_easy_multi_bin:
    """
    Run Semibin2 on easy multi binning mode
    """
    input:
        db = "results/{genera}/testing/6_binning/semibin2/{assembler}_individual_assembly/generate_concatenated_db/concatenated.fa",
        bams = expand("results/{genera}/testing/6_binning/semibin2/{assembler}_individual_assembly/align_to_concatenated_db/{sample}_aligned_sorted.bam", 
        genera=config["genera"], 
        sample=SAMPLES,
        assembler=config["assembler"])
    output:
        outdir = directory("results/{genera}/testing/6_binning/semibin2/{assembler}_individual_assembly/multi_binning/output_bins"),
        check = "results/{genera}/testing/6_binning/semibin2/{assembler}_individual_assembly/multi_binning/check.tsv"
    params:
        base = "results/{genera}/testing/6_binning/semibin2/{assembler}_individual_assembly/multi_binning/",
        outdir = "results/{genera}/testing/6_binning/semibin2/{assembler}_individual_assembly/multi_binning/output_bins",
        compress = "none",
        minlen = 1500,
        threads = 1
    log:
        stdout = "logs/{genera}/testing/6_binning/semibin2/{assembler}_individual_assembly/multi_binning/multi_binning.out",
        stderr = "logs/{genera}/testing/6_binning/semibin2/{assembler}_individual_assembly/multi_binning/multi_binning.err"
    shell:
        """
        module unload miniconda 
        source activate /vast/palmer/pi/turner/flg9/conda_envs/semibin

        SemiBin2 multi_easy_bin \
        -i {input.db} \
        -b {input.bams} \
        --self-supervised \
        --min-len {params.minlen} \
        --threads {params.threads} \
        -o {params.base} \
        --compression {params.compress} \
        1>> {log.stdout} 2>> {log.stderr}

        ls {params.outdir}/*.fa > {output.check}
        """

def get_semibin2_multibin(wc):
    ckpt = checkpoints.semibin_easy_multi_bin.get(genera=wc.genera)
    import os, glob
    bins_dir = ckpt.output.outdir
    return sorted(glob.glob(os.path.join(bins_dir, "*.fa")))