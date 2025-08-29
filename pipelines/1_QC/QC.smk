import os
import glob
import pandas as pd

configfile: "config/assembly.yaml"

samples_df = pd.read_csv("tsv/beluga_raw_reads.tsv", sep="\t")
SAMPLES = samples_df["sample_id"].tolist()
READS = {row.sample_id: {"r1": row.r1, "r2": row.r2} for row in samples_df.itertuples()}

rule all:
    input:
        # 1. QC pipeline results # Will change 1_assembly outdir to 1_QC for later versions of this pipeline
        expand("results/{genera}/1_assembly/quality_control/{sample}/{sample}_val_1.fq.gz", sample=SAMPLES, genera=config["genera"]),
        expand("results/{genera}/1_assembly/quality_control/{sample}/{sample}_val_2.fq.gz", sample=SAMPLES, genera=config["genera"]),
        expand("results/{genera}/1_assembly/quality_control/{sample}/{sample}_val_1_fastqc.html", sample=SAMPLES, genera=config["genera"]),
        expand("results/{genera}/1_assembly/quality_control/{sample}/{sample}_val_2_fastqc.html", sample=SAMPLES, genera=config["genera"]),
        expand("results/{genera}/1_assembly/aggregate_qc_data/multiqc_report.html", genera=config["genera"]),
        expand("results/{genera}/1_assembly/mask_beluga/beluga_genome_sequence_masked.fa", genera=config["genera"]),
        expand("results/{genera}/1_assembly/map_beluga/{sample}/mapped.fq", sample=SAMPLES, genera=config["genera"]),
        expand("results/{genera}/1_assembly/map_beluga/{sample}/unmapped.fq", sample=SAMPLES, genera=config["genera"]),
        expand("results/{genera}/1_assembly/map_beluga/{sample}/unmapped_R1.fq", sample=SAMPLES, genera=config["genera"]),
        expand("results/{genera}/1_assembly/map_beluga/{sample}/unmapped_R2.fq", sample=SAMPLES, genera=config["genera"]),
        expand("results/{genera}/1_assembly/mask_human/human_genome_sequence_masked.fa", genera=config["genera"]),
        expand("results/{genera}/1_assembly/map_human/{sample}/mapped.fq", sample=SAMPLES, genera=config["genera"]),
        expand("results/{genera}/1_assembly/map_human/{sample}/unmapped.fq", sample=SAMPLES, genera=config["genera"]),
        expand("results/{genera}/1_assembly/map_human/{sample}/unmapped_R1.fq", sample=SAMPLES, genera=config["genera"]),
        expand("results/{genera}/1_assembly/map_human/{sample}/unmapped_R2.fq", sample=SAMPLES, genera=config["genera"]),
        expand("results/{genera}/1_assembly/dedup_reads/{sample}/{sample}_host_removed_dedup_R1.fastq", sample=SAMPLES, genera=config["genera"]),
        expand("results/{genera}/1_assembly/dedup_reads/{sample}/{sample}_host_removed_dedup_R2.fastq", sample=SAMPLES, genera=config["genera"]),
        expand("results/{genera}/1_assembly/co_assembly_reads/all_samples_R1.fq", genera=config["genera"]),
        expand("results/{genera}/1_assembly/co_assembly_reads/all_samples_R2.fq", genera=config["genera"]),
        expand("results/{genera}/1_assembly/co_assembly_reads/all_samples_R1_norm.fq", genera=config["genera"]),
        expand("results/{genera}/1_assembly/co_assembly_reads/all_samples_R2_norm.fq", genera=config["genera"]),

        # 2. Assembly pipeline results 
        expand("results/{genera}/1_assembly/metagenome_assembly/{sample}/contigs.fasta", sample=SAMPLES, genera=config["genera"]),
        expand("results/{genera}/1_assembly/metagenome_assembly_megahit/{sample}/final.contigs.fa", sample=SAMPLES, genera=config["genera"]),
        expand("results/{genera}/1_assembly/dedup_contigs_spades/{sample}/{sample}_DEDUP95.fasta", sample=SAMPLES, genera=config["genera"]),
        expand("results/{genera}/1_assembly/dedup_contigs_megahit/{sample}/{sample}_DEDUP95.fasta", sample=SAMPLES, genera=config["genera"]),
        expand("results/{genera}/1_assembly/contig_index_spades/{sample}/{sample}_indexed_contig.1.bt2", sample=SAMPLES, genera=config["genera"]),
        expand("results/{genera}/1_assembly/contig_index_spades/{sample}/{sample}_indexed_contig.2.bt2", sample=SAMPLES, genera=config["genera"]),
        expand("results/{genera}/1_assembly/contig_index_spades/{sample}/{sample}_indexed_contig.3.bt2", sample=SAMPLES, genera=config["genera"]),
        expand("results/{genera}/1_assembly/contig_index_spades/{sample}/{sample}_indexed_contig.4.bt2", sample=SAMPLES, genera=config["genera"]),
        expand("results/{genera}/1_assembly/contig_index_spades/{sample}/{sample}_indexed_contig.rev.1.bt2", sample=SAMPLES, genera=config["genera"]),
        expand("results/{genera}/1_assembly/contig_index_spades/{sample}/{sample}_indexed_contig.rev.2.bt2", sample=SAMPLES, genera=config["genera"]),
        expand("results/{genera}/1_assembly/contig_read_alignment_spades/{sample}_aligned_sorted.bam", sample=SAMPLES, genera=config["genera"]),
        expand("results/{genera}/1_assembly/contig_index_megahit/{sample}/{sample}_indexed_contig.1.bt2", sample=SAMPLES, genera=config["genera"]),
        expand("results/{genera}/1_assembly/contig_index_megahit/{sample}/{sample}_indexed_contig.2.bt2", sample=SAMPLES, genera=config["genera"]),
        expand("results/{genera}/1_assembly/contig_index_megahit/{sample}/{sample}_indexed_contig.3.bt2", sample=SAMPLES, genera=config["genera"]),
        expand("results/{genera}/1_assembly/contig_index_megahit/{sample}/{sample}_indexed_contig.4.bt2", sample=SAMPLES, genera=config["genera"]),
        expand("results/{genera}/1_assembly/contig_index_megahit/{sample}/{sample}_indexed_contig.rev.1.bt2", sample=SAMPLES, genera=config["genera"]),
        expand("results/{genera}/1_assembly/contig_index_megahit/{sample}/{sample}_indexed_contig.rev.2.bt2", sample=SAMPLES, genera=config["genera"]),
        expand("results/{genera}/1_assembly/contig_read_alignment_megahit/{sample}_aligned_sorted.bam", sample=SAMPLES, genera=config["genera"]),
        expand("results/{genera}/1_assembly/assembly_eval/{sample}/metaspades_assembly_DEDUP95_m1500.fasta", sample=SAMPLES, genera=config["genera"]),
        expand("results/{genera}/1_assembly/assembly_eval/{sample}/assembly_stats.csv", sample=SAMPLES, genera=config["genera"]),
        expand("results/{genera}/1_assembly/assembly_eval/{sample}/report.html", sample=SAMPLES, genera=config["genera"]),

        # 3. Binning pipeline results 
        expand("results/{genera}/2_binning/concoct/{sample}/CONCOCT.*.fa", genera=config["genera"], sample=SAMPLES),
        expand("results/{genera}/2_binning/concoct/{sample}/concoct_output/clustering_merged.csv", genera=config["genera"], sample=SAMPLES),
        expand("results/{genera}/2_binning/metabat/{sample}/METABAT.txt", genera=config["genera"], sample=SAMPLES),
        expand("results/{genera}/2_binning/metabat/{sample}/METABAT.*.fa", genera=config["genera"], sample=SAMPLES),
        expand("results/{genera}/2_binning/maxbin/maxbin.txt", genera=config["genera"], sample=SAMPLES),
        expand("results/{genera}/2_binning/maxbin/{sample}/MAXBIN.*.fa", genera=config["genera"], sample=SAMPLES),
        expand("results/{genera}/2_binning/maxbin/{sample}/MAXBIN.summary", genera=config["genera"], sample=SAMPLES),
        expand("results/{genera}/2_binning/semibin2/data_split.csv", genera=config["genera"], sample=SAMPLES),
        expand("results/{genera}/2_binning/semibin2/data.csv", genera=config["genera"], sample=SAMPLES),
        expand("results/{genera}/2_binning/semibin2/model.pt", genera=config["genera"], sample=SAMPLES),
        expand("results/{genera}/2_binning/semibin2/{sample}/bin.*.fa", genera=config["genera"], sample=SAMPLES)
       
rule quality_control:
    """
    Trim adapters and low-quality reads
    """
    input:
        fwd = lambda wildcards: READS[wildcards.sample]["r1"],
        rev = lambda wildcards: READS[wildcards.sample]["r2"]
    output:
        "results/{genera}/1_assembly/quality_control/{sample}/{sample}_val_1.fq.gz",
        "results/{genera}/1_assembly/quality_control/{sample}/{sample}_val_2.fq.gz",
        "results/{genera}/1_assembly/quality_control/{sample}/{sample}_val_1_fastqc.html",
        "results/{genera}/1_assembly/quality_control/{sample}/{sample}_val_2_fastqc.html"
    params:
        genera = config["genera"],
        outdir = "results/{genera}/1_assembly/quality_control/{sample}"
    log:
        stdout = "logs/{genera}/1_assembly/quality_control/{sample}/qc.out",
        stderr = "logs/{genera}/1_assembly/quality_control/{sample}/qc.err"
    shell:
        """
        module unload miniconda
        module load Trim_Galore/0.6.7-GCCcore-10.2.0

        trim_galore \
        --paired --fastqc --length 125 --cores 8 -q 25 {input.fwd} {input.rev} \
        --output_dir {params.outdir} --basename {wildcards.sample} \
        1> {log.stdout} 2> {log.stderr}
        """

rule aggregate_qc_data:
    """
    Aggregate QC data to create a single report across our many samples
    """
    input:
        expand("results/{genera}/1_assembly/quality_control/{sample}/{sample}_val_1_fastqc.html", genera=config["genera"], sample=SAMPLES),
        expand("results/{genera}/1_assembly/quality_control/{sample}/{sample}_val_2_fastqc.html", genera=config["genera"], sample=SAMPLES)
    output:
        "results/{genera}/1_assembly/aggregate_qc_data/multiqc_report.html"
    params:
        input_dirs = lambda wildcards: f"results/{wildcards.genera}/1_assembly/quality_control",
        outdir = "results/{genera}/1_assembly/aggregate_qc_data"
    log:
        stdout = "logs/{genera}/1_assembly/aggregate_qc_data/multiqc.out",
        stderr = "logs/{genera}/1_assembly/aggregate_qc_data/multiqc.err"
    shell:
        """
        module unload miniconda 
        module load MultiQC/1.10.1-foss-2020b-Python-3.8.6

        multiqc \
        {params.input_dirs} -o {params.outdir} \
        1> {log.stdout} 2> {log.stderr}
        """

rule mask_beluga_host_genome:
    """
    Mask low-complexity or microbial contaminant regions in a beluga reference sequence, reducing false positives during host read filtering
    """
    input:
        beluga_refseq = config["beluga_refseq"]
    output:
        beluga_ref_masked = "results/{genera}/1_assembly/mask_beluga/beluga_genome_sequence_masked.fa"
    params:
        genera=config["genera"],
        threads = 4
    log:
        stdout = "logs/{genera}/1_assembly/mask_beluga/bbmask.out",
        stderr = "logs/{genera}/1_assembly/mask_beluga/bbmask.err"
    shell:
        """
        module unload miniconda
        module load BBMap/38.90-GCCcore-10.2.0

        bbmask.sh \
        in={input.beluga_refseq} out={output.beluga_ref_masked} threads={params.threads} \
        1> {log.stdout} 2> {log.stderr}
        """

rule align_to_beluga_host_genome:
    """
    Filter reads against masked beluga reference genome
    """
    input:
        beluga_ref_masked = "results/{genera}/1_assembly/mask_beluga/beluga_genome_sequence_masked.fa",
        r1 = "results/{genera}/1_assembly/quality_control/{sample}/{sample}_val_1.fq.gz",
        r2 = "results/{genera}/1_assembly/quality_control/{sample}/{sample}_val_2.fq.gz"
    output:
        mapped = "results/{genera}/1_assembly/map_beluga/{sample}/mapped.fq",
        unmapped = "results/{genera}/1_assembly/map_beluga/{sample}/unmapped.fq"
    params:
        index_path = "results/{genera}/1_assembly/map_beluga/index/",
        genera=config["genera"],
        id=0.95,
        indel=3,
        ratio=0.16,
        bandwidth=12,
        match="fast",
        hits=2,
        mem="25g",
        quickmatch="t",
        fast="t"
    log:
        stdout = "logs/{genera}/1_assembly/map_beluga/{sample}/bbmap_beluga.out",
        stderr = "logs/{genera}/1_assembly/map_beluga/{sample}/bbmap_beluga.err"
    shell:
        """
        module unload miniconda
        module load BBMap/38.90-GCCcore-10.2.0

        bbmap.sh \
        minid={params.id} maxindel={params.indel} bwr={params.ratio} bw={params.bandwidth} \
        quickmatch={params.quickmatch} fast={params.fast} minhits={params.hits} \
        ref={input.beluga_ref_masked} in1={input.r1} in2={input.r2} outm={output.mapped} \
        outu={output.unmapped} path={params.index_path} -Xmx{params.mem} \
        1> {log.stdout} 2> {log.stderr}
        """

rule split_unmapped_reads_beluga:
    """
    Split interleaved reads that did not map to masked beluga reference into separate fastqs
    """
    input:
        unmapped = "results/{genera}/1_assembly/map_beluga/{sample}/unmapped.fq"
    output:
        r1 = "results/{genera}/1_assembly/map_beluga/{sample}/unmapped_R1.fq",
        r2 = "results/{genera}/1_assembly/map_beluga/{sample}/unmapped_R2.fq"
    params:
        genera=config["genera"]
    log:
        stdout = "logs/{genera}/1_assembly/map_beluga/{sample}/reformat_beluga.out",
        stderr = "logs/{genera}/1_assembly/map_beluga/{sample}/reformat_beluga.err"
    shell:
        """
        module unload miniconda
        module load BBMap/38.90-GCCcore-10.2.0

        reformat.sh \
        in={input.unmapped} out1={output.r1} out2={output.r2} \
        1> {log.stdout} 2> {log.stderr}
        """

rule mask_human_host_genome:
    """
    Mask low-complexity or microbial contaminant regions in a human reference sequence, reducing false positives during host read filtering
    """
    input:
        human_refseq = config["human_refseq"]
    output:
        human_ref_masked = "results/{genera}/1_assembly/mask_human/human_genome_sequence_masked.fa"
    params:
        genera=config["genera"],
        threads=4
    log:
        stdout = "logs/{genera}/1_assembly/mask_human/bbmask.out",
        stderr = "logs/{genera}/1_assembly/mask_human/bbmask.err"
    shell:
        """
        module unload miniconda
        module load BBMap/38.90-GCCcore-10.2.0

        bbmask.sh \
        in={input.human_refseq} out={output.human_ref_masked} threads={params.threads} \
        1> {log.stdout} 2> {log.stderr}
        """

rule align_to_human_host_genome:
    """
    Filter reads against masked human reference genome
    """
    input:
        human_ref_masked = "results/{genera}/1_assembly/mask_human/human_genome_sequence_masked.fa",
        r1 = "results/{genera}/1_assembly/map_beluga/{sample}/unmapped_R1.fq",
        r2 = "results/{genera}/1_assembly/map_beluga/{sample}/unmapped_R2.fq"
    output:
        mapped = "results/{genera}/1_assembly/map_human/{sample}/mapped.fq",
        unmapped = "results/{genera}/1_assembly/map_human/{sample}/unmapped.fq"
    params:
        index_path = "results/{genera}/1_assembly/map_human/index/",
        genera=config["genera"],
        id=0.95,
        indel=3,
        ratio=0.16,
        bandwidth=12,
        match="fast",
        hits=2,
        mem="25g",
        quickmatch="t",
        fast="t"
    log:
        stdout = "logs/{genera}/1_assembly/map_human/{sample}/bbmap_human.out",
        stderr = "logs/{genera}/1_assembly/map_human/{sample}/bbmap_human.err"
    shell:
        """
        module unload miniconda
        module load BBMap/38.90-GCCcore-10.2.0

        bbmap.sh \
        minid={params.id} maxindel={params.indel} bwr={params.ratio} bw={params.bandwidth} \
        quickmatch={params.quickmatch} fast={params.fast} minhits={params.hits} \
        ref={input.human_ref_masked} in1={input.r1} in2={input.r2} outm={output.mapped} \
        outu={output.unmapped} path={params.index_path} -Xmx{params.mem} \
        1> {log.stdout} 2> {log.stderr}
        """

rule split_unmapped_reads_human:
    """
    Split interleaved reads that did not map to masked human reference into separate fastqs
    """
    input:
        unmapped = "results/{genera}/1_assembly/map_human/{sample}/unmapped.fq"
    output:
        r1 = "results/{genera}/1_assembly/map_human/{sample}/unmapped_R1.fq",
        r2 = "results/{genera}/1_assembly/map_human/{sample}/unmapped_R2.fq"
    params:
        genera=config["genera"]
    log:
        stdout = "logs/{genera}/1_assembly/map_human/{sample}/reformat_human.out",
        stderr = "logs/{genera}/1_assembly/map_human/{sample}/reformat_human.err"
    shell:
        """
        module unload miniconda
        module load BBMap/38.90-GCCcore-10.2.0

        reformat.sh \
        in={input.unmapped} out1={output.r1} out2={output.r2} \
        1> {log.stdout} 2> {log.stderr}
        """

rule deduplicate_reads:
    """
    Run CD-HIT to deduplicate reads that may have been introduced during NexteraXT PCR amplification
    input:
        r1 = "results/{genera}/1_assembly/map_human/{sample}/unmapped_R1.fq",
        r2 = "results/{genera}/1_assembly/map_human/{sample}/unmapped_R2.fq"
    output:
        r1_dedup = "results/{genera}/1_assembly/dedup_reads/{sample}/{sample}_host_removed_dedup_R1.fastq",
        r2_dedup = "results/{genera}/1_assembly/dedup_reads/{sample}/{sample}_host_removed_dedup_R2.fastq"
    params:
        genera=config["genera"],
        outdir = "results/{genera}/1_assembly/dedup_reads/{sample}"
    log:
        stdout = "logs/{genera}/1_assembly/dedup_reads/{sample}/dedup_reads.out",
        stderr = "logs/{genera}/1_assembly/dedup_reads/{sample}/dedup_reads.err"
    shell:
        """
        module unload miniconda
        source activate /home/flg9/.conda/envs/cd-hit_ycrc
        mkdir -p {params.outdir}

        cd-hit-dup \
        -i {input.r1} -i2 {input.r2} -o {output.r1_dedup} -o2 {output.r2_dedup} \
        1> {log.stdout} 2> {log.stderr}
        """

rule normalization:
    """
    Concatenate and normalize read pairs prior to co-assembly
    """
    input:
        r1 = "results/{genera}/1_assembly/map_human/{sample}/unmapped_R1.fq",
        r2 = "results/{genera}/1_assembly/map_human/{sample}/unmapped_R2.fq"
    output:
        r1 = "results/{genera}/1_assembly/co_assembly_reads/all_samples_R1.fq",
        r2 = "results/{genera}/1_assembly/co_assembly_reads/all_samples_R2.fq",
        norm1 = "results/{genera}/1_assembly/co_assembly_reads/all_samples_R1_norm.fq",
        norm2 = "results/{genera}/1_assembly/co_assembly_reads/all_samples_R2_norm.fq"
    params:
        genera=config["genera"],
        target=70,
        mindepth=2,
        threads=4,
        prefilter"t"
    log:
        stdout = "logs/{genera}/1_assembly/normalization/{sample}/dedup_reads.out",
        stderr = "logs/{genera}/1_assembly/normalization/{sample}/dedup_reads.err"
    shell:
        """
        module unload miniconda 
        module load BBMap/38.90-GCCcore-10.2.0

        cat {input.fwd} >> {output.r1}
        cat {input.rev} >> {output.r2}

        bbnorm.sh in={output.r1} in2={output.r2} \
        out={output.norm1} out2={output.norm2} \
        target={params.target} mindepth={params.mindepth} \
        threads={params.threads} prefilter={params.prefilter} \
        1> {log.stdout} 2> {log.stderr}
        """

include: "pipelines/assembly.smk"
include: "pipelines/binning.smk"
#include: pipelines/taxonomic_classification.smk
#include: pipelines/virome_characterization.smk