import os
import glob
import pandas as pd

configfile: "config/assembly.yaml"

samples_df = pd.read_csv("tsv/beluga_raw_reads.tsv", sep="\t")
SAMPLES = samples_df["sample_id"].tolist()
READS = {row.sample_id: {"r1": row.r1, "r2": row.r2} for row in samples_df.itertuples()}

rule all:
    input:
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
        expand("results/{genera}/1_assembly/metagenome_assembly/{sample}/contigs.fasta", sample=SAMPLES, genera=config["genera"]),
        expand("results/{genera}/1_assembly/dedup_contigs/{sample}/{sample}_DEDUP95.fasta", sample=SAMPLES, genera=config["genera"]),
        expand("results/{genera}/1_assembly/contig_index/{sample}/{sample}_indexed_contig.1.bt2", sample=SAMPLES, genera=config["genera"]),
        expand("results/{genera}/1_assembly/contig_index/{sample}/{sample}_indexed_contig.2.bt2", sample=SAMPLES, genera=config["genera"]),
        expand("results/{genera}/1_assembly/contig_index/{sample}/{sample}_indexed_contig.3.bt2", sample=SAMPLES, genera=config["genera"]),
        expand("results/{genera}/1_assembly/contig_index/{sample}/{sample}_indexed_contig.4.bt2", sample=SAMPLES, genera=config["genera"]),
        expand("results/{genera}/1_assembly/contig_index/{sample}/{sample}_indexed_contig.rev.1.bt2", sample=SAMPLES, genera=config["genera"]),
        expand("results/{genera}/1_assembly/contig_index/{sample}/{sample}_indexed_contig.rev.2.bt2", sample=SAMPLES, genera=config["genera"]),
        expand("results/{genera}/1_assembly/contig_read_alignment/{sample}_aligned_sorted.bam", sample=SAMPLES, genera=config["genera"]),
        expand("results/{genera}/1_assembly/contig_aln_stats/{sample}_mappingstats.txt", sample=SAMPLES, genera=config["genera"]),
        expand("results/{genera}/1_assembly/assembly_eval/{sample}/metaspades_assembly_DEDUP95_m1500.fasta", sample=SAMPLES, genera=config["genera"]),
        expand("results/{genera}/1_assembly/assembly_eval/{sample}/assembly_stats.csv", sample=SAMPLES, genera=config["genera"]),
        expand("results/{genera}/1_assembly/assembly_eval/{sample}/report.html", sample=SAMPLES, genera=config["genera"])

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
    """
    input:
        #r1 = "results/{genera}/1_assembly/bowtie_align/{sample}/{sample}_host_removed_R1.fastq",
        #r2 = "results/{genera}/1_assembly/bowtie_align/{sample}/{sample}_host_removed_R2.fastq"
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

rule metagenome_assembly:
    """
    Assemble contigs
    """
    input:
        fwd = "results/{genera}/1_assembly/dedup_reads/{sample}/{sample}_host_removed_dedup_R1.fastq",
        rev = "results/{genera}/1_assembly/dedup_reads/{sample}/{sample}_host_removed_dedup_R2.fastq"
    output:
        "results/{genera}/1_assembly/metagenome_assembly/{sample}/contigs.fasta"
    params:
        genera=config["genera"],
        outdir = "results/{genera}/1_assembly/metagenome_assembly/{sample}"
    log:
        stdout = "logs/{genera}/1_assembly/metagenome_assembly/{sample}/assembly.out",
        stderr = "logs/{genera}/1_assembly/metagenome_assembly/{sample}/assembly.err"
    resources:
        mem_mb=200000,
        threads=4
    shell:
        """
        module unload miniconda
        module load SPAdes/3.15.5-GCC-12.2.0

        spades.py \
        --meta --threads {resources.threads} -1 {input.fwd} -2 {input.rev} -o {params.outdir} \
        1> {log.stdout} 2> {log.stderr}
        """

rule deduplicate_contigs:
    """
    Run CD-HIT to deduplicate contigs that may have been introduced during metagenome assembly
    """
    input:
        "results/{genera}/1_assembly/metagenome_assembly/{sample}/contigs.fasta"
    output:
        dedups = "results/{genera}/1_assembly/dedup_contigs/{sample}/{sample}_DEDUP95.fasta"
    params:
        genera=config["genera"],
        outdir = "results/{genera}/1_assembly/dedup_contigs/{sample}"
    log:
        stdout = "logs/{genera}/1_assembly/dedup_contigs/{sample}/dedup_contigs.out",
        stderr = "logs/{genera}/1_assembly/dedup_contigs/{sample}/dedup_contigs.err"
    shell:
        """
        module unload miniconda
        module load CD-HIT/4.8.1-GCC-10.2.0
        mkdir -p {params.outdir}

        cd-hit-est \
        -i {input} -o {output.dedups} -T 4 -c 0.95 \
        1> {log.stdout} 2> {log.stderr}
        """

rule align_reads_to_contigs:
    """
    Align raw reads back to assembled contigs to assess coverage, assembly quality, and obtain BAM files for binning
    """
    input:
        contigs = "results/{genera}/1_assembly/dedup_contigs/{sample}/{sample}_DEDUP95.fasta",
        r1 = "results/{genera}/1_assembly/dedup_reads/{sample}/{sample}_host_removed_dedup_R1.fastq",
        r2 = "results/{genera}/1_assembly/dedup_reads/{sample}/{sample}_host_removed_dedup_R2.fastq"
    output:
        "results/{genera}/1_assembly/contig_index/{sample}/{sample}_indexed_contig.1.bt2",
        "results/{genera}/1_assembly/contig_index/{sample}/{sample}_indexed_contig.2.bt2",
        "results/{genera}/1_assembly/contig_index/{sample}/{sample}_indexed_contig.3.bt2",
        "results/{genera}/1_assembly/contig_index/{sample}/{sample}_indexed_contig.4.bt2",
        "results/{genera}/1_assembly/contig_index/{sample}/{sample}_indexed_contig.rev.1.bt2",
        "results/{genera}/1_assembly/contig_index/{sample}/{sample}_indexed_contig.rev.2.bt2",
        "results/{genera}/1_assembly/contig_read_alignment/{sample}_aligned_sorted.bam"
    params:
        genera=config["genera"],
        outdir = "results/{genera}/1_assembly/contig_index/{sample}"
    log:
        stdout = "logs/{genera}/1_assembly/contig_read_alignment/{sample}_aln.out",
        stderr = "logs/{genera}/1_assembly/contig_read_alignment/{sample}_aln.err"
    shell:
        """
        module unload miniconda 
        module load Bowtie2/2.5.1-GCC-12.2.0
        module load SAMtools/1.21-GCC-12.2.0

        # 1. Build contig index
        bowtie2-build \
        -f {input.contigs} {params.outdir}/{wildcards.sample}_indexed_contig \
        1>> {log.stdout} 2>> {log.stderr}

        # 2. Align reads back to assembled contigs
        bowtie2 \
        -x {params.outdir}/{wildcards.sample}_indexed_contig -1 {input.r1} -2 {input.r2} | samtools view -b -F 4 -F 2048 | samtools sort -o {output[6]} \
        1>> {log.stdout} 2>> {log.stderr}
        """

rule contig_quality:
    """
    Assess contig alignments
    """
    input:
        "results/{genera}/1_assembly/contig_read_alignment/{sample}_aligned_sorted.bam"
    output:
        "results/{genera}/1_assembly/contig_aln_stats/{sample}_mappingstats.txt"
    params:
        genera=config["genera"]
    log:
        stdout = "logs/{genera}/1_assembly/contig_aln_stats/{sample}_contig_qc.out",
        stderr = "logs/{genera}/1_assembly/contig_aln_stats/{sample}_contig_qc.err"
    shell:
        """
        module unload miniconda 
        module load SAMtools/1.21-GCC-12.2.0

        samtools flagstat \
        {input} > {output} \
        1> {log.stdout} 2> {log.stderr}
        """

rule assembly_eval:
    """
    Evaluate assemblies
    """
    input:
        contigs = "results/{genera}/1_assembly/dedup_contigs/{sample}/{sample}_DEDUP95.fasta"
    output:
        filtered_contigs = "results/{genera}/1_assembly/assembly_eval/{sample}/metaspades_assembly_DEDUP95_m1500.fasta",
        whole_assembly_stats = "results/{genera}/1_assembly/assembly_eval/{sample}/assembly_stats.csv",
        metaquast_output = "results/{genera}/1_assembly/assembly_eval/{sample}/report.html"
    params:
        len = 1500,
        threads = 4,
        outdir = "results/{genera}/1_assembly/assembly_eval/{sample}",
        outfile_label = "SPAdes.m1500"
    log:
        stdout = "logs/{genera}/1_assembly/assembly_eval/{sample}/assembly_eval.out",
        stderr = "logs/{genera}/1_assembly/assembly_eval/{sample}/assembly_eval.err"
    shell:
        """
        module unload miniconda
        module load SeqKit/2.8.1
        module load BBMap/38.90-GCCcore-10.2.0
        source activate /home/flg9/.conda/envs/quast

        mkdir -p {params.outdir}

        set -x
        echo "Running seqkit..." 1>> {log.stdout} 2>> {log.stderr}

        #1. Filter out contigs <1.5kb
        seqkit seq \
        {input.contigs} --min-len {params.len} --threads {params.threads} -o {output.filtered_contigs} 

        echo "Running stats.sh..." 1>> {log.stdout} 2>> {log.stderr}

        #2. Generate specific assembly data covering quality of contigs generated
        stats.sh \
        in={output.filtered_contigs} out={output.whole_assembly_stats}

        echo "Running quast.py..." 1>> {log.stdout} 2>> {log.stderr}

        #3. Lastly, generate comprehensive evaluation of assembly data using metaquast
        metaquast.py \
        -t {params.threads} --labels {params.outfile_label} --output-dir {params.outdir} \
        {output.filtered_contigs}

        echo "Tarring results ..." 1>> {log.stdout} 2>> {log.stderr}
        tar \
        -cvzf metaquast_results.tar.gz {params.outdir}
        """

#include: pipelines/binning.smk
#include: pipelines/taxonomic_classification.smk
#include: pipelines/virome_characterization.smk