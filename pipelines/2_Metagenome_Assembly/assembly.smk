import os
import glob
import pandas as pd

rule metagenome_assembly_spades:
    """
    Assemble contigs from individual samples using SPAdes
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
        --meta --threads {resources.threads} \
        -1 {input.fwd} -2 {input.rev} -o {params.outdir} \
        1> {log.stdout} 2> {log.stderr}
        """

rule metagenome_assembly_megahit:
    """
    Assemble contigs from individual samples with Megahit
    """
    input:
        fwd = "results/{genera}/1_assembly/dedup_reads/{sample}/{sample}_host_removed_dedup_R1.fastq",
        rev = "results/{genera}/1_assembly/dedup_reads/{sample}/{sample}_host_removed_dedup_R2.fastq"
    output:
        "results/{genera}/1_assembly/metagenome_assembly_megahit/{sample}/final.contigs.fa"
    params:
        genera=config["genera"],
        preset="meta-large"
        outdir = "results/{genera}/1_assembly/metagenome_assembly_megahit/{sample}"
    log:
        stdout = "logs/{genera}/1_assembly/metagenome_assembly_megahit/{sample}/assembly.out",
        stderr = "logs/{genera}/1_assembly/metagenome_assembly_megahit/{sample}/assembly.err"
    shell:
        """
        module unload miniconda
        source activate XXXXXX

        megahit \
        -1 {input.fwd} -2 {input.rev} -o {params.outdir} \
        --presets {params.preset} \
        1> {log.stdout} 2> {log.stderr}
        """

rule metagenome_coassembly_spades:
    """
    """
    input:
    output:
    params:
    log:
    shell:
        """
        """

rule metagenome_coassembly_megahit:
    """
    """
    input:
    output:
    params:
    log:
    shell:
        """
        """

rule deduplicate_contigs:
    """
    Run CD-HIT to deduplicate contigs that may have been introduced during metagenome assembly
    """
    input:
        SPAdes_assembly = "results/{genera}/1_assembly/metagenome_assembly/{sample}/contigs.fasta",
        megahit_assembly = "results/{genera}/1_assembly/metagenome_assembly_megahit/{sample}/final.contigs.fa"
    output:
        SPAdes_dedups = "results/{genera}/1_assembly/dedup_contigs_spades/{sample}/{sample}_DEDUP95.fasta",
        megahit_dedups = "results/{genera}/1_assembly/dedup_contigs_megahit/{sample}/{sample}_DEDUP95.fasta"
    params:
        genera=config["genera"],
        out1 = "results/{genera}/1_assembly/dedup_contigs_spades/{sample}",
        out2 = "results/{genera}/1_assembly/dedup_contigs_megahit/{sample}"
    log:
        stdout = "logs/{genera}/1_assembly/dedup_contigs/{sample}/dedup_contigs.out",
        stderr = "logs/{genera}/1_assembly/dedup_contigs/{sample}/dedup_contigs.err"
    shell:
        """
        module unload miniconda
        module load CD-HIT/4.8.1-GCC-10.2.0

        mkdir -p {params.out1}
        mkdir -p {params.out2}

        # 1. Deduplicate contigs from SPAdes assembly
        cd-hit-est \
        -i {inpu.SPAdes_assembly} -o {out1.dedups} -T 4 -c 0.95 \
        1>> {log.stdout} 2>> {log.stderr}

        # 2. Deduplicate contigs from Megahit assembly
        cd-hit-est \
        -i {input.megahit_assembly} -o {out2.dedups} -T 4 -c 0.95 \
        1>> {log.stdout} 2>> {log.stderr}
        """

rule align_reads_to_contigs_spades:
    """
    Align raw reads back to assembled contigs to assess coverage, assembly quality, and obtain BAM files for binning from SPAdes assembly
    """
    input:
        contigs = "results/{genera}/1_assembly/dedup_contigs_spades/{sample}/{sample}_DEDUP95.fasta",
        r1 = "results/{genera}/1_assembly/dedup_reads/{sample}/{sample}_host_removed_dedup_R1.fastq",
        r2 = "results/{genera}/1_assembly/dedup_reads/{sample}/{sample}_host_removed_dedup_R2.fastq"
    output:
        "results/{genera}/1_assembly/contig_index_spades/{sample}/{sample}_indexed_contig.1.bt2",
        "results/{genera}/1_assembly/contig_index_spades/{sample}/{sample}_indexed_contig.2.bt2",
        "results/{genera}/1_assembly/contig_index_spades/{sample}/{sample}_indexed_contig.3.bt2",
        "results/{genera}/1_assembly/contig_index_spades/{sample}/{sample}_indexed_contig.4.bt2",
        "results/{genera}/1_assembly/contig_index_spades/{sample}/{sample}_indexed_contig.rev.1.bt2",
        "results/{genera}/1_assembly/contig_index_spades/{sample}/{sample}_indexed_contig.rev.2.bt2",
        "results/{genera}/1_assembly/contig_read_alignment_spades/{sample}_aligned_sorted.bam"
    params:
        genera=config["genera"],
        outdir = "results/{genera}/1_assembly/contig_index_spades/{sample}"
    log:
        stdout = "logs/{genera}/1_assembly/contig_read_alignment_spades/{sample}_aln.out",
        stderr = "logs/{genera}/1_assembly/contig_read_alignment_spades/{sample}_aln.err"
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

rule align_reads_to_contigs_megahit:
    """
    Align raw reads back to assembled contigs to assess coverage, assembly quality, and obtain BAM files for binning from megahit assembly
    """
    input:
        contigs = "results/{genera}/1_assembly/dedup_contigs_megahit/{sample}/{sample}_DEDUP95.fasta",
        r1 = "results/{genera}/1_assembly/dedup_reads/{sample}/{sample}_host_removed_dedup_R1.fastq",
        r2 = "results/{genera}/1_assembly/dedup_reads/{sample}/{sample}_host_removed_dedup_R2.fastq"
    output:
        "results/{genera}/1_assembly/contig_index_megahit/{sample}/{sample}_indexed_contig.1.bt2",
        "results/{genera}/1_assembly/contig_index_megahit/{sample}/{sample}_indexed_contig.2.bt2",
        "results/{genera}/1_assembly/contig_index_megahit/{sample}/{sample}_indexed_contig.3.bt2",
        "results/{genera}/1_assembly/contig_index_megahit/{sample}/{sample}_indexed_contig.4.bt2",
        "results/{genera}/1_assembly/contig_index_megahit/{sample}/{sample}_indexed_contig.rev.1.bt2",
        "results/{genera}/1_assembly/contig_index_megahit/{sample}/{sample}_indexed_contig.rev.2.bt2",
        "results/{genera}/1_assembly/contig_read_alignment_megahit/{sample}_aligned_sorted.bam"
    params:
        genera=config["genera"],
        outdir = "results/{genera}/1_assembly/contig_index_megahit/{sample}"
    log:
        stdout = "logs/{genera}/1_assembly/contig_read_alignment_megahit/{sample}_aln.out",
        stderr = "logs/{genera}/1_assembly/contig_read_alignment_megahit/{sample}_aln.err"
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

rule assembly_eval:
    """
    Evaluate assemblies
    """
    input:
        spades_contigs = "results/{genera}/1_assembly/dedup_contigs_spades/{sample}/{sample}_DEDUP95.fasta",
        megahit_contigs = "results/{genera}/1_assembly/dedup_contigs_megahit/{sample}/{sample}_DEDUP95.fasta"
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