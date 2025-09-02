import os
import glob
import pandas as pd

# This pipeline contains 3 assemblies - 
    # 1. Individual assemblies where each sample undergoes its own assembly
        # individual_metagenome_assembly_spades & individual_metagenome_assembly_megahit
    # 2. Individual co-assemblies where reads are assembled by individual
        # individual_metagenome_coassembly_spades & individual_metagenome_coassembly_megahit
    # 3. Master co-assembly where all of our reads from all individuals are grouped and co-assembled
        # master_metagenome_coassembly_spades & master_metagenome_coassembly_megahit

rule individual_metagenome_assembly_spades:
    """
    Assemble contigs from individual samples using SPAdes
    """
    input:
        r1 = "results/{genera}/2_assembly/dedup_reads/{sample}/{sample}_host_removed_dedup_R1.fastq",
        r2 = "results/{genera}/2_assembly/dedup_reads/{sample}/{sample}_host_removed_dedup_R2.fastq"
    output:
        "results/{genera}/2_assembly/individual_metagenome_assembly_spades/{sample}/contigs.fasta"
    params:
        genera=config["genera"],
        outdir = "results/{genera}/2_assembly/individual_metagenome_assembly_spades/{sample}"
    log:
        stdout = "logs/{genera}/2_assembly/individual_metagenome_assembly_spades/{sample}/assembly.out",
        stderr = "logs/{genera}/2_assembly/individual_metagenome_assembly_spades/{sample}/assembly.err"
    resources:
        mem_mb=200000,
        threads=4
    shell:
        """
        module unload miniconda
        module load SPAdes/3.15.5-GCC-12.2.0

        spades.py \
        --meta --threads {resources.threads} \
        -1 {input.r1} -2 {input.r2} -o {params.outdir} \
        1> {log.stdout} 2> {log.stderr}
        """

rule individual_metagenome_assembly_megahit:
    """
    Assemble contigs from individual samples with Megahit
    """
    input:
        r1 = "results/{genera}/2_assembly/dedup_reads/{sample}/{sample}_host_removed_dedup_R1.fastq",
        r2 = "results/{genera}/2_assembly/dedup_reads/{sample}/{sample}_host_removed_dedup_R2.fastq"
    output:
        "results/{genera}/2_assembly/individual_metagenome_assembly_megahit/{sample}/final.contigs.fa"
    params:
        genera=config["genera"],
        preset="meta-large"
        outdir = "results/{genera}/2_assembly/individual_metagenome_assembly_megahit/{sample}"
    log:
        stdout = "logs/{genera}/2_assembly/individual_metagenome_assembly_megahit/{sample}/assembly.out",
        stderr = "logs/{genera}/2_assembly/individual_metagenome_assembly_megahit/{sample}/assembly.err"
    shell:
        """
        module unload miniconda
        source activate /home/flg9/.conda/envs/megahit

        megahit \
        -1 {input.r1} -2 {input.r2} -o {params.outdir} \
        --presets {params.preset} \
        1> {log.stdout} 2> {log.stderr}
        """

rule master_metagenome_coassembly_spades:
    """
    Generate a co-assembly from merged PE files of all samples using SPAdes
    """
    input:
        r1 = "results/{genera}/1_assembly/co_assembly_reads/all_samples_R1_norm.fq",
        r2 = "results/{genera}/1_assembly/co_assembly_reads/all_samples_R2_norm.fq"
    output:
        "results/{genera}/2_assembly/master_metagenome_coassembly_spades/{sample}/contigs.fasta"
    params:
        genera=config["genera"],
        outdir = "results/{genera}/2_assembly/master_metagenome_coassembly_spades/{sample}"
    log:
        stdout = "logs/{genera}/2_assembly/master_metagenome_coassembly_spades/{sample}/assembly.out",
        stderr = "logs/{genera}/2_assembly/master_metagenome_coassembly_spades/{sample}/assembly.err"
    resources:
        mem_mb=200000,
        threads=4
    shell:
        """
        module unload miniconda
        module load SPAdes/3.15.5-GCC-12.2.0

        spades.py \
        --meta --threads {resources.threads} \
        -1 {input.r1} -2 {input.r2} -o {params.outdir} \
        1> {log.stdout} 2> {log.stderr}
        """

rule master_metagenome_coassembly_megahit:
    """
    Generate a co-assembly from merged PE files of all samples using Megahit
    """
    input:
        r1 = "results/{genera}/1_assembly/co_assembly_reads/all_samples_R1_norm.fq",
        r2 = "results/{genera}/1_assembly/co_assembly_reads/all_samples_R2_norm.fq"
    output:
        "results/{genera}/2_assembly/master_metagenome_coassembly_megahit/{sample}/final.contigs.fa"
    params:
        genera=config["genera"],
        preset="meta-large"
        outdir = "results/{genera}/2_assembly/master_metagenome_coassembly_megahit/{sample}"
    log:
        stdout = "logs/{genera}/2_assembly/master_metagenome_coassembly_megahit/{sample}/assembly.out",
        stderr = "logs/{genera}/2_assembly/master_metagenome_coassembly_megahit/{sample}/assembly.err"
    shell:
        """
        module unload miniconda
        source activate /home/flg9/.conda/envs/megahit

        megahit \
        -1 {input.r1} -2 {input.r2} -o {params.outdir} \
        --presets {params.preset} \
        1> {log.stdout} 2> {log.stderr}
        """

rule individual_metagenome_coassembly_spades:
    """
    Generate a co-assembly using spades
    """
    input:
        norm1 = "results/{genera}/1_QC/individual_coassembly_normalization/{individual}_ME_R1_norm.fq.gz",
        norm2 = "results/{genera}/1_QC/individual_coassembly_normalization/{individual}_ME_R2_norm.fq.gz"
    output:
        "results/{genera}/2_assembly/individual_metagenome_coassembly_spades/{sample}/contigs.fasta"
    params:
        genera=config["genera"],
        outdir = "results/{genera}/2_assembly/individual_metagenome_coassembly_spades/{sample}"
    log:
        stdout = "logs/{genera}/2_assembly/individual_metagenome_coassembly_spades/{sample}/assembly.out",
        stderr = "logs/{genera}/2_assembly/individual_metagenome_coassembly_spades/{sample}/assembly.err"
    resources:
        mem_mb=200000,
        threads=4
    shell:
        """
        module unload miniconda
        module load SPAdes/3.15.5-GCC-12.2.0

        spades.py \
        --meta --threads {resources.threads} \
        -1 {input.r1} -2 {input.r2} -o {params.outdir} \
        1> {log.stdout} 2> {log.stderr}
        """

rule indvidual_metagenome_coassembly_megahit:
    """
    """
    input:
        norm1 = "results/{genera}/1_QC/individual_coassembly_normalization/{individual}_ME_R1_norm.fq.gz",
        norm2 = "results/{genera}/1_QC/individual_coassembly_normalization/{individual}_ME_R2_norm.fq.gz"
    output:
        "results/{genera}/2_assembly/indvidual_metagenome_coassembly_megahit/{sample}/final.contigs.fa"
    params:
        genera=config["genera"],
        preset="meta-large"
        outdir = "results/{genera}/2_assembly/indvidual_metagenome_coassembly_megahit/{sample}"
    log:
        stdout = "logs/{genera}/2_assembly/indvidual_metagenome_coassembly_megahit/{sample}/assembly.out",
        stderr = "logs/{genera}/2_assembly/indvidual_metagenome_coassembly_megahit/{sample}/assembly.err"
    shell:
        """
        module unload miniconda
        source activate /home/flg9/.conda/envs/megahit

        megahit \
        -1 {input.r1} -2 {input.r2} -o {params.outdir} \
        --presets {params.preset} \
        1> {log.stdout} 2> {log.stderr}
        """

rule deduplicate_contigs:
    """
    Run CD-HIT to deduplicate contigs that may have been introduced during metagenome assembly
    """
    input:
        SPAdes_master = "results/{genera}/2_assembly/master_metagenome_coassembly_spades/{sample}/contigs.fasta",
        SPAdes_indiv = "results/{genera}/2_assembly/individual_metagenome_coassembly_spades/{sample}/contigs.fasta",
        SPAdes_single = "results/{genera}/2_assembly/individual_metagenome_assembly_spades/{sample}/contigs.fasta",
        megahit_master = "results/{genera}/2_assembly/master_metagenome_coassembly_megahit/{sample}/final.contigs.fa",
        megahit_indiv = "results/{genera}/2_assembly/indvidual_metagenome_coassembly_megahit/{sample}/final.contigs.fa",
        megahit_single = "results/{genera}/2_assembly/individual_metagenome_assembly_megahit/{sample}/final.contigs.fa",
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

        #SPAdes deduplication 

        # 1. Deduplicate contigs from SPAdes single assemblies
        cd-hit-est \
        -i {inpu.SPAdes_single} -o {out1.dedups} -T 4 -c 0.95 \
        1>> {log.stdout} 2>> {log.stderr}

        # 2. Deduplicate contigs from SPAdes per individual assemblies
        cd-hit-est \
        -i {inpu.SPAdes_indiv} -o {out1.dedups} -T 4 -c 0.95 \
        1>> {log.stdout} 2>> {log.stderr}

        # 3. Deduplicate contigs from SPAdes master assembly
        cd-hit-est \
        -i {inpu.SPAdes_master} -o {out1.dedups} -T 4 -c 0.95 \
        1>> {log.stdout} 2>> {log.stderr}

        # Megahit deduplication 

        # 1. Deduplicate contigs from Megahit assembly
        cd-hit-est \
        -i {input.megahit_assembly} -o {out2.dedups} -T 4 -c 0.95 \
        1>> {log.stdout} 2>> {log.stderr}

        # 2. Deduplicate contigs from Megahit assembly
        cd-hit-est \
        -i {input.megahit_assembly} -o {out2.dedups} -T 4 -c 0.95 \
        1>> {log.stdout} 2>> {log.stderr}

        # 3. Deduplicate contigs from Megahit assembly
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