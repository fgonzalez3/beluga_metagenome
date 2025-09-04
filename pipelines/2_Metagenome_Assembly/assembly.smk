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

rule individual_metagenome_assembly_spades: # done
    """
    Assemble contigs from individual samples using SPAdes
    """
    input:
        r1 = "results/{genera}/1_pre_processing/dedup_reads/{sample}/{sample}_host_removed_dedup_R1.fastq",
        r2 = "results/{genera}/1_pre_processing/dedup_reads/{sample}/{sample}_host_removed_dedup_R2.fastq"
    output:
        "results/{genera}/2_assembly/SPAdes/individual_metagenome_assembly/{sample}/contigs.fasta"
    params:
        genera=config["genera"],
        outdir = "results/{genera}/2_assembly/SPAdes/individual_metagenome_assembly/{sample}"
    log:
        stdout = "logs/{genera}/2_assembly/SPAdes/individual_metagenome_assembly/{sample}/assembly.out",
        stderr = "logs/{genera}/2_assembly/SPAdes/individual_metagenome_assembly{sample}/assembly.err"
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

rule individual_metagenome_assembly_megahit: # test
    """
    Assemble contigs from individual samples with Megahit
    """
    input:
        r1 = "results/{genera}/1_pre_processing/dedup_reads/{sample}/{sample}_host_removed_dedup_R1.fastq",
        r2 = "results/{genera}/1_pre_processing/dedup_reads/{sample}/{sample}_host_removed_dedup_R2.fastq"
    output:
        "results/{genera}/2_assembly/megahit/individual_metagenome_assembly/{sample}/final.contigs.fa"
    params:
        genera=config["genera"],
        preset="meta-large"
        outdir = "results/{genera}/2_assembly/megahit/individual_metagenome_assembly/{sample}"
    log:
        stdout = "logs/{genera}/2_assembly/megahit/individual_metagenome_assembly/{sample}/assembly.out",
        stderr = "logs/{genera}/2_assembly/megahit/individual_metagenome_assembly/{sample}/assembly.err"
    shell:
        """
        module unload miniconda
        source activate /home/flg9/.conda/envs/megahit

        megahit \
        -1 {input.r1} -2 {input.r2} -o {params.outdir} \
        --presets {params.preset} \
        1> {log.stdout} 2> {log.stderr}
        """

rule master_metagenome_coassembly_spades: # test
    """
    Generate a master co-assembly from merged PE files of all samples using SPAdes
    """
    input:
        r1 = "results/{genera}/1_pre_processing/master_coassembly_normalization/all_samples_ME_R1_norm.fq.gz",
        r2 = "results/{genera}/1_pre_processing/master_coassembly_normalization/all_samples_ME_R2_norm.fq.gz"
    output:
        "results/{genera}/2_assembly/SPAdes/master_metagenome_coassembly/{sample}/contigs.fasta"
    params:
        genera=config["genera"],
        outdir = "results/{genera}/2_assembly/SPAdes/master_metagenome_coassembly/{sample}"
    log:
        stdout = "logs/{genera}/2_assembly/SPAdes/master_metagenome_coassembly/{sample}/assembly.out",
        stderr = "logs/{genera}/2_assembly/SPAdes/master_metagenome_coassembly/{sample}/assembly.err"
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

rule master_metagenome_coassembly_megahit: # test
    """
    Generate a co-assembly from merged PE files of all samples using Megahit
    """
    input:
        r1 = "results/{genera}/1_pre_processing/master_coassembly_normalization/all_samples_ME_R1_norm.fq.gz",
        r2 = "results/{genera}/1_pre_processing/master_coassembly_normalization/all_samples_ME_R2_norm.fq.gz"
    output:
        "results/{genera}/2_assembly/megahit/master_metagenome_coassembly/{sample}/final.contigs.fa"
    params:
        genera=config["genera"],
        preset="meta-large"
        outdir = "results/{genera}/megahit/2_assembly/master_metagenome_coassembly/{sample}"
    log:
        stdout = "logs/{genera}/2_assembly/megahit/master_metagenome_coassembly/{sample}/assembly.out",
        stderr = "logs/{genera}/2_assembly/megahit/master_metagenome_coassembly/{sample}/assembly.err"
    shell:
        """
        module unload miniconda
        source activate /home/flg9/.conda/envs/megahit

        megahit \
        -1 {input.r1} -2 {input.r2} -o {params.outdir} \
        --presets {params.preset} \
        1> {log.stdout} 2> {log.stderr}
        """

# Define individuals for this next rule 

INDIVIDUALS = ["JUNO", "KELA", "NATTY"]

def indiv_r1(wildcards):
    return sorted(glob.glob(f"results/{wildcards.genera}/1_pre_processing/dedup_reads/{sample}/Sample_{wildcards.individual}_*_host_removed_dedup_R1.fastq"))

def indiv_r2(wildcards):
    return sorted(glob.glob(f"results/{wildcards.genera}/1_pre_processing/dedup_reads/{sample}/Sample_{wildcards.individual}_*_host_removed_dedup_R2.fastq"))

rule individual_metagenome_coassembly_spades: # test
    """
    Generate a co-assembly using spades per individual
    """
    input:
        r1 = "results/{genera}/1_pre_processing/individual_coassembly_normalization/{individual}_ME_R1_norm.fq.gz",
        r2 = "results/{genera}/1_pre_processing/individual_coassembly_normalization/{individual}_ME_R2_norm.fq.gz"
    output:
        "results/{genera}/2_assembly/SPAdes/individual_metagenome_coassembly/{sample}/contigs.fasta"
    params:
        genera=config["genera"],
        outdir = "results/{genera}/2_assembly/SPAdes/individual_metagenome_coassembly/{sample}"
    log:
        stdout = "logs/{genera}/2_assembly/SPAdes/individual_metagenome_coassembly/{sample}/assembly.out",
        stderr = "logs/{genera}/2_assembly/SPAdes/individual_metagenome_coassembly/{sample}/assembly.err"
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

rule indvidual_metagenome_coassembly_megahit: # test
    """
    Generate a co-assembly using megahit per individual
    """
    input:
        r1 = "results/{genera}/1_pre_processing/individual_coassembly_normalization/{individual}_ME_R1_norm.fq.gz",
        r2 = "results/{genera}/1_pre_processing/individual_coassembly_normalization/{individual}_ME_R2_norm.fq.gz"
    output:
        "results/{genera}/2_assembly/megahit/indvidual_metagenome_coassembly/{sample}/final.contigs.fa"
    params:
        genera=config["genera"],
        preset="meta-large"
        outdir = "results/{genera}/2_assembly/megahit/indvidual_metagenome_coassembly/{sample}"
    log:
        stdout = "logs/{genera}/2_assembly/megahit/indvidual_metagenome_coassembly/{sample}/assembly.out",
        stderr = "logs/{genera}/2_assembly/megahit/indvidual_metagenome_coassembly/{sample}/assembly.err"
    shell:
        """
        module unload miniconda
        source activate /home/flg9/.conda/envs/megahit

        megahit \
        -1 {input.r1} -2 {input.r2} -o {params.outdir} \
        --presets {params.preset} \
        1> {log.stdout} 2> {log.stderr}
        """

rule deduplicate_contigs: # test 
    """
    Run CD-HIT to deduplicate overrepresented contigs that may have been introduced during metagenome assembly
    """
    input:
        SPAdes_master = "results/{genera}/2_assembly/SPAdes/master_metagenome_coassembly/{sample}/contigs.fasta",
        SPAdes_whales = "results/{genera}/2_assembly/SPAdes/individual_metagenome_coassembly/{sample}/contigs.fasta",
        SPAdes_single = "results/{genera}/2_assembly/SPAdes/individual_metagenome_assembly/{sample}/contigs.fasta",
        megahit_master = "results/{genera}/2_assembly/megahit/master_metagenome_coassembly/{sample}/final.contigs.fa",
        megahit_whales = "results/{genera}/2_assembly/megahit/indvidual_metagenome_coassembly/{sample}/final.contigs.fa",
        megahit_single = "results/{genera}/2_assembly/megahit/individual_metagenome_assembly/{sample}/final.contigs.fa",
    output:
        SPAdes_single_dedups = "results/{genera}/2_assembly/dedup_contigs/SPAdes_single/{sample}/{sample}_DEDUP95.fasta",
        SPAdes_master_dedups = "results/{genera}/2_assembly/dedup_contigs/SPAdes_master/{sample}/{sample}_DEDUP95.fasta",
        SPAdes_whales_dedups = "results/{genera}/2_assembly/dedup_contigs/SPAdes_whales/{sample}/{sample}_DEDUP95.fasta",
        megahit_single_dedups = "results/{genera}/2_assembly/dedup_contigs/megahit_single/{sample}/{sample}_DEDUP95.fasta",
        megahit_master_dedups = "results/{genera}/2_assembly/dedup_contigs/megahit_master/{sample}/{sample}_DEDUP95.fasta",
        megahit_whales_dedups = "results/{genera}/2_assembly/dedup_contigs/megahit_whales/{sample}/{sample}_DEDUP95.fasta"
    params:
        genera=config["genera"],
        out1 = "results/{genera}/2_assembly/dedup_contigs_spades/{sample}",
        out2 = "results/{genera}/2_assembly/dedup_contigs/SPAdes_master/{sample}",
        out3 = "results/{genera}/2_assembly/dedup_contigs/SPAdes_whales/{sample}",
        out4 = "results/{genera}/2_assembly/dedup_contigs/megahit_single/{sample}",
        out5 = "results/{genera}/2_assembly/dedup_contigs/megahit_master/{sample}",
        out6 = "results/{genera}/2_assembly/dedup_contigs/megahit_whales/{sample}"
    log:
        stdout = "logs/{genera}/2_assembly/dedup_contigs/SPAdes_single/{sample}/dedup_contigs.out",
        stderr = "logs/{genera}/2_assembly/dedup_contigs/{sample}/dedup_contigs.err"
    shell:
        """
        module unload miniconda
        module load CD-HIT/4.8.1-GCC-10.2.0

        # Create output directories 
        mkdir -p {params.out1}
        mkdir -p {params.out2}
        mkdir -p {params.out3}
        mkdir -p {params.out4}
        mkdir -p {params.out5}
        mkdir -p {params.out6}

        # SPAdes de-duplication 

        # 1. Deduplicate contigs from SPAdes single assemblies
        cd-hit-est \
        -i {input.SPAdes_single} -o {output.SPAdes_single_dedups} -T 4 -c 0.95 \
        1>> {log.stdout} 2>> {log.stderr}

        # 2. Deduplicate contigs from SPAdes per whale assemblies
        cd-hit-est \
        -i {inpu.SPAdes_indiv} -o {output.SPAdes_whales_dedups} -T 4 -c 0.95 \
        1>> {log.stdout} 2>> {log.stderr}

        # 3. Deduplicate contigs from SPAdes master assembly
        cd-hit-est \
        -i {inpu.SPAdes_master} -o {output.SPAdes_master_dedups} -T 4 -c 0.95 \
        1>> {log.stdout} 2>> {log.stderr}

        # Megahit de-duplication 

        # 1. Deduplicate contigs from Megahit single assemblies
        cd-hit-est \
        -i {input.megahit_assembly} -o {output.megahit_single_dedups} -T 4 -c 0.95 \
        1>> {log.stdout} 2>> {log.stderr}

        # 2. Deduplicate contigs from Megahit per whale assemblies
        cd-hit-est \
        -i {input.megahit_assembly} -o {out2.megahit_whales_dedups} -T 4 -c 0.95 \
        1>> {log.stdout} 2>> {log.stderr}

        # 3. Deduplicate contigs from Megahit master assemblies
        cd-hit-est \
        -i {input.megahit_assembly} -o {out2.megahit_master_dedups} -T 4 -c 0.95 \
        1>> {log.stdout} 2>> {log.stderr}
        """

rule align_reads_to_contigs_spades: # test
    """
    Align raw reads back to assembled contigs to assess coverage, assembly quality, and obtain BAM files for binning from SPAdes assembly.
    BAM files will consist of read alignments to individual sample assemblies and co-assemblies. 
    """
    input:
        # 1. Individual metagenome assemblies
        c1 = "results/{genera}/2_assembly/dedup_contigs_spades/{sample}/{sample}_DEDUP95.fasta",
        r1 = "results/{genera}/1_pre_processing/dedup_reads/{sample}/{sample}_host_removed_dedup_R1.fastq",
        r2 = "results/{genera}/1_pre_processing/dedup_reads/{sample}/{sample}_host_removed_dedup_R2.fastq",

        # 2. Per whale metagenome co-assemblies
        c2 = "results/{genera}/2_assembly/dedup_contigs/SPAdes_whales/{sample}/{sample}_DEDUP95.fasta",
        r3 = "results/{genera}/1_pre_processing/individual_coassembly_normalization/{individual}_ME_R1_norm.fq.gz",
        r4 = "results/{genera}/1_pre_processing/individual_coassembly_normalization/{individual}_ME_R2_norm.fq.gz",

        # 3. Master co-assembly
        c3 = "results/{genera}/2_assembly/SPAdes/master_metagenome_coassembly/{sample}/contigs.fasta",
        r5 = "results/{genera}/1_pre_processing/master_coassembly_normalization/all_samples_ME_R1_norm.fq.gz",
        r6 = "results/{genera}/1_pre_processing/master_coassembly_normalization/all_samples_ME_R2_norm.fq.gz"
    output:
        # 1. Individual metagenome assemblies
        "results/{genera}/2_assembly/align_reads_to_contigs/contig_index_spades_individual_assemblies/{sample}/{sample}_indexed_contig.1.bt2",
        "results/{genera}/2_assembly/align_reads_to_contigs/contig_index_spades_individual_assemblies/{sample}/{sample}_indexed_contig.2.bt2",
        "results/{genera}/2_assembly/align_reads_to_contigs/contig_index_spades_individual_assemblies/{sample}/{sample}_indexed_contig.3.bt2",
        "results/{genera}/2_assembly/align_reads_to_contigs/contig_index_spades_individual_assemblies/{sample}/{sample}_indexed_contig.4.bt2",
        "results/{genera}/2_assembly/align_reads_to_contigs/contig_index_spades_individual_assemblies/{sample}/{sample}_indexed_contig.rev.1.bt2",
        "results/{genera}/2_assembly/align_reads_to_contigs/contig_index_spades_individual_assemblies/{sample}/{sample}_indexed_contig.rev.2.bt2",
        "results/{genera}/2_assembly/align_reads_to_contigs/contig_read_alignment_spades/{sample}_aligned_sorted.bam",

        # 2. Per whale metagenome co-assemblies
        "results/{genera}/2_assembly/align_reads_to_contigs/contig_index_spades_per_whale_coassemblies/{sample}/{sample}_indexed_contig.1.bt2",
        "results/{genera}/2_assembly/align_reads_to_contigs/contig_index_spades_per_whale_coassemblies/{sample}/{sample}_indexed_contig.2.bt2",
        "results/{genera}/2_assembly/align_reads_to_contigs/contig_index_spades_per_whale_coassemblies/{sample}/{sample}_indexed_contig.3.bt2",
        "results/{genera}/2_assembly/align_reads_to_contigs/contig_index_spades_per_whale_coassemblies/{sample}/{sample}_indexed_contig.4.bt2",
        "results/{genera}/2_assembly/align_reads_to_contigs/contig_index_spades_per_whale_coassemblies/{sample}/{sample}_indexed_contig.rev.1.bt2",
        "results/{genera}/2_assembly/align_reads_to_contigs/contig_index_spades_per_whale_coassemblies/{sample}/{sample}_indexed_contig.rev.2.bt2",
        "results/{genera}/2_assembly/align_reads_to_contigs/contig_index_spades_per_whale_coassemblies/{sample}_aligned_sorted.bam",

        # 3. Master co-assembly
        "results/{genera}/2_assembly/align_reads_to_contigs/contig_index_spades_master_coassemblies/{sample}/{sample}_indexed_contig.1.bt2",
        "results/{genera}/2_assembly/align_reads_to_contigs/contig_index_spades_master_coassemblies/{sample}/{sample}_indexed_contig.2.bt2",
        "results/{genera}/2_assembly/align_reads_to_contigs/contig_index_spades_master_coassemblies/{sample}/{sample}_indexed_contig.3.bt2",
        "results/{genera}/2_assembly/align_reads_to_contigs/contig_index_spades_master_coassemblies/{sample}/{sample}_indexed_contig.4.bt2",
        "results/{genera}/2_assembly/align_reads_to_contigs/contig_index_spades_master_coassemblies/{sample}/{sample}_indexed_contig.rev.1.bt2",
        "results/{genera}/2_assembly/align_reads_to_contigs/contig_index_spades_master_coassemblies/{sample}/{sample}_indexed_contig.rev.2.bt2",
        "results/{genera}/2_assembly/align_reads_to_contigs/contig_index_spades_master_coassemblies/{sample}_aligned_sorted.bam"
    params:
        genera=config["genera"],
        outdir1 = "results/{genera}/2_assembly/align_reads_to_contigs/contig_index_spades_individual_assemblies/{sample}",
        outdir2 = "results/{genera}/2_assembly/align_reads_to_contigs/contig_index_spades_per_whale_coassemblies/{sample}",
        outdir3 = "results/{genera}/2_assembly/align_reads_to_contigs/contig_index_spades_master_coassemblies/{sample}"
    log:
        stdout1 = "logs/{genera}/2_assembly/align_reads_to_contigs/contig_index_spades_individual_assemblies/{sample}_aln.out",
        stderr1 = "logs/{genera}/2_assembly/align_reads_to_contigs/contig_index_spades_individual_assemblies/{sample}_aln.err",

        stdout2 = "logs/{genera}/2_assembly/align_reads_to_contigs/contig_index_spades_per_whale_coassemblies/{sample}_aln.out",
        stderr2 = "logs/{genera}/2_assembly/align_reads_to_contigs/contig_index_spades_per_whale_coassemblies/{sample}_aln.err",

        stdout3 = "logs/{genera}/2_assembly/align_reads_to_contigs/contig_index_spades_master_coassemblies/{sample}_aln.out",
        stderr3 = "logs/{genera}/2_assembly/align_reads_to_contigs/contig_index_spades_master_coassemblies/{sample}_aln.err"
    shell:
        """
        module unload miniconda 
        module load Bowtie2/2.5.1-GCC-12.2.0
        module load SAMtools/1.21-GCC-12.2.0

        # 1. Individual metagenome assemblies

        # a. Build contig index
        bowtie2-build \
        -f {input.c1} {params.outdir1}/{wildcards.sample}_indexed_contig \
        1>> {log.stdout1} 2>> {log.stderr1}

        # b. Align reads back to assembled contigs
        bowtie2 \
        -x {params.outdir1}/{wildcards.sample}_indexed_contig -1 {input.r1} -2 {input.r2} | samtools view -b -F 4 -F 2048 | samtools sort -o {output[6]} \
        1>> {log.stdout1} 2>> {log.stderr1}

        # 2. Per whale metagenome co-assemblies

        # a. Build contig index
        bowtie2-build \
        -f {input.c2} {params.outdir2}/{wildcards.sample}_indexed_contig \
        1>> {log.stdout2} 2>> {log.stderr2}

        # b. Align reads back to assembled contigs
        bowtie2 \
        -x {params.outdir2}/{wildcards.sample}_indexed_contig -1 {input.r3} -2 {input.r4} | samtools view -b -F 4 -F 2048 | samtools sort -o {output[13]} \
        1>> {log.stdout2} 2>> {log.stderr2}

        # 3. Master co-assembly

        # a. Build contig index
        bowtie2-build \
        -f {input.c3} {params.outdir3}/{wildcards.sample}_indexed_contig \
        1>> {log.stdout3} 2>> {log.stderr3}

        # b. Align reads back to assembled contigs
        bowtie2 \
        -x {params.outdir3}/{wildcards.sample}_indexed_contig -1 {input.r5} -2 {input.r6} | samtools view -b -F 4 -F 2048 | samtools sort -o {output[20]} \
        1>> {log.stdout3} 2>> {log.stderr3}
        """

rule align_reads_to_contigs_megahit: # test
    """
    Align raw reads back to assembled contigs to assess coverage, assembly quality, and obtain BAM files for binning from megahit assembly
    """
    input:
        contigs = "results/{genera}/2_assembly/dedup_contigs_megahit/{sample}/{sample}_DEDUP95.fasta",
        r1 = "results/{genera}/2_assembly/dedup_reads/{sample}/{sample}_host_removed_dedup_R1.fastq",
        r2 = "results/{genera}/2_assembly/dedup_reads/{sample}/{sample}_host_removed_dedup_R2.fastq"
    output:
        "results/{genera}/2_assembly/contig_index_megahit/{sample}/{sample}_indexed_contig.1.bt2",
        "results/{genera}/2_assembly/contig_index_megahit/{sample}/{sample}_indexed_contig.2.bt2",
        "results/{genera}/2_assembly/contig_index_megahit/{sample}/{sample}_indexed_contig.3.bt2",
        "results/{genera}/2_assembly/contig_index_megahit/{sample}/{sample}_indexed_contig.4.bt2",
        "results/{genera}/2_assembly/contig_index_megahit/{sample}/{sample}_indexed_contig.rev.1.bt2",
        "results/{genera}/2_assembly/contig_index_megahit/{sample}/{sample}_indexed_contig.rev.2.bt2",
        "results/{genera}/2_assembly/contig_read_alignment_megahit/{sample}_aligned_sorted.bam"
    params:
        genera=config["genera"],
        outdir = "results/{genera}/2_assembly/contig_index_megahit/{sample}"
    log:
        stdout = "logs/{genera}/2_assembly/contig_read_alignment_megahit/{sample}_aln.out",
        stderr = "logs/{genera}/2_assembly/contig_read_alignment_megahit/{sample}_aln.err"
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
        spades_contigs = "results/{genera}/2_assembly/dedup_contigs_spades/{sample}/{sample}_DEDUP95.fasta",
        megahit_contigs = "results/{genera}/2_assembly/dedup_contigs_megahit/{sample}/{sample}_DEDUP95.fasta"
    output:
        filtered_contigs = "results/{genera}/2_assembly/assembly_eval/{sample}/metaspades_assembly_DEDUP95_m1500.fasta",
        whole_assembly_stats = "results/{genera}/2_assembly/assembly_eval/{sample}/assembly_stats.csv",
        metaquast_output = "results/{genera}/2_assembly/assembly_eval/{sample}/report.html"
    params:
        len = 1500,
        threads = 4,
        outdir = "results/{genera}/2_assembly/assembly_eval/{sample}",
        outfile_label = "SPAdes.m1500"
    log:
        stdout = "logs/{genera}/2_assembly/assembly_eval/{sample}/assembly_eval.out",
        stderr = "logs/{genera}/2_assembly/assembly_eval/{sample}/assembly_eval.err"
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