# This pipeline goes through the following steps - 
    # 1-6. Indexing and alignment of PE read files against de-duplicated contigs from each assembly
    # 7. Quality evaluations of all assemblies

rule align_reads_to_individual_assemblies:
    """
    Align reads to contigs assembled from individual assemblies
    """
    input:
        contigs = "results/{genera}/testing/3_dedup_contigs/{assembler}/individual_metagenome_assembly/{sample}/{sample}_DEDUP95.fasta",
        r1 = "results/{genera}/1_pre_processing/dedup_reads/{sample}/{sample}_host_removed_dedup_R1.fastq",
        r2 = "results/{genera}/1_pre_processing/dedup_reads/{sample}/{sample}_host_removed_dedup_R2.fastq"
    output:
        "results/{genera}/testing/4_align_reads_to_contigs/{assembler}/contig_index_individual_assemblies/{sample}/{sample}_indexed_contig.1.bt2",
        "results/{genera}/testing/4_align_reads_to_contigs/{assembler}/contig_index_individual_assemblies/{sample}/{sample}_indexed_contig.2.bt2",
        "results/{genera}/testing/4_align_reads_to_contigs/{assembler}/contig_index_individual_assemblies/{sample}/{sample}_indexed_contig.3.bt2",
        "results/{genera}/testing/4_align_reads_to_contigs/{assembler}/contig_index_individual_assemblies/{sample}/{sample}_indexed_contig.4.bt2",
        "results/{genera}/testing/4_align_reads_to_contigs/{assembler}/contig_index_individual_assemblies/{sample}/{sample}_indexed_contig.rev.1.bt2",
        "results/{genera}/testing/4_align_reads_to_contigs/{assembler}/contig_index_individual_assemblies/{sample}/{sample}_indexed_contig.rev.2.bt2",
        "results/{genera}/testing/4_align_reads_to_contigs/{assembler}/contig_read_alignment_individual_assemblies/{sample}_aligned_sorted.bam",
        "results/{genera}/testing/4_align_reads_to_contigs/{assembler}/contig_read_alignment_individual_assemblies/{sample}_aligned_sorted.bam.bai"
    params:
        outdir = "results/{genera}/testing/4_align_reads_to_contigs/{assembler}/contig_index_individual_assemblies/{sample}",
        genera=config["genera"]
    log:
        stdout = "logs/{genera}/testing/4_align_reads_to_contigs/{assembler}/contig_read_alignment_individual_assemblies/{sample}_aln.out",
        stderr = "logs/{genera}/testing/4_align_reads_to_contigs/{assembler}/contig_read_alignment_individual_assemblies/{sample}_aln.err"
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

        # 3. Index BAM files to produce the bam.bai file that concoct requires later on
        samtools index \
        {output[6]} -o {output[7]} \
        1>> {log.stdout} 2>> {log.stderr}
        """