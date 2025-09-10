# This pipeline goes through the following steps - 
    # 1-6. Indexing and alignment of PE read files against de-duplicated contigs from each assembly
    # 7. Quality evaluations of all assemblies

rule align_reads_to_contigs_spades_individual_assemblies:
    """
    Align reads to contigs assembled from individual metagenome assemblies
    """
    input:
        # Individual metagenome assemblies
        contigs = "results/{genera}/3_dedup_contigs/SPAdes_single/{sample}/{sample}_DEDUP95.fasta",
        r1 = "results/{genera}/1_pre_processing/dedup_reads/{sample}/{sample}_host_removed_dedup_R1.fastq",
        r2 = "results/{genera}/1_pre_processing/dedup_reads/{sample}/{sample}_host_removed_dedup_R2.fastq"
    output:
        "results/{genera}/4_align_reads_to_contigs/contig_index_spades_individual_assemblies/{sample}/{sample}_indexed_contig.1.bt2",
        "results/{genera}/4_align_reads_to_contigs/contig_index_spades_individual_assemblies/{sample}/{sample}_indexed_contig.2.bt2",
        "results/{genera}/4_align_reads_to_contigs/contig_index_spades_individual_assemblies/{sample}/{sample}_indexed_contig.3.bt2",
        "results/{genera}/4_align_reads_to_contigs/contig_index_spades_individual_assemblies/{sample}/{sample}_indexed_contig.4.bt2",
        "results/{genera}/4_align_reads_to_contigs/contig_index_spades_individual_assemblies/{sample}/{sample}_indexed_contig.rev.1.bt2",
        "results/{genera}/4_align_reads_to_contigs/contig_index_spades_individual_assemblies/{sample}/{sample}_indexed_contig.rev.2.bt2",
        "results/{genera}/4_align_reads_to_contigs/contig_read_alignment_individual_assemblies_spades/{sample}_aligned_sorted.bam"
    params:
        genera=config["genera"],
        outdir = "results/{genera}/4_align_reads_to_contigs/contig_index_spades_individual_assemblies/{sample}"
    log:
        stdout = "logs/{genera}/4_align_reads_to_contigs/contig_read_alignment_individual_assemblies_spades/{sample}_aln.out",
        stderr = "logs/{genera}/4_align_reads_to_contigs/contig_read_alignment_individual_assemblies_spades/{sample}_aln.err"
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

rule align_reads_to_contigs_spades_per_whale_coassemblies:
    """
    Align reads to contigs assembled from per whale metagenome co-assemblies
    """
    input:
        # Individual whale metagenome co-assemblies
        contigs = "results/{genera}/3_dedup_contigs/SPAdes_whales/{individual}/{individual}_DEDUP95.fasta",
        r3 = "results/{genera}/1_pre_processing/individual_coassembly_normalization/{individual}_ME_R1_norm.fq.gz",
        r4 = "results/{genera}/1_pre_processing/individual_coassembly_normalization/{individual}_ME_R2_norm.fq.gz"
    output:
        "results/{genera}/4_align_reads_to_contigs/contig_index_spades_per_whale_assemblies/{individual}/{individual}_indexed_contig.1.bt2",
        "results/{genera}/4_align_reads_to_contigs/contig_index_spades_per_whale_assemblies/{individual}/{individual}_indexed_contig.2.bt2",
        "results/{genera}/4_align_reads_to_contigs/contig_index_spades_per_whale_assemblies/{individual}/{individual}_indexed_contig.3.bt2",
        "results/{genera}/4_align_reads_to_contigs/contig_index_spades_per_whale_assemblies/{individual}/{individual}_indexed_contig.4.bt2",
        "results/{genera}/4_align_reads_to_contigs/contig_index_spades_per_whale_assemblies/{individual}/{individual}_indexed_contig.rev.1.bt2",
        "results/{genera}/4_align_reads_to_contigs/contig_index_spades_per_whale_assemblies/{individual}/{individual}_indexed_contig.rev.2.bt2",
        "results/{genera}/4_align_reads_to_contigs/contig_read_alignment_per_whale_assemblies_spades/{individual}_aligned_sorted.bam"
    params:
        genera=config["genera"],
        outdir = "results/{genera}/4_align_reads_to_contigs/contig_index_spades_per_whale_assemblies/{individual}"
    log:
        stdout = "logs/{genera}/4_align_reads_to_contigs/contig_read_alignment_per_whale_assemblies_spades/{individual}_aln.out",
        stderr = "logs/{genera}/4_align_reads_to_contigs/contig_read_alignment_per_whale_assemblies_spades/{individual}_aln.err"
    shell:
        """
        module unload miniconda 
        module load Bowtie2/2.5.1-GCC-12.2.0
        module load SAMtools/1.21-GCC-12.2.0

        # 1. Build contig index
        bowtie2-build \
        -f {input.contigs} {params.outdir}/{wildcards.individual}_indexed_contig \
        1>> {log.stdout} 2>> {log.stderr}

        # 2. Align reads back to assembled contigs
        bowtie2 \
        -x {params.outdir}/{wildcards.individual}_indexed_contig -1 {input.r1} -2 {input.r2} | samtools view -b -F 4 -F 2048 | samtools sort -o {output[6]} \
        1>> {log.stdout} 2>> {log.stderr}
        """

rule align_reads_to_contigs_spades_master_coassemblies:
    """
    Align reads to contigs assembled from individual metagenome assemblies
    """
    input:
        # Master metagenome co-assemblies
        contigs = "results/{genera}/3_dedup_contigs/SPAdes_master/DEDUP95.fasta",
        r1 = "results/{genera}/1_pre_processing/master_coassembly_normalization/all_samples_ME_R1_norm.fq.gz",
        r2 = "results/{genera}/1_pre_processing/master_coassembly_normalization/all_samples_ME_R2_norm.fq.gz"
    output:
        "results/{genera}/4_align_reads_to_contigs/contig_index_spades_master_coassemblies/indexed_contig.1.bt2",
        "results/{genera}/4_align_reads_to_contigs/contig_index_spades_master_coassemblies/indexed_contig.2.bt2",
        "results/{genera}/4_align_reads_to_contigs/contig_index_spades_master_coassemblies/indexed_contig.3.bt2",
        "results/{genera}/4_align_reads_to_contigs/contig_index_spades_master_coassemblies/indexed_contig.4.bt2",
        "results/{genera}/4_align_reads_to_contigs/contig_index_spades_master_coassemblies/indexed_contig.rev.1.bt2",
        "results/{genera}/4_align_reads_to_contigs/contig_index_spades_master_coassemblies/indexed_contig.rev.2.bt2",
        "results/{genera}/4_align_reads_to_contigs/contig_read_alignment_master_coassemblies_spades/aligned_sorted.bam"
    params:
        genera=config["genera"],
        outdir = "results/{genera}/4_align_reads_to_contigs/contig_index_spades_individual_assemblies"
    log:
        stdout = "logs/{genera}/4_align_reads_to_contigs/contig_read_alignment_master_coassemblies_spades/aln.out",
        stderr = "logs/{genera}/4_align_reads_to_contigs/contig_read_alignment_master_coassemblies_spades/aln.err"
    shell:
        """
        module unload miniconda 
        module load Bowtie2/2.5.1-GCC-12.2.0
        module load SAMtools/1.21-GCC-12.2.0

        # 1. Build contig index
        bowtie2-build \
        -f {input.contigs} {params.outdir}/indexed_contig \
        1>> {log.stdout} 2>> {log.stderr}

        # 2. Align reads back to assembled contigs
        bowtie2 \
        -x {params.outdir}/indexed_contig -1 {input.r1} -2 {input.r2} | samtools view -b -F 4 -F 2048 | samtools sort -o {output[6]} \
        1>> {log.stdout} 2>> {log.stderr}
        """

rule align_reads_to_contigs_megahit_individual_assemblies:
    """
    Align reads to contigs assembled from individual metagenome assemblies
    """
    input:
        # Individual metagenome assemblies
        contigs = "results/{genera}/3_dedup_contigs/megahit_single/{sample}/{sample}_DEDUP95.fasta",
        r1 = "results/{genera}/1_pre_processing/dedup_reads/{sample}/{sample}_host_removed_dedup_R1.fastq",
        r2 = "results/{genera}/1_pre_processing/dedup_reads/{sample}/{sample}_host_removed_dedup_R2.fastq"
    output:
        "results/{genera}/4_align_reads_to_contigs/contig_index_megahit_individual_assemblies/{sample}/{sample}_indexed_contig.1.bt2",
        "results/{genera}/4_align_reads_to_contigs/contig_index_megahit_individual_assemblies/{sample}/{sample}_indexed_contig.2.bt2",
        "results/{genera}/4_align_reads_to_contigs/contig_index_megahit_individual_assemblies/{sample}/{sample}_indexed_contig.3.bt2",
        "results/{genera}/4_align_reads_to_contigs/contig_index_megahit_individual_assemblies/{sample}/{sample}_indexed_contig.4.bt2",
        "results/{genera}/4_align_reads_to_contigs/contig_index_megahit_individual_assemblies/{sample}/{sample}_indexed_contig.rev.1.bt2",
        "results/{genera}/4_align_reads_to_contigs/contig_index_megahit_individual_assemblies/{sample}/{sample}_indexed_contig.rev.2.bt2",
        "results/{genera}/4_align_reads_to_contigs/contig_read_alignment_individual_assemblies_megahit/{sample}_aligned_sorted.bam"
    params:
        genera=config["genera"],
        outdir = "results/{genera}/4_align_reads_to_contigs/contig_index_megahit_individual_assemblies/{sample}"
    log:
        stdout = "logs/{genera}/4_align_reads_to_contigs/contig_read_alignment_individual_assemblies_megahit/{sample}_aln.out",
        stderr = "logs/{genera}/4_align_reads_to_contigs/contig_read_alignment_individual_assemblies_megahit/{sample}_aln.err"
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

rule align_reads_to_contigs_megahit_per_whale_coassemblies:
    """
    Align reads to contigs assembled from per whale metagenome co-assemblies
    """
    input:
        # Individual whale metagenome co-assemblies
        contigs = "results/{genera}/3_dedup_contigs/megahit_whales/{individual}/{individual}_DEDUP95.fasta",
        r3 = "results/{genera}/1_pre_processing/individual_coassembly_normalization/{individual}_ME_R1_norm.fq.gz",
        r4 = "results/{genera}/1_pre_processing/individual_coassembly_normalization/{individual}_ME_R2_norm.fq.gz"
    output:
        "results/{genera}/4_align_reads_to_contigs/contig_index_megahit_per_whale_assemblies/{individual}/{individual}_indexed_contig.1.bt2",
        "results/{genera}/4_align_reads_to_contigs/contig_index_megahit_per_whale_assemblies/{individual}/{individual}_indexed_contig.2.bt2",
        "results/{genera}/4_align_reads_to_contigs/contig_index_megahit_per_whale_assemblies/{individual}/{individual}_indexed_contig.3.bt2",
        "results/{genera}/4_align_reads_to_contigs/contig_index_megahit_per_whale_assemblies/{individual}/{individual}_indexed_contig.4.bt2",
        "results/{genera}/4_align_reads_to_contigs/contig_index_megahit_per_whale_assemblies/{individual}/{individual}_indexed_contig.rev.1.bt2",
        "results/{genera}/4_align_reads_to_contigs/contig_index_megahit_per_whale_assemblies/{individual}/{individual}_indexed_contig.rev.2.bt2",
        "results/{genera}/4_align_reads_to_contigs/contig_read_alignment_per_whale_assemblies_megahit/{individual}_aligned_sorted.bam"
    params:
        genera=config["genera"],
        outdir = "results/{genera}/4_align_reads_to_contigs/contig_index_megahit_per_whale_assemblies/{individual}"
    log:
        stdout = "logs/{genera}/4_align_reads_to_contigs/contig_read_alignment_per_whale_assemblies_megahit/{individual}_aln.out",
        stderr = "logs/{genera}/4_align_reads_to_contigs/contig_read_alignment_per_whale_assemblies_megahit/{individual}_aln.err"
    shell:
        """
        module unload miniconda 
        module load Bowtie2/2.5.1-GCC-12.2.0
        module load SAMtools/1.21-GCC-12.2.0

        # 1. Build contig index
        bowtie2-build \
        -f {input.contigs} {params.outdir}/{wildcards.individual}_indexed_contig \
        1>> {log.stdout} 2>> {log.stderr}

        # 2. Align reads back to assembled contigs
        bowtie2 \
        -x {params.outdir}/{wildcards.individual}_indexed_contig -1 {input.r1} -2 {input.r2} | samtools view -b -F 4 -F 2048 | samtools sort -o {output[6]} \
        1>> {log.stdout} 2>> {log.stderr}
        """

rule align_reads_to_contigs_megahit_master_coassemblies:
    """
    Align reads to contigs assembled from individual metagenome assemblies
    """
    input:
        # Master metagenome co-assemblies
        contigs = "results/{genera}/3_dedup_contigs/megahit_master/DEDUP95.fasta",
        r1 = "results/{genera}/1_pre_processing/master_coassembly_normalization/all_samples_ME_R1_norm.fq.gz",
        r2 = "results/{genera}/1_pre_processing/master_coassembly_normalization/all_samples_ME_R2_norm.fq.gz"
    output:
        "results/{genera}/4_align_reads_to_contigs/contig_index_megahit_master_coassemblies/indexed_contig.1.bt2",
        "results/{genera}/4_align_reads_to_contigs/contig_index_megahit_master_coassemblies/indexed_contig.2.bt2",
        "results/{genera}/4_align_reads_to_contigs/contig_index_megahit_master_coassemblies/indexed_contig.3.bt2",
        "results/{genera}/4_align_reads_to_contigs/contig_index_megahit_master_coassemblies/indexed_contig.4.bt2",
        "results/{genera}/4_align_reads_to_contigs/contig_index_megahit_master_coassemblies/indexed_contig.rev.1.bt2",
        "results/{genera}/4_align_reads_to_contigs/contig_index_megahit_master_coassemblies/indexed_contig.rev.2.bt2",
        "results/{genera}/4_align_reads_to_contigs/contig_read_alignment_master_coassemblies_megahit/aligned_sorted.bam"
    params:
        genera=config["genera"],
        outdir = "results/{genera}/4_align_reads_to_contigs/contig_index_megahit_individual_assemblies"
    log:
        stdout = "logs/{genera}/4_align_reads_to_contigs/contig_read_alignment_master_coassemblies_megahit/aln.out",
        stderr = "logs/{genera}/4_align_reads_to_contigs/contig_read_alignment_master_coassemblies_megahit/aln.err"
    shell:
        """
        module unload miniconda 
        module load Bowtie2/2.5.1-GCC-12.2.0
        module load SAMtools/1.21-GCC-12.2.0

        # 1. Build contig index
        bowtie2-build \
        -f {input.contigs} {params.outdir}/indexed_contig \
        1>> {log.stdout} 2>> {log.stderr}

        # 2. Align reads back to assembled contigs
        bowtie2 \
        -x {params.outdir}/indexed_contig -1 {input.r1} -2 {input.r2} | samtools view -b -F 4 -F 2048 | samtools sort -o {output[6]} \
        1>> {log.stdout} 2>> {log.stderr}
        """