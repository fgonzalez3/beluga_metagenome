# This pipeline goes through the following steps - 
    # 1. Adapter and low read quality trimming
    # 2. Aggregation of quality reports produced from step 1 where you can view and compare all reports for PE read file
    # 3-6. Masking of both human and beluga whale reference genomes, as well as PE read alignment to masked genomes to filter out host reads
    # 7. Deduplication of PE read files 
    # 8. Normalization of PE read files for PE reads that will undergo co-assembly of any kind

rule adapter_trimming: # done 
    """
    Trim adapters and low-quality reads
    """
    input:
        fwd = lambda wildcards: READS[wildcards.sample]["r1"],
        rev = lambda wildcards: READS[wildcards.sample]["r2"]
    output:
        "results/{genera}/1_pre_processing/adapter_trimming/{sample}/{sample}_val_1.fq.gz",
        "results/{genera}/1_pre_processing/adapter_trimming/{sample}/{sample}_val_2.fq.gz",
        "results/{genera}/1_pre_processing/adapter_trimming/{sample}/{sample}_val_1_fastqc.html",
        "results/{genera}/1_pre_processing/adapter_trimming/{sample}/{sample}_val_2_fastqc.html"
    params:
        genera = config["genera"],
        outdir = "results/{genera}/1_pre_processing/adapter_trimming/{sample}"
    log:
        stdout = "logs/{genera}/1_pre_processing/adapter_trimming/{sample}/qc.out",
        stderr = "logs/{genera}/1_pre_processing/adapter_trimming/{sample}/qc.err"
    shell:
        """
        module unload miniconda
        module load Trim_Galore/0.6.7-GCCcore-10.2.0

        trim_galore \
        --paired --fastqc --length 125 --cores 8 -q 25 {input.fwd} {input.rev} \
        --output_dir {params.outdir} --basename {wildcards.sample} \
        1> {log.stdout} 2> {log.stderr}
        """

rule aggregate_qc_data: # done
    """
    Aggregate QC data to create a single report across our many samples
    """
    input:
        expand("results/{genera}/1_pre_processing/adapter_trimming/{sample}/{sample}_val_1_fastqc.html", genera=config["genera"], sample=SAMPLES),
        expand("results/{genera}/1_pre_processing/adapter_trimming/{sample}/{sample}_val_2_fastqc.html", genera=config["genera"], sample=SAMPLES)
    output:
        "results/{genera}/1_pre_processing/aggregate_qc_data/multiqc_report.html"
    params:
        input_dirs = lambda wildcards: f"results/{wildcards.genera}/1_pre_processing/adapter_trimming",
        outdir = "results/{genera}/1_pre_processing/aggregate_qc_data"
    log:
        stdout = "logs/{genera}/1_pre_processing/aggregate_qc_data/multiqc.out",
        stderr = "logs/{genera}/1_pre_processing/aggregate_qc_data/multiqc.err"
    shell:
        """
        module unload miniconda 
        module load MultiQC/1.10.1-foss-2020b-Python-3.8.6

        multiqc \
        {params.input_dirs} -o {params.outdir} \
        1> {log.stdout} 2> {log.stderr}
        """

rule mask_beluga_host_genome: # done
    """
    Mask low-complexity or microbial contaminant regions in a beluga reference sequence, reducing false positives during host read filtering
    """
    input:
        beluga_refseq = config["beluga_refseq"]
    output:
        beluga_ref_masked = "results/{genera}/1_pre_processing/mask_beluga_host_genome/beluga_genome_sequence_masked.fa"
    params:
        genera=config["genera"],
        threads = 4
    log:
        stdout = "logs/{genera}/1_pre_processing/mask_beluga_host_genome/bbmask.out",
        stderr = "logs/{genera}/1_pre_processing/mask_beluga_host_genome/bbmask.err"
    shell:
        """
        module unload miniconda
        module load BBMap/38.90-GCCcore-10.2.0

        bbmask.sh \
        in={input.beluga_refseq} out={output.beluga_ref_masked} threads={params.threads} \
        1> {log.stdout} 2> {log.stderr}
        """

rule align_to_beluga_host_genome: # done 
    """
    Filter reads against masked beluga reference genome
    """
    input:
        beluga_ref_masked = "results/{genera}/1_pre_processing/mask_beluga_host_genome/beluga_genome_sequence_masked.fa",
        r1 = "results/{genera}/1_pre_processing/adapter_trimming/{sample}/{sample}_val_1.fq.gz",
        r2 = "results/{genera}/1_pre_processing/adapter_trimming/{sample}/{sample}_val_2.fq.gz"
    output:
        mapped = "results/{genera}/1_pre_processing/align_to_beluga_host_genome/{sample}/mapped.fq",
        unmapped = "results/{genera}/1_pre_processing/align_to_beluga_host_genome/{sample}/unmapped.fq"
    params:
        index_path = "results/{genera}/1_pre_processing/align_to_beluga_host_genome/index/",
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
        stdout = "logs/{genera}/1_pre_processing/align_to_beluga_host_genome/{sample}/bbmap_beluga.out",
        stderr = "logs/{genera}/1_pre_processing/align_to_beluga_host_genome/{sample}/bbmap_beluga.err"
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

rule split_unmapped_reads_beluga: # done
    """
    Split interleaved reads that did not map to masked beluga reference into separate fastqs
    """
    input:
        unmapped = "results/{genera}/1_pre_processing/align_to_beluga_host_genome/{sample}/unmapped.fq"
    output:
        r1 = "results/{genera}/1_pre_processing/split_unmapped_reads_beluga/{sample}/unmapped_R1.fq",
        r2 = "results/{genera}/1_pre_processing/split_unmapped_reads_beluga/{sample}/unmapped_R2.fq"
    params:
        genera=config["genera"]
    log:
        stdout = "logs/{genera}/1_pre_processing/split_unmapped_reads_beluga/{sample}/reformat_beluga.out",
        stderr = "logs/{genera}/1_pre_processing/split_unmapped_reads_beluga/{sample}/reformat_beluga.err"
    shell:
        """
        module unload miniconda
        module load BBMap/38.90-GCCcore-10.2.0

        reformat.sh \
        in={input.unmapped} \
        out1={output.r1} \
        out2={output.r2} \
        1> {log.stdout} 2> {log.stderr}
        """

rule mask_human_host_genome: # done
    """
    Mask low-complexity or microbial contaminant regions in a human reference sequence, reducing false positives during host read filtering
    """
    input:
        human_refseq = config["human_refseq"]
    output:
        human_ref_masked = "results/{genera}/1_pre_processing/mask_human_host_genome/human_genome_sequence_masked.fa"
    params:
        genera=config["genera"],
        threads=4
    log:
        stdout = "logs/{genera}/1_pre_processing/mask_human_host_genome/bbmask.out",
        stderr = "logs/{genera}/1_pre_processing/mask_human_host_genome/bbmask.err"
    shell:
        """
        module unload miniconda
        module load BBMap/38.90-GCCcore-10.2.0

        bbmask.sh \
        in={input.human_refseq} \
        out={output.human_ref_masked} \
        threads={params.threads} \
        1> {log.stdout} 2> {log.stderr}
        """

rule align_to_human_host_genome: # done
    """
    Filter reads against masked human reference genome
    """
    input:
        human_ref_masked = "results/{genera}/1_pre_processing/mask_human_host_genome/human_genome_sequence_masked.fa",
        r1 = "results/{genera}/1_pre_processing/split_unmapped_reads_beluga/{sample}/unmapped_R1.fq",
        r2 = "results/{genera}/1_pre_processing/split_unmapped_reads_beluga/{sample}/unmapped_R2.fq"
    output:
        mapped = "results/{genera}/1_pre_processing/align_to_human_host_genome/{sample}/mapped.fq",
        unmapped = "results/{genera}/1_pre_processing/align_to_human_host_genome/{sample}/unmapped.fq"
    params:
        index_path = "results/{genera}/1_pre_processing/align_to_human_host_genome/index/",
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
        stdout = "logs/{genera}/1_pre_processing/align_to_human_host_genome/{sample}/bbmap_human.out",
        stderr = "logs/{genera}/1_pre_processing/align_to_human_host_genome/{sample}/bbmap_human.err"
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

rule split_unmapped_reads_human: # done
    """
    Split interleaved reads that did not map to masked human reference into separate fastqs
    """
    input:
        unmapped = "results/{genera}/1_pre_processing/align_to_human_host_genome/{sample}/unmapped.fq"
    output:
        r1 = "results/{genera}/1_pre_processing/split_unmapped_reads_human/{sample}/unmapped_R1.fq",
        r2 = "results/{genera}/1_pre_processing/split_unmapped_reads_human/{sample}/unmapped_R2.fq"
    params:
        genera=config["genera"]
    log:
        stdout = "logs/{genera}/1_pre_processing/split_unmapped_reads_human/{sample}/reformat_human.out",
        stderr = "logs/{genera}/1_pre_processing/split_unmapped_reads_human/{sample}/reformat_human.err"
    shell:
        """
        module unload miniconda
        module load BBMap/38.90-GCCcore-10.2.0

        reformat.sh \
        in={input.unmapped} \
        out1={output.r1} \
        out2={output.r2} \
        1> {log.stdout} 2> {log.stderr}
        """

rule dedup_reads: # done
    """
    Run CD-HIT to deduplicate reads that may have been introduced during NexteraXT PCR amplification
    input:
        r1 = "results/{genera}/1_pre_processing/split_unmapped_reads_human/{sample}/unmapped_R1.fq",
        r2 = "results/{genera}/1_pre_processing/split_unmapped_reads_human/{sample}/unmapped_R2.fq"
    output:
        r1_dedup = "results/{genera}/1_pre_processing/dedup_reads/{sample}/{sample}_host_removed_dedup_R1.fastq",
        r2_dedup = "results/{genera}/1_pre_processing/dedup_reads/{sample}/{sample}_host_removed_dedup_R2.fastq"
    params:
        genera=config["genera"],
        outdir = "results/{genera}/1_pre_processing/dedup_reads/{sample}"
    log:
        stdout = "logs/{genera}/1_pre_processing/dedup_reads/{sample}/dedup_reads.out",
        stderr = "logs/{genera}/1_pre_processing/dedup_reads/{sample}/dedup_reads.err"
    shell:
        """
        module unload miniconda
        source activate /home/flg9/.conda/envs/cd-hit_ycrc
        mkdir -p {params.outdir}

        cd-hit-dup \
        -i {input.r1} -i2 {input.r2} -o {output.r1_dedup} -o2 {output.r2_dedup} \
        1> {log.stdout} 2> {log.stderr}
        """