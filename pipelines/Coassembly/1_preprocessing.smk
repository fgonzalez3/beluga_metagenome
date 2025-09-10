
rule master_coassembly_normalization: # test
    """
    Concatenate all forward and reverse PE files into two separate files for co-assembly of all samples.
    Further, normalize these merged PE files to reduce redundancy.
    """
    input:
        r1 = "results/{genera}/1_pre_processing/dedup_reads/{sample}/{sample}_host_removed_dedup_R1.fastq",
        r2 = "results/{genera}/1_pre_processing/dedup_reads/{sample}/{sample}_host_removed_dedup_R2.fastq"
    output:
        r1 = "results/{genera}/1_pre_processing/master_coassembly_normalization/all_samples_R1.fq",
        r2 = "results/{genera}/1_pre_processing/master_coassembly_normalization/all_samples_R2.fq",
        norm1 = "results/{genera}/1_pre_processing/master_coassembly_normalization/all_samples_ME_R1_norm.fq.gz",
        norm2 = "results/{genera}/1_pre_processing/master_coassembly_normalization/all_samples_ME_R2_norm.fq.gz"
    params:
        genera=config["genera"],
        target=70,
        mindepth=2,
        threads=4,
        prefilter="t"
    log:
        stdout = "logs/{genera}/1_pre_processing/master_coassembly_normalization/{sample}/dedup_reads.out",
        stderr = "logs/{genera}/1_pre_processing/master_coassembly_normalization/{sample}/dedup_reads.err"
    shell:
        """
        module unload miniconda 
        module load BBMap/38.90-GCCcore-10.2.0

        # 1. Merge all raw reads into one master co-assembly
        cat {input.fwd} >> {output.r1}
        cat {input.rev} >> {output.r2}

        # 2. Normalize master co-assembly
        bbnorm.sh \
        in={output.r1} in2={output.r2} \
        out={output.norm1} out2={output.norm2} \
        target={params.target} mindepth={params.mindepth} \
        threads={params.threads} prefilter={params.prefilter} \
        1>> {log.stdout} 2>> {log.stderr}
        """

# Define individuals for this next rule 

INDIVIDUALS = ["JUNO", "KELA", "NATTY"]

def indiv_r1(wildcards):
    return sorted(glob.glob(f"results/{wildcards.genera}/1_pre_processing/dedup_reads/{sample}/Sample_{wildcards.individual}_*_host_removed_dedup_R1.fastq"))

def indiv_r2(wildcards):
    return sorted(glob.glob(f"results/{wildcards.genera}/1_pre_processing/dedup_reads/{sample}/Sample_{wildcards.individual}_*_host_removed_dedup_R2.fastq"))

rule individual_coassembly_normalization: # test
    """
    Concatenate all forward and reverse PE files from each individual for co-assembly, then normalize
    """
    input:
        r1 = indiv_r1,
        r2 = indiv_r2
    output:
        r1 = "results/{genera}/1_pre_processing/individual_coassembly_normalization/{individual}_ME_R1.fastq.gz",
        r2 = "results/{genera}/1_pre_processing/individual_coassembly_normalization/{individual}_ME_R2.fastq.gz",
        norm1 = "results/{genera}/1_pre_processing/individual_coassembly_normalization/{individual}_ME_R1_norm.fq.gz",
        norm2 = "results/{genera}/1_pre_processing/individual_coassembly_normalization/{individual}_ME_R2_norm.fq.gz"
    params:
        target=70,
        mindepth=2,
        threads=4,
        prefilter="t"
    log:
        stdout = "logs/{genera}/1_pre_processing/individual_coassembly_normalization/{individual}/dedup_reads.out",
        stderr = "logs/{genera}/1_pre_processing/individual_coassembly_normalization/{individual}/dedup_reads.err"
    shell:
        """
        module unload miniconda 
        module load BBMap/38.90-GCCcore-10.2.0

        # 1. Merge by individual 
        cat {input.r1} > {output.r1}
        cat {input.r2} > {output.r2}

        # 2. Normalize by individual
        bbnorm.sh \
        in={output.r1} in2={output.r2} \
        out={output.norm1} out2={output.norm2} \
        target={params.target} mindepth={params.mindepth} \
        threads={params.threads} prefilter={params.prefilter} \
        1>> {log.stdout} 2>> {log.stderr}
        """