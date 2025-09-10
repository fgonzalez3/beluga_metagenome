INDIVIDUALS = ["JUNO", "KELA", "NATTY"]

def indiv_r1(wildcards):
    return sorted(glob.glob(f"results/{wildcards.genera}/1_pre_processing/dedup_reads/{sample}/Sample_{wildcards.individual}_*_host_removed_dedup_R1.fastq"))

def indiv_r2(wildcards):
    return sorted(glob.glob(f"results/{wildcards.genera}/1_pre_processing/dedup_reads/{sample}/Sample_{wildcards.individual}_*_host_removed_dedup_R2.fastq"))

rule master_metagenome_coassembly_spades: # test
    """
    Generate a master co-assembly from merged PE files of all samples using SPAdes
    """
    input:
        r1 = "results/{genera}/1_pre_processing/master_coassembly_normalization/all_samples_ME_R1_norm.fq.gz",
        r2 = "results/{genera}/1_pre_processing/master_coassembly_normalization/all_samples_ME_R2_norm.fq.gz"
    output:
        "results/{genera}/2_assembly/SPAdes/master_metagenome_coassembly/contigs.fasta"
    params:
        genera=config["genera"],
        outdir = "results/{genera}/2_assembly/SPAdes/master_metagenome_coassembly/"
    log:
        stdout = "logs/{genera}/2_assembly/SPAdes/master_metagenome_coassembly/assembly.out",
        stderr = "logs/{genera}/2_assembly/SPAdes/master_metagenome_coassembly/assembly.err"
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
        "results/{genera}/2_assembly/megahit/master_metagenome_coassembly/final.contigs.fa"
    params:
        genera=config["genera"],
        preset="meta-large"
        outdir = "results/{genera}/2_assembly/megahit/master_metagenome_coassembly"
    log:
        stdout = "logs/{genera}/2_assembly/megahit/master_metagenome_coassembly/assembly.out",
        stderr = "logs/{genera}/2_assembly/megahit/master_metagenome_coassembly/assembly.err"
    shell:
        """
        module unload miniconda
        source activate /home/flg9/.conda/envs/megahit

        megahit \
        -1 {input.r1} -2 {input.r2} -o {params.outdir} \
        --presets {params.preset} \
        1> {log.stdout} 2> {log.stderr}
        """

rule individual_metagenome_coassembly_spades: # test
    """
    Generate a co-assembly using spades per individual
    """
    input:
        r1 = "results/{genera}/1_pre_processing/individual_coassembly_normalization/{individual}_ME_R1_norm.fq.gz",
        r2 = "results/{genera}/1_pre_processing/individual_coassembly_normalization/{individual}_ME_R2_norm.fq.gz"
    output:
        "results/{genera}/2_assembly/SPAdes/individual_metagenome_coassembly/{individual}/contigs.fasta"
    params:
        genera=config["genera"],
        outdir = "results/{genera}/2_assembly/SPAdes/individual_metagenome_coassembly/{individual}"
    log:
        stdout = "logs/{genera}/2_assembly/SPAdes/individual_metagenome_coassembly/{individual}/assembly.out",
        stderr = "logs/{genera}/2_assembly/SPAdes/individual_metagenome_coassembly/{individual}/assembly.err"
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
        "results/{genera}/2_assembly/megahit/indvidual_metagenome_coassembly/{individual}/final.contigs.fa"
    params:
        genera=config["genera"],
        preset="meta-large"
        outdir = "results/{genera}/2_assembly/megahit/indvidual_metagenome_coassembly/{individual}"
    log:
        stdout = "logs/{genera}/2_assembly/megahit/indvidual_metagenome_coassembly/{individual}/assembly.out",
        stderr = "logs/{genera}/2_assembly/megahit/indvidual_metagenome_coassembly/{individual}/assembly.err"
    shell:
        """
        module unload miniconda
        source activate /home/flg9/.conda/envs/megahit

        megahit \
        -1 {input.r1} -2 {input.r2} -o {params.outdir} \
        --presets {params.preset} \
        1> {log.stdout} 2> {log.stderr}
        """