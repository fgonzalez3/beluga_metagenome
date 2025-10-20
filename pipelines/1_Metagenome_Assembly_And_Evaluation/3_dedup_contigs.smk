# This pipeline goes through the following steps - 
    # 1-3. De-duplication of contigs assembled from each of our assemblies
       # 1. Individual assemblies from each PE read file from both SPAdes and Megahit assemblies
       # 2. Co-assemblies produced from concatenating PE read files for each corresponding whale
       # 3. Co-assemblies produced from concatenating all PE read file regardless of individual

rule deduplicate_contigs_individual_assemblies:
    """
    Deduplicate overrepresented contigs that may have been introduced during metagenome assembly
    """
    input:
        lambda wildcards: (
            f"results/{wildcards.genera}/testing/2_assembly/{wildcards.assembler}/individual_metagenome_assembly/{wildcards.sample}/"
            + ("contigs.fasta" if wildcards.assembler == "spades" else "final.contigs.fa")
        )
    output:
        "results/{genera}/testing/3_dedup_contigs/{assembler}/individual_metagenome_assembly/{sample}/{sample}_DEDUP95.fasta"
    params:
        outdir = "results/{genera}/testing/3_dedup_contigs/{assembler}/individual_metagenome_assembly/{sample}"
    log:
        stdout = "logs/{genera}/testing/3_dedup_contigs/{assembler}/individual_metagenome_assembly/{sample}/dedup_contigs.out",
        stderr = "logs/{genera}/testing/3_dedup_contigs/{assembler}/individual_metagenome_assembly/{sample}/dedup_contigs.err"
    shell:
        """
        mkdir -p {params.outdir}

        module unload miniconda
        module load CD-HIT/4.8.1-GCC-10.2.0

        cd-hit-est \
            -i {input} \
            -o {output} \
            -T 4 -c 0.95 \
            1>> {log.stdout} 2>> {log.stderr}
        """