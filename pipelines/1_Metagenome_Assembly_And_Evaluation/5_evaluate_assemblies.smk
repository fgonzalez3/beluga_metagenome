
rule filter_indivudal_assemblies: # test
    """
    Filter for contigs that are a minimum 1.5kb in length. 
    This is the standard length for binning, so we'll evaluate them at this length prior to binning.
    """
    input:
        # SPAdes assemblies 
        c1 = "results/{genera}/3_dedup_contigs/SPAdes_single/{sample}/{sample}_DEDUP95.fasta",
        # Megahit assemblies
        c2 = "results/{genera}/3_dedup_contigs/megahit_single/{sample}/{sample}_DEDUP95.fasta"
    output:
        # List of filtered contigs 
        SPAdes_filter = "results/{genera}/5_evaluate_assemblies/filter_individual_assemblies/{sample}/metaspades_assembly_DEDUP95_m1500.fasta",
        Megahit_filter = "results/{genera}/5_evaluate_assemblies/filter_individual_assemblies/{sample}/megahit_assembly_DEDUP95_m1500.fasta"
    params:
        len = 1500,
        threads = 4,
        outdir = "results/{genera}/5_evaluate_assemblies/filter_individual_assemblies/{sample}",
    log:
        stdout = "logs/{genera}/5_evaluate_assemblies/filter_individual_assemblies/{sample}/filter_assemblies.out",
        stderr = "logs/{genera}/5_evaluate_assemblies/filter_individual_assemblies/{sample}/filter_assemblies.err"
    shell:
        """
        module unload miniconda
        module load SeqKit/2.8.1

        # 1. Create outdir
        mkdir -p {params.outdir}

        set -x
        echo "Running seqkit..." 1>> {log.stdout} 2>> {log.stderr}

        # 2. Filter out contigs <1.5kb
        seqkit seq \
        {input.c1} --min-len {params.len} --threads {params.threads} -o {output.SPAdes_filter} 

        seqkit seq \
        {input.c2} --min-len {params.len} --threads {params.threads} -o {output.Megahit_filter} 

        echo "Running stats.sh..." 1>> {log.stdout} 2>> {log.stderr}
        """

rule metaquast_individual_assemblies: # test
    """
    Run MetaQUAST on all filtered assemblies from all samples for a global comparison.
    """
    input:
        c1 = expand("results/{genera}/5_evaluate_assemblies/filter_individual_assemblies/{sample}/metaspades_assembly_DEDUP95_m1500.fasta", sample=SAMPLES, genera=config["genera"]),
        c2 = expand("results/{genera}/5_evaluate_assemblies/filter_individual_assemblies/{sample}/megahit_assembly_DEDUP95_m1500.fasta", sample=SAMPLES, genera=config["genera"])
    output:
        whole_assembly_stats = "results/{genera}/5_evaluate_assemblies/individual_assembly_eval/assembly_stats.csv",
        report = "results/{genera}/5_evaluate_assemblies/individual_assembly_eval/report.html"
    params:
        outdir = "results/{genera}/5_evaluate_assemblies/indvidual_assembly_eval",
        labels = "SPAdes.individual.m1500,Megahit.individual.m1500",
        threads = 4
    log:
        stdout = "logs/{genera}/5_evaluate_assemblies/individual_assembly_eval/assembly_eval.out",
        stderr = "logs/{genera}/5_evaluate_assemblies/individual_assembly_eval/assembly_eval.err"
    shell:
        """
        module unload miniconda
        module load BBMap/38.90-GCCcore-10.2.0
        source activate /home/flg9/.conda/envs/quast

        mkdir -p {params.outdir}

        set -x

        #1. Generate specific assembly data covering quality of contigs generated

        echo "Running stats.sh..." 1>> {log.stdout} 2>> {log.stderr}

        #stats.sh \
        #in={output.filtered_contigs} out={output.whole_assembly_stats}

        statswrapper.sh \
        in={input.c1},{input.c2} out={output.whole_assembly_stats}

        # 2. Generate comprehensive evaluation of assembly data using metaquast

        echo "Running metaquast.py..." 1>> {log.stdout} 2>> {log.stderr}

        metaquast.py \
        -t {params.threads} --labels {params.outfile_labels} \
        --output-dir {params.outdir} {output.report} \
        1>> {log.stdout} 2>> {log.stderr}
        """

rule filter_master_coassemblies:
    """
    Filter for contigs that are a minimum 1.5kb in length. 
    This is the standard length for binning, so we'll evaluate them at this length prior to binning.
    """
    input:
        # SPAdes assemblies 
        c1 = "results/{genera}/3_dedup_contigs/dedup_contigs/SPAdes_master/DEDUP95.fasta",
        c2 = "results/{genera}/3_dedup_contigs/dedup_contigs/megahit_master/DEDUP95.fasta"
    output:
        # List of filtered contigs 
        SPAdes_filter = "results/{genera}/5_evaluate_assemblies/filter_master_assemblies/metaspades_assembly_DEDUP95_m1500.fasta",
        Megahit_filter = "results/{genera}/5_evaluate_assemblies/filter_master_assemblies/megahit_assembly_DEDUP95_m1500.fasta"
    params:
        len = 1500,
        threads = 4,
        outdir = "results/{genera}/5_evaluate_assemblies/filter_master_assemblies",
    log:
        stdout = "logs/{genera}/5_evaluate_assemblies/filter_master_assemblies/filter_assemblies.out",
        stderr = "logs/{genera}/5_evaluate_assemblies/filter_master_assemblies/filter_assemblies.err"
    shell:
        """
        module unload miniconda
        module load SeqKit/2.8.1

        # 1. Create outdir
        mkdir -p {params.outdir}

        set -x
        echo "Running seqkit..." 1>> {log.stdout} 2>> {log.stderr}

        # 2. Filter out contigs <1.5kb
        seqkit seq \
        {input.c1} --min-len {params.len} --threads {params.threads} -o {output.SPAdes_filter} 

        seqkit seq \
        {input.c2} --min-len {params.len} --threads {params.threads} -o {output.Megahit_filter} 

        echo "Running stats.sh..." 1>> {log.stdout} 2>> {log.stderr}
        """

rule metaquast_master_coassemblies:
    """
    Run MetaQUAST on all filtered assemblies from all samples for a global comparison.
    """
    input:
        c1 = expand("results/{genera}/5_evaluate_assemblies/filter_master_assemblies/metaspades_assembly_DEDUP95_m1500.fasta", sample=SAMPLES, genera=config["genera"]),
        c2 = expand("results/{genera}/5_evaluate_assemblies/filter_master_assemblies/megahit_assembly_DEDUP95_m1500.fasta", sample=SAMPLES, genera=config["genera"])
    output:
        whole_assembly_stats = "results/{genera}/5_evaluate_assemblies/master_assembly_eval/assembly_stats.csv",
        report = "results/{genera}/5_evaluate_assemblies/master_assembly_eval/report.html"
    params:
        outdir = "results/{genera}/5_evaluate_assemblies/master_assembly_eval",
        labels = "SPAdes.master.m1500,Megahit.master.m1500",
        threads = 4
    log:
        stdout = "logs/{genera}/5_evaluate_assemblies/master_assembly_eval/assembly_eval.out",
        stderr = "logs/{genera}/5_evaluate_assemblies/master_assembly_eval/assembly_eval.err"
    shell:
        """
        module unload miniconda
        module load BBMap/38.90-GCCcore-10.2.0
        source activate /home/flg9/.conda/envs/quast

        mkdir -p {params.outdir}

        set -x

        #1. Generate specific assembly data covering quality of contigs generated

        echo "Running stats.sh..." 1>> {log.stdout} 2>> {log.stderr}

        #stats.sh \
        #in={output.filtered_contigs} out={output.whole_assembly_stats}

        statswrapper.sh \
        in={input.c1},{input.c2} out={output.whole_assembly_stats}

        # 2. Generate comprehensive evaluation of assembly data using metaquast

        echo "Running metaquast.py..." 1>> {log.stdout} 2>> {log.stderr}

        metaquast.py \
        -t {params.threads} --labels {params.outfile_labels} \
        --output-dir {params.outdir} {output.report} \
        1>> {log.stdout} 2>> {log.stderr}
        """

rule filter_whale_coassemblies:
    """
    Filter for contigs that are a minimum 1.5kb in length. 
    This is the standard length for binning, so we'll evaluate them at this length prior to binning.
    """
    input:
        # SPAdes assemblies 
        c1 = "results/{genera}/3_dedup_contigs/dedup_contigs/SPAdes_master/DEDUP95.fasta",
        c2 = "results/{genera}/3_dedup_contigs/dedup_contigs/megahit_master/DEDUP95.fasta"
    output:
        # List of filtered contigs 
        SPAdes_filter = "results/{genera}/5_evaluate_assemblies/filter_master_assemblies/metaspades_assembly_DEDUP95_m1500.fasta",
        Megahit_filter = "results/{genera}/5_evaluate_assemblies/filter_master_assemblies/megahit_assembly_DEDUP95_m1500.fasta"
    params:
        len = 1500,
        threads = 4,
        outdir = "results/{genera}/5_evaluate_assemblies/filter_master_assemblies",
    log:
        stdout = "logs/{genera}/5_evaluate_assemblies/filter_master_assemblies/filter_assemblies.out",
        stderr = "logs/{genera}/5_evaluate_assemblies/filter_master_assemblies/filter_assemblies.err"
    shell:
        """
        module unload miniconda
        module load SeqKit/2.8.1

        # 1. Create outdir
        mkdir -p {params.outdir}

        set -x
        echo "Running seqkit..." 1>> {log.stdout} 2>> {log.stderr}

        # 2. Filter out contigs <1.5kb
        seqkit seq \
        {input.c1} --min-len {params.len} --threads {params.threads} -o {output.SPAdes_filter} 

        seqkit seq \
        {input.c2} --min-len {params.len} --threads {params.threads} -o {output.Megahit_filter} 

        echo "Running stats.sh..." 1>> {log.stdout} 2>> {log.stderr}
        """

rule metaquast_whale_coassemblies:
    """
    Run MetaQUAST on all filtered assemblies from all samples for a global comparison.
    """
    input:
        c1 = expand("results/{genera}/5_evaluate_assemblies/filter_master_assemblies/metaspades_assembly_DEDUP95_m1500.fasta", sample=SAMPLES, genera=config["genera"]),
        c2 = expand("results/{genera}/5_evaluate_assemblies/filter_master_assemblies/megahit_assembly_DEDUP95_m1500.fasta", sample=SAMPLES, genera=config["genera"])
    output:
        whole_assembly_stats = "results/{genera}/5_evaluate_assemblies/master_assembly_eval/assembly_stats.csv",
        report = "results/{genera}/5_evaluate_assemblies/master_assembly_eval/report.html"
    params:
        outdir = "results/{genera}/5_evaluate_assemblies/master_assembly_eval",
        labels = "SPAdes.master.m1500,Megahit.master.m1500",
        threads = 4
    log:
        stdout = "logs/{genera}/5_evaluate_assemblies/master_assembly_eval/assembly_eval.out",
        stderr = "logs/{genera}/5_evaluate_assemblies/master_assembly_eval/assembly_eval.err"
    shell:
        """
        module unload miniconda
        module load BBMap/38.90-GCCcore-10.2.0
        source activate /home/flg9/.conda/envs/quast

        mkdir -p {params.outdir}

        set -x

        #1. Generate specific assembly data covering quality of contigs generated

        echo "Running stats.sh..." 1>> {log.stdout} 2>> {log.stderr}

        #stats.sh \
        #in={output.filtered_contigs} out={output.whole_assembly_stats}

        statswrapper.sh \
        in={input.c1},{input.c2} out={output.whole_assembly_stats}

        # 2. Generate comprehensive evaluation of assembly data using metaquast

        echo "Running metaquast.py..." 1>> {log.stdout} 2>> {log.stderr}

        metaquast.py \
        -t {params.threads} --labels {params.outfile_labels} \
        --output-dir {params.outdir} {output.report} \
        1>> {log.stdout} 2>> {log.stderr}
        """