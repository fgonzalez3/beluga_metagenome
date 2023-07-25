# Assessment of orthologous presence and absence

BUSCO will assess orthologue presence absence using blastn, a rapid method of finding close matches in large databases. It uses blastn to make sure that it does not miss any part of any possible coding sequences. To run the program, give it: 

```
#!/bin/bash
#SBATCH --job-name=busco-clostridum
#SBATCH --nodes=1
#SBATCH --mem=50G
#SBATCH --time=24:00:00
#SBATCH --output=busco.out
#SBATCH --error=busco.err

ml miniconda 
conda activate genome_annotation

busco  -i Clostridia.fasta -o Busco_Clostridia -l bacteria_odb10 -m geno
```

Once this finishes running, an output file showing annotation completeness using orthologous will appear. This *.txt file will note the total number of orthologous found, the number expected, and the number missing, giving an indication of genome completeness. 

# Genome annotation 

Prokka provides rapid annotation of prokaryotic genomes. To run this program, give it: 

```
#!/bin/bash
#SBATCH --job-name=prokka-clostridum
#SBATCH --nodes=1
#SBATCH --mem=50G
#SBATCH --time=24:00:00
#SBATCH --output=prokka.out
#SBATCH --error=prokka.err

ml miniconda 
conda activate genome_annotation

prokka --kingdom Bacteria --genus Clostridium --outdir output_directory Clostridia.fasta
```

This will provide a .gff file in the output directory that we will use to visualize genome annotations. 
