# Abundance 

I next wanted to take a look at the Krona output closer. To do this, I had to first call upon Bracken to create a new report that contained abundance metrics for various microbe levels. 

```
#!/bin/bash
#SBATCH --job-name=bracken
#SBATCH --nodes=1
#SBATCH --mem=50G
#SBATCH --time=1:00:00
#SBATCH --output=bracken.out
#SBATCH --error=bracken.err

module load miniconda
conda activate abundance

# Run Bracken for each taxonomic level individually
bracken -d KRAKEN2_DB -i contigs.report -o bracken_level_1.txt -l P -t 10
bracken -d KRAKEN2_DB -i contigs.report -o bracken_level_2.txt -l C -t 10
bracken -d KRAKEN2_DB -i contigs.report -o bracken_level_3.txt -l O -t 10
bracken -d KRAKEN2_DB -i contigs.report -o bracken_level_4.txt -l F -t 10
bracken -d KRAKEN2_DB -i contigs.report -o bracken_level_5.txt -l G -t 10
bracken -d KRAKEN2_DB -i contigs.report -o bracken_level_6.txt -l S -t 10

# Concatenate the results into a single output file
cat bracken_level_*.txt > bracken.txt
```

Next, I uploaded the Bracken output to R and visualized it. My first visualization is seen below. You can file the markdown file for that [here](). 

<img width="417" alt="Screen Shot 2023-07-21 at 3 47 28 PM" src="https://github.com/fgonzalez3/beluga_metagenome/assets/51669806/4cf9f5dd-905b-4f6f-86a9-0a0bbe244641">
