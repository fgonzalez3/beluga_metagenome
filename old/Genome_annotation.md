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

This will provide a .gff file in the output directory that we will use to visualize genome annotations within the Integrative Genomics Viewer (IGV). 

# Anvi'o 

To annotate gene function for downstream R analysis, I used a combination of the Anvi'o [metagenomics workshop pipeline](https://github.com/fgonzalez3/beluga_metagenome/blob/main/README.pdf) and their [metagenomics workflow](https://merenlab.org/2016/06/22/anvio-tutorial-v2/). Below are some rough scripts that I used for that. 

```
anvi-interactive -d Clostridia_module.txt \
                 -p 03_PROFILE/anvio \
                 --manual-mode \
                 --dry-run



bowtie2-build FMT_HIG_FITNESS_KC_MAG_00120.fasta FMT_HIG_FITNESS_KC_MAG_00120
 
# Make a directory to house the mapping results
  mkdir -p 02_MAPPING
  # use a for loop to map the recipient gut metagenomes from PRE and POST FMT metagenom
  # against our MAG reference
  for FASTA in anvio; do
      # 1. perform read recruitment with bowtie2 to get a SAM file:
      echo -e "Mapping: "${FASTA}""
      bowtie2 --threads 4 \
              -x Clostridia \
              -1 whale_feces_S170_R1_001.fastq \
              -2 whale_feces_S170_R2_001.fastq \
              --no-unal \
              -S 02_MAPPING/"${FASTA}".sam
      # 2. covert the resulting SAM file to a BAM file:
      samtools view -F 4 -bS  02_MAPPING/"${FASTA}".sam > 02_MAPPING/"${FASTA}
      # 3. sort the BAM file:
      samtools sort 02_MAPPING/"${FASTA}"-RAW.bam -o 02_MAPPING/"${FASTA}".bam
      # 4. index the BAM file:
      samtools index 02_MAPPING/"${FASTA}".bam
  done


mkdir -p 02_MAPPING

for FASTA in Beluga do; 
echo -e "Mapping: "${FASTA}""
bowtie2 --threads 4 \
              -x Clostridia \
              -1 whale_feces_S170_R1_001.fastq \
              -2 whale_feces_S170_R2_001.fastq \
              --no-unal \
              -S 02_MAPPING/"${FASTA}".sam

      samtools view -F 4 -bS  02_MAPPING/"${FASTA}".sam > 02_MAPPING/"${FASTA}
      samtools sort 02_MAPPING/"${FASTA}"-RAW.bam -o 02_MAPPING/"${FASTA}".bam
      samtools index 02_MAPPING/"${FASTA}".bam

done

bowtie2 --threads 4 -x Clostridia -1 whale_feces_S170_R1_001.fastq -2 whale_feces_S170_R2_001.fastq --no-unal -S 02_MAPPING/alignment.sam

samtools view -F 4 -bS  02_MAPPING/alignment.sam > 02_MAPPING/alignment-RAW.bam
samtools sort 02_MAPPING/alignment-RAW.bam -o 02_MAPPING/alignment.bam
samtools index 02_MAPPING/alignment.bam

anvi-profile -c Clostridia.db -i 02_MAPPING/alignment.bam --num-threads 4 -o 03_PROFILE --cluster-contigs

anvi-interactive -p 03_PROFILE/PROFILE.db -c Clostridia.db

  anvi-estimate-metabolism -e external-genomes.txt -O Clostridia


  anvi-estimate-metabolism -e external-genomes.txt -O Clostridia --matrix-format
    anvi-matrix-to-newick Clostridia-completeness-MATRIX.txt

    anvi-interactive -A Clostridia-completeness-MATRIX.txt -t Clostridia-completeness-MATRIX.txt.newick -p Clostridia-completeness-MATRIX_profile.db --manual-mode


bowtie2 --threads 4 -x 04_MAPPING/contigs -1 01_QC/Sample_01-QUALITY_PASSED_R1.fastq -2 01_QC/Sample_01-QUALITY_PASSED_R2.fastq -S 04_MAPPING/Sample_01.sam

scp *fastq.gz transfer-mccleary.ycrc.yale.edu:/gpfs/gibbs/project/turner/flg9/TurnerLab/beluga_feces/china_metagenome


 
anvi-estimate-metabolism -c FMT_HIG_FITNESS_KC_MAG_00120.db -O KC_MAG_00120.txt
```

