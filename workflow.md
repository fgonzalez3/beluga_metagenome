# Quality Control 

1. We received trimmed raw reads from SeqCenter, which I analyzed for this project. The first step was to run [FASTQC](https://github.com/s-andrews/FastQC) on these samples.

```
#!/bin/bash
#SBATCH --job-name=fastqc
#SBATCH --nodes=1
#SBATCH --time=24:00:00
#SBATCH --output=fastqc.out
#SBATCH --error=fastqc.err

module load miniconda 
conda activate beluga 

set -e

fastqc *.fastq

mkdir fastqc_trimmed_reads

mv *.zip fastqc_trimmed_reads
mv *.html fastqc_trimmed_reads

cd fastqc_trimmed_reads

for filename in *.zip 
do 
unzip $filename
done 

mkdir docs
cat */summary.txt > docs/fastqc_summaries.txt
```

2. Since we did not need to trim, this step was taken. Though if we were to, a similar process could be taken as the one below using [Trimmomatic](https://github.com/usadellab/Trimmomatic). 

```
#!/bin/bash
#SBATCH --job-name=trimmomatic
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --time=24:00:00
#SBATCH --output=trimmomatic.out
#SBATCH --error=trimmomatic.err

module load miniconda 
conda activate beluga 

trimmomatic PE -threads 4 whale_feces_S170_R1_001.fastq whale_feces_S170_R2_001.fastq \
whale_feces_S170_R1_001.trimmed.fastq whale_feces_S170_R1_001.orphaned.fastq \
whale_feces_S170_R2_001.trimmed.fastq whale_feces_S170_R2_001.orphaned.fastq \
ILLUMINACLIP:adapter_sequences.fasta SLIDINGWINDOW:4:20

# this is what we would run if we had adapters or needed to trim, but the .fastq files are already trimmed 
# so let's just worry about filtering out host background 
```

# Filtering

1. Since the reads were already trimmed, we could skip ahead to the filtering step. Here, first mapped paired, trimmed reads to the Beluga genome and removed all reads that mapped. This left behind only unmapped reads from microbes that we would analyze downstream. I downloaded the Beluga genome from [here](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_002288925.2/). Using [Bowtie2](https://github.com/BenLangmead/bowtie2) and [Samtools](https://github.com/samtools/samtools), I was able to get unmapped reads that would be used for assembly. These were the same steps outlined [here](https://www.metagenomics.wiki/tools/short-read/remove-host-sequences). 

```
#!/bin/bash
#SBATCH --job-name=filtering
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --mem=10G
#SBATCH --partition=week
#SBATCH --time=48:00:00
#SBATCH --output=filtering.out
#SBATCH --error=filtering.err

module load miniconda 
conda activate primer_design # conda package i used for another project that already had bowtie and samtools downloaded

# 1) bowtie2 mapping against host genome: write all (mapped and unmapped) reads to a single .bam file
# a) create read to use reference genome using beluga .fna file 

bowtie2-build beluga_genome_sequence.fna host_DB

# b) bowtie2 mapping against host sequence database, keep both aligned and unaligned reads (paired-end reads)

bowtie2 -p 8 -x host_DB -1 whale_feces_S170_R1_001.fastq -2 whale_feces_S170_R2_001.fastq -S SAMPLE_mapped_and_unmapped.sam

# c) convert file .sam to .bam

samtools view -bS SAMPLE_mapped_and_unmapped.sam > SAMPLE_mapped_and_unmapped.bam

# 2) filter required unmapped reads 
# SAMtools SAM-flag filter: get unmapped pairs (both reads R1 and R2 unmapped)

# -f extract only alignments with both reads unmapped 
# -F do not extract alignments which are not primary alignment

samtools view -b -f 12 -F 256 SAMPLE_mapped_and_unmapped.bam > SAMPLE_BOTH_reads_unmapped.bam

# 3) split paired-end reads into separated fastq files
# sort bam file by read name -n to have paired reads next to each other 

samtools sort -n -m 5G -@8 SAMPLE_BOTH_reads_unmapped.bam -o SAMPLE_BOTH_reads_unmapped.sorted.bam 
samtools fastq -@8 SAMPLE_BOTH_reads_unmapped.sorted.bam -1 SAMPLE_host_removed_R1.fastq.gz -2 SAMPLE_host_removed_R2.fastq.gz -0 /dev/null -s /dev/null -n
```

# Metagenome assembly 

1. Once reads underwent filtering, I assembled a metagenome. [SPAdes](https://github.com/ablab/spades) is best for this data. 

```
#!/bin/bash
#SBATCH --job-name=assembly
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --mem=10G
#SBATCH --partition=week
#SBATCH --time=48:00:00
#SBATCH --output=assembly.out
#SBATCH --error=assembly.err

module load SPAdes/3.15.1-GCCcore-10.2.0-Python-3.8.6

mkdir results

metaspades.py -1 SAMPLE_host_removed_R1.fastq.gz -2 SAMPLE_host_removed_R2.fastq.gz -o results
```

# Read mapping 

1. Here, map reads to contigs to get coverage data. I converted my gunzip files to fasta files using [Seqtk](https://github.com/lh3/seqtk). 

```
#!/bin/bash
#SBATCH --job-name=alignment
#SBATCH --nodes=2
#SBATCH --time=24:00:00
#SBATCH --output=alignment.out
#SBATCH --error=alignment.err

ml miniconda 
conda activate primer_design

# bowtie2 mapping host_removed reads to contigs for binning input 

seqtk seq -a SAMPLE_host_removed_R1.fastq.gz > out1.fa
seqtk seq -a SAMPLE_host_removed_R2.fastq.gz > out2.fa

bowtie2-build -f results/contigs.fasta index_prefix
bowtie2 -x index_prefix -f -1 out1.fa -2 out2.fa -S contig_alignments.bam

samtools sort contig_alignments.bam -o sample.sorted.bam
samtools index sample.sorted.bam -o indexed.bam
```

# Binning 

1. Next, bin your MAGs. I used [MaxBin2](https://anaconda.org/bioconda/maxbin2), but there are others you can use as well. I first tried this on the HCC, a little simpler this way. 

```
#!/bin/bash
#SBATCH --job-name=binning
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --time=24:00:00
#SBATCH --output=binning.out
#SBATCH --error=binning.err

module load MaxBin/2.2.7-gompi-2020b 

run_MaxBin.pl -thread 8 -contig results/contigs.fasta -reads SAMPLE_host_removed_R1.fastq.gz -reads2 SAMPLE_host_removed_R2.fastq.gz -out MAXBIN
```
1a. Conda didn't automatically point to the Perl libraries for these functions so we initially tried to include a few lines to ensure that my Conda environment was looking in the right place for those modules. This didn't end up working out, which is why I had McCleary install MaxBin instead. Either way, this is still useful if this issue arises again in the future. 

```
#!/bin/bash
#SBATCH --job-name=binning
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --time=24:00:00
#SBATCH --output=binning.out
#SBATCH --error=binning.err

ml miniconda 
conda activate binning 

mkdir -p $CONDA_PREFIX/etc/conda/activate.d

echo 'export PERL5LIB='/gpfs/gibbs/project/turner/flg9/conda_envs/binning/lib/site_perl/5.26.2/:/gpfs/gibbs/project/turner/flg9/conda_envs/binning/lib/5.26.2 >> $CONDA_PREFIX/etc/conda/activate.d/env_vars.sh

conda deactivate binning
conda activate binning

run_MaxBin.pl -thread 8 -contig results/contigs.fasta -reads SAMPLE_host_removed_R1.fastq.gz -reads2 SAMPLE_host_removed_R2.fastq.gz -out MAXBIN
```

2. Then, assess the quality of bins and MAGs. I used [CHECKM](https://github.com/Ecogenomics/CheckM). 

```
#!/bin/bash
#SBATCH --job-name=checkm
#SBATCH --nodes=1
#SBATCH --time=12:00:00
#SBATCH --ntasks-per-node=4
#SBATCH --mem=60G
#SBATCH --output=checkm.out
#SBATCH --error=checkm.err

ml miniconda 
conda activate qc_binning 

# checkm calls genes internally using prodigal 

mkdir binning/CHECKM
checkm lineage_wf -t 4 -x fasta binning binning/CHECKM

# run qa to make tables and plots 

checkm qa binning/CHECKM/lineage.ms binning/CHECKM/ --file binning/CHECKM/quality.tsv --tab_table -o 2
```

# Taxonomy

1. Next, it's time to assign taxonomy. There are two ways to do this a) assign to paired reads or b) assign to contigs. I did both, but we first have to create a database using [Kraken2](https://github.com/DerrickWood/kraken2) if we are running on Yale's cluster. HCC already has this downloaded, so I only needed to preload it. 

```
#!/bin/bash
#SBATCH --job-name=kraken_DB
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=24:00:00
#SBATCH --output=kraken_DB.out
#SBATCH --error=kraken_DB.err

module load kraken2/2.0.8-beta	 

kraken2-build --standard --threads 24 --db KRAKEN2_DB
```

2. Now assign taxonomy to contigs. I first tried this on the HCC since they already have Kraken2 installed - much more straightforward than running it through Yale's cluster. 

```
#!/bin/bash
#SBATCH --job-name=kraken
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=12
#SBATCH --mem=60G
#SBATCH --time=48:00:00
#SBATCH --output=kraken.out
#SBATCH --error=kraken.err

module load kraken2

kraken2 --preload --db $KRAKEN2_DB

mkdir TAXONOMY_MAG

# run kraken2 with the full input data and monitor memory usage
/usr/bin/time -v kraken2 --db $KRAKEN2_DB --threads 12 --memory-mapping --output TAXONOMY_MAG/contigs.kraken --report TAXONOMY_MAG/contigs.report contigs.fasta 2> memory_usage.txt 
```

3. Then assign taxonomy to reads. 

```
#!/bin/bash
#SBATCH --job-name=kraken
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=12
#SBATCH --mem=60G
#SBATCH --time=48:00:00
#SBATCH --output=kraken.out
#SBATCH --error=kraken.err

module load kraken2

kraken2 --preload --db $KRAKEN2_DB

mkdir TAXONOMY_MAG

# run kraken2 with the full input data and monitor memory usage
/usr/bin/time -v kraken2 --db $KRAKEN2_DB --threads 12 --memory-mapping --output TAXONOMY_MAG/contigs.kraken --report TAXONOMY_MAG/contigs.report --paired SAMPLE_host_removed_R1.fastq.gz SAMPLE_host_removed_R2.fastq.gz 2> memory_usage.txt
```

4. To run on McCleary, a few extra steps (similar to getting MaxBin to run) had to be taken.

```
ml miniconda
conda activate taxonomy

cd project/conda_envs/taxonomy/

####finds directory that has FTP.pm inside based off of error provided####
find . -type d -name 'Net'

#####enter directory where Net/FTP.pm is located##########
cd lib/perl5/core_perl/

####get full path location to library#########
readlink -f .

####creates file that is read everytime conda environment is activated#########
mkdir -p $CONDA_PREFIX/etc/conda/activate.d

#####saves library location from readlink -f . to file in activate.d that will tell conda where perl5 library is everytime###########
echo 'export PERL5LIB=/gpfs/gibbs/project/turner/flg9/conda_envs/taxonomy/lib/perl5/core_perl' >> $CONDA_PREFIX/etc/conda/activate.d/env_vars.sh

conda deactivate
```

Now, my Conda environment should have the correct location for the Kraken2 libraries meaning we can run the workflow in the same way we did for the HCC. Building databases take up quite a bit of RAM, so it's necessary to load in ~150-200G of memory to avoid crashes. 

```
#!/bin/bash
#SBATCH --job-name=kraken
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --mem=200G
#SBATCH --time=24:00:00
#SBATCH --output=kraken.out
#SBATCH --error=kraken.err

module load miniconda 
conda activate taxonomy

# first build a database using the standard db

kraken2-build --standard --threads 4 --db KRAKEN2_DB

# then run on contigs 

kraken2 --db KRAKEN2_DB --threads 4 --output TAXONOMY_MAG/contigs.kraken --report TAXONOMY_MAG/contigs.report contigs.fasta 
```

# Visualization

1. Now let's visualize these outputs. I used [Krona](https://github.com/marbl/Krona/wiki) for this. 

```
#!/bin/bash
#SBATCH --job-name=krona
#SBATCH --nodes=1
#SBATCH --mem=8G
#SBATCH --time=6:00:00
#SBATCH --output=krona.out
#SBATCH --error=krona.err

module load krona 

cut -f2,3 contigs.kraken > krona.input
ktImportTaxonomy krona.input -o krona.out.html
```
This looks great, now let's figure out how to import this into R to manipulate it some more. 
![snapshot](https://github.com/fgonzalez3/beluga_metagenome/assets/51669806/b05ada59-bdc1-4671-a81e-dc983f3f99e4)

# Abundance 

To get abundance metrics useful for R visualization, we first run [BRACKEN](https://github.com/jenniferlu717/Bracken) on our KRAKEN report file. I created a new Conda environment using Bioconda, link [here](https://anaconda.org/bioconda/bracken). This will require quite a bit of CPU and RAM to get it running properly. Also, multiple threads or it will take multiple days with only one.  

```
#!/bin/bash
#SBATCH --job-name=bracken
#SBATCH --nodes=2
#SBATCH --mem=200G
#SBATCH --ntasks=10
#SBATCH --cpus-per-task=10
#SBATCH --time=24:00:00
#SBATCH --output=bracken.out
#SBATCH --error=bracken.err

module load miniconda 
conda activate abundance

mkdir BRACKEN

# first build a BRACKEN db

bracken-build -d /gpfs/gibbs/project/turner/flg9/TurnerLab/beluga_feces/taxonomy/KRAKEN2_DB -t 15

# then run BRACKEN for abundance estimation 

bracken -d KRAKEN2_DB -i contigs.report -o BRACKEN/contigs.bracken -t 15
```
