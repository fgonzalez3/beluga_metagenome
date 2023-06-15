### Quality Control 

1. We received trimmed raw reads from SeqCenter, which I analyzed for this project. The first step was to run FASTQC on these samples:

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

### Filtering

1. Since the reads were already trimmed, we could skip ahead to the filtering step. Here, first mapped paired, trimmed reads to the Beluga genome and removed all reads that mapped. This left behind only unmapped reads from microbes that we would analyze downstream. I downloaded the Beluga genome from [here](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_002288925.2/). Using Bowtie and Samtools, I was able to get unmapped reads that would be used for assembly. 

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

### Metagenome assembly 

1. Once reads underwent filtering, I assembled a metagenome. SPAdes is best for this data. 

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

### Binning 

1. Next I binned my MAGs. There are two programs you could use here. 

```
#!/bin/bash
#SBATCH --job-name=binning
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --time=24:00:00
#SBATCH --output=binning.out
#SBATCH --error=binning.err

module load miniconda 
conda activate binning 

mkdir MAXBIN

run_MaxBin.pl -thread 8 -contig contigs.fasta -reads SAMPLE_host_removed_R1.fastq.gz -reads2 SAMPLE_host_removed_R2.fastq.gz -o MAXBIN
```

OR

```
#!/bin/bash
#SBATCH --job-name=metabat
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --time=24:00:00
#SBATCH --output=metabat.out
#SBATCH --error=metabat.err

module load miniconda 
conda activate binning 

mkdir METABAT

runMetaBat.sh -i contigs.fasta SAMPLE_BOTH_reads_unmapped.sorted.bam -o METABAT

runMetaBat.sh contigs.fasta SAMPLE_BOTH_reads_unmapped.sorted.bam -o METABAT

jgi_summarize_bam_contig_depths --outputDepth depth.txt *.bam

metabat2 -i contigs.fasta -a depth.txt -o METABAT
```

### Taxonomy

1. Next, it's time to assign taxonomy. There are two ways to do this a) assign to paired reads or b) assign to contigs. I did both, but we first have to create a database using Kraken2. 

```
#!/bin/bash
#SBATCH --job-name=kraken_DB
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=24:00:00
#SBATCH --output=kraken_DB.out
#SBATCH --error=kraken_DB.err

module load kraken2/2.0.8-beta	 

kraken2-build --standard --threads 24 --db KRAKEN_DB
```

2. Now assign taxonomy to contigs. 

```
#!/bin/bash
#SBATCH --job-name=kraken
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=12
#SBATCH --time=24:00:00
#SBATCH --output=kraken.out
#SBATCH --error=kraken.err

module load miniconda 
conda activate beluga

kraken --preload --db $KRAKEN_DB

mkdir TAXONOMY_MAG

kraken2 --db $KRAKEN_DB --threads 12 -input contigs.fasta --output TAXONOMY_MAG/contigs.kraken --report TAXONOMY_MAG/contigs.report
```

3. Then assign taxonomy to reads. 

```
#!/bin/bash
#SBATCH --job-name=kraken
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=12
#SBATCH --time=48:00:00
#SBATCH --output=kraken.out
#SBATCH --error=kraken.err

module load kraken2

mkdir TAXONOMY_MAG

kraken2 --db $KRAKEN2_DB --threads 12 --paired SAMPLE_host_removed_R1.fastq.gz SAMPLE_host_removed_R2.fastq.gz --output TAXONOMY_MAG/contigs.kraken --report TAXONOMY_MAG/contigs.report
```
