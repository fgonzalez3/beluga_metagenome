# China metagenome

I first navigated to the data used in the Bai et al paper to obtain the Beluga raw reads from that paper. I chose to only focus on the Beluga samples, though others were taken from dolphins and seals. Their data is freely available on NCBI following this [link](https://www.ncbi.nlm.nih.gov/biosample/SAMN20056375). After looking back at the paper, this was amplicon-based sequencing. No further steps were taken beyond QC. 

# De-merging paired-end reads 

The Beluga paired-end reads (SRR15037280_ANDY2.fastq.gz  SRR15037281_ANDY1.fastq.gz  SRR15037286_TINA3.fastq.gz  SRR15037287_TINA2.fastq.gz  SRR15037293_TINA1.fastq.gz  SRR15037304_ANDY3.fastq.gz) were merged, so I had to take an extra step here to de-merge them and get the forward/reverse reads necessary for downstream analysis: 

a) I first made a new Conda environment containing [BBMAP](https://anaconda.org/bioconda/bbmap), a part of the BBTools suite, which provides a tool called reformat.sh that can be used to split interleaved FastQ files into separate R1 and R2 files. 

b) I next ran the following script: 

```
#!/bin/bash
#SBATCH --job-name=unmerge_pairs
#SBATCH --nodes=1
#SBATCH --mem=50G
#SBATCH --time=01:00:00
#SBATCH --output=unmerge_pairs.out
#SBATCH --error=unmerge_pairs.err

module load miniconda 
conda activate unmerge 

reformat.sh in=SRR15037293_TINA1.fastq.gz out1=TINA1_R1.fastq out2=TINA1_R2.fastq
reformat.sh in=SRR15037287_TINA2.fastq.gz out1=TINA2_R1.fastq out2=TINA2_R2.fastq
reformat.sh in=SRR15037286_TINA3.fastq.gz out1=TINA3_R1.fastq out2=TINA3_R2.fastq
reformat.sh in=SRR15037281_ANDY1.fastq.gz out1=ANDY1_R1.fastq out2=ANDY1_R2.fastq
reformat.sh in=SRR15037280_ANDY2.fastq.gz out1=ANDY2_R1.fastq out2=ANDY2_R2.fastq
reformat.sh in=SRR15037304_ANDY3.fastq.gz out1=ANDY3_R1.fastq out2=ANDY3_R2.fastq
```

# Quality control 

I next did some quality control analyses on these samples. I mainly wanted to verify that the adapter sequences had been previously trimmed before proceeding. 

```
#!/bin/bash
#SBATCH --job-name=fastqc
#SBATCH --nodes=1
#SBATCH --time=12:00:00
#SBATCH --output=fastqc.out
#SBATCH --error=fastqc.err

module load miniconda 
conda activate beluga 

samples=("ANDY1" "ANDY2" "ANDY3" "TINA1" "TINA2" "TINA3")

for sample_name in "${samples[@]}"
do
    mkdir "$sample_name"
    fastqc "${sample_name}_R1.fastq" "${sample_name}_R2.fastq" -o "$sample_name"

    cd "$sample_name" # Enter the sample directory

    # Unzip FastQC results
    for filename in *.zip 
    do 
        unzip "$filename"
    done 

    # Create 'docs' directory and combine summary files
    mkdir docs
    cat */summary.txt > docs/fastqc_summaries.txt

    cd .. # Go back to the main directory
done
```

# Further Quality Control 

Running FASTQC showed numerous quality metric failures. Specifically, Per base sequence content, Per sequence GC content, Sequence duplication levels, and Overrepresented sequences metrics failed. To troubleshoot, I ran some additional quality control in order to make downstream analysis easier. I (a) first trimmed low-quality bases and (b) deduplicated reads. Afterward, I ran FASTQC on this parsed-through data. For my trimming of low-quality bases, I followed the procedure of the Beluga metagenome paper as closely as possible to replicate their results. That would be encompassed by using a sliding window of 5 and trimming bases with a <20 quality score, removing reads shorter than 200bp, and removing low-bases at the beginning and end of the reads (LEADING:3 and TRAILING:3). 

```
#!/bin/bash
#SBATCH --job-name=china_qc
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --time=12:00:00
#SBATCH --output=china_qc.out
#SBATCH --error=china_qc.err

ml miniconda 
conda activate beluga 

trimmomatic PE -threads 4 -phred33 ANDY1_R1.fastq ANDY1_R2.fastq ANDY1_R1.filtered.fastq ANDY1_R1.unfiltered.fastq ANDY1_R2.filtered.fastq ANDY1_R2.unfiltered.fastq SLIDINGWINDOW:5:20 MINLEN:200 LEADING:3 TRAILING:3

mkdir cd-hit 
mv *.filtered.fastq cd-hit 
cd cd-hit

conda activate morbillivirus 

cd-hit-est -i ANDY1_R1.filtered.fastq -o ANDY1_R1.filtered.DEDUP.fastq -T 4 -c 0.95 -n 8
cd-hit-est -i ANDY1_R2.filtered.fastq -o ANDY1_R2.filtered.DEDUP.fastq -T 4 -c 0.95 -n 8

conda activate beluga 

mkdir FASTQC 
mv *.DEDUP.fastq FASTQC
cd FASTQC

fastqc ANDY1_R1.filtered.DEDUP.fastq ANDY1_R2.filtered.DEDUP.fastq
```
