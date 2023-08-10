# China metagenome

I first navigated to the data used in the Bai et al paper to obtain the Beluga raw reads from that paper. I chose to only focus on the Beluga samples, though others were taken from dolphins and seals. Their data is freely available on NCBI following this [link](https://www.ncbi.nlm.nih.gov/biosample/SAMN20056375). 

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
