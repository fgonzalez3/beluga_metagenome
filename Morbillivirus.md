# De-duplicating contigs

First, de-duplicate assembly. This will hurry any efforts in the future. 

```
#!/bin/bash
#SBATCH --job-name=CD-HIT
#SBATCH --nodes=2
#SBATCH --time=24:00:00
#SBATCH --mem=30gb
#SBATCH --output=CD-HIT.out
#SBATCH --error=CD-HIT.err

ml miniconda 
conda activate morbillivirus

cd-hit-est -i contigs.fasta -o contigs_DEDUP.fasta -M 0 -T 0
```

# Make a BLAST database

There were only 7 Morbillivirus reference sequences, so this DB is pretty small but should yield some hits if they exist. 

```
#!/bin/bash
#SBATCH --job-name=morbillivirus_DB
#SBATCH --nodes=1
#SBATCH --time=24:00:00
#SBATCH --mem=60gb
#SBATCH --output=morbillivirus_DB.out
#SBATCH --error=morbillivirus_DB.err

module load miniconda
conda activate kraken

makeblastdb -in genbank_morbillivirus.fasta -dbtype nucl -parse_seqids -out morbillivirus_nt_db
```

# BLAST

Next, I ran BLAST using my custom DB. I did this in the command line: 

```
blastn -word_size 10 -evalue 0.001 -query contigs_DEDUP.fasta -db morbillivirus_nt_db -outfmt '6 qseqid nident pident length evalue bitscore sgi sacc stitle'  -max_target_seqs 10 -out morbillivirus_blast_nt.txt

cat morbillivirus_blast_nt.txt | awk '{print $1}' | sort | uniq > morbillivirus_unique_contigs.txt
cat morbillivirus_blast_nt.txt | awk -F\t '($4>99 && $5<0.00001)' > morbillivirus_blast_nt_100len5eval.txt
cat morbillivirus_blast_nt.txt | awk -F\t '($4>99 && $6>99)' > morbillivirus_blast_nt_100len100bit.txt
cat morbillivirus_blast_nt_100len5eval.txt | awk '{print $1}' | sort | uniq > morbillivirus_unique_contigs_feces_nt_hiqual.txt
cat morbillivirus_unique_contigs_feces_nt_hiqual.txt | awk -F\_ '{print $1"_"$2}' | sort | uniq > morbillivirus_unique_sampleID_feces_nt_hiqual.txt
```

This yielded hits ~25bp in length, likely due to some weird sequence homology or viral integration. We'll need to run RNAseq to determine if there are any true hits though. 

