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

