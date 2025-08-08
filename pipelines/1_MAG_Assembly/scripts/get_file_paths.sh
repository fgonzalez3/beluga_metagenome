#!/bin/bash

# Directory where the raw reads are stored
raw_reads_dir="/vast/palmer/pi/turner/flg9/TLab/Sequencing_Data/Beluga_Metagenome_RawReads/100M_PE/concatenated_reads/*"

# Output TSV file
output_file="beluga_raw_reads.tsv"

# Write the header line
echo -e "sample_id\tr1\tr2" > $output_file

# Loop over the directories in raw_reads_dir
for sample_dir in $(ls -d $raw_reads_dir); do
    # Extract the base sample name
    base_sample_name=$(basename $sample_dir)

    # Construct the paths to the R1 and R2 files
    r1_path=$(ls $sample_dir/*R1_combined.fastq.gz)
    r2_path=$(ls $sample_dir/*R2_combined.fastq.gz)

    # Check if both R1 and R2 files exist
    if [[ -e $r1_path && -e $r2_path ]]; then
        # If both files exist, write a line to the TSV file
        echo -e "$base_sample_name\t$r1_path\t$r2_path" >> $output_file
    fi
done


# create empty output file
f = open("tsv/beluga_raw_reads.tsv", "w")

# create empty variable that for loop can write to each time it loops throgh a sample
file_path = ""

def get_amp_paths(files):
    """
    Create a file that contains the paths to the raw reads for each sample
    """

    with open("tsv/beluga_raw_reads.tsv", "w") as tsv_file, \
    open("raw_read_dir/*.fasta", "r") as short_reads:

        short_read_files = short_reads.read()        

    for line in short_read_files:

        # Extract the base sample name from the directory name
        base_name = sample_dir

        # Construct the paths to the R1 and R2 files
        r1_path = raw_reads_dir + sample_dir + "/*R1_001.fastq.gz"
        r2_path = raw_reads_dir + sample_dir + "/*R2_001.fastq.gz"

return base_name, r1_path, r2_path

raw_reads_dir = "/vast/palmer/pi/turner/flg9/turner_lab/seqs/"
print(get_amp_paths(raw_reads_dir))