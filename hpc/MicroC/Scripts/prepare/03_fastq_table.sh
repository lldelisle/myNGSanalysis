#!/bin/bash

# Script to generate the samplesFastqTable where:
# first column is the sample name
# second column is the relative fastq1 path
# third column is the relative fastq2 path

#################
#### SET UP #####
#################

# Define the parent directory containing sample folders
mkdir -p $microcPilot/fastq/
pathToFastq="$microcPilot/fastq/"

#################
#### SCRIPT #####
#################

# Initialize the output file
output_file="$pathToFastq/samplesFastqTable.txt"

cd $pathToFastq

# Iterate over each sample directory
for sample_dir in */; do
    # Extract the sample name (folder name)
    sample_name=$(basename "$sample_dir")

    # Find the FASTQ files ending with 1.fq.gz and 2.fq.gz
    read1_file=$(find "$sample_dir" -type f -name '*1.fq.gz')
    read2_file=$(find "$sample_dir" -type f -name '*2.fq.gz')

    # Append the sample information to the output file
    echo -e "$sample_name\t$read1_file\t$read2_file" >> "$output_file"
done