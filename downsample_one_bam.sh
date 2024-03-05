#!/bin/bash
set -eu

# Check if the correct number of arguments is given
if [ "$#" -ne 3 ]; then
    echo "Usage: $0 [Input BAM] [Downsample Percentage] [Output Directory]"
    exit 1
fi

# Assign variables from arguments
input_bam="$1"
downsample_percentage="$2"
output_dir="$3"  # Assign the third argument to output_dir

# Function to be called when an error occurs
error_handling() {
    echo "An error occurred. Exiting the script."
    exit 1
}

# Trap to call error_handling function on error
trap error_handling ERR

# Create the output directory
mkdir -p "$output_dir"

# Downsampling and indexing for the BAM file
output_bam="${input_bam%.bam}_${downsample_percentage}_downsampled.bam"
echo "Downsampling $input_bam..."
samtools view -b -s $downsample_percentage $input_bam > "$output_dir/$output_bam"
echo "$input_bam downsampled to $downsample_percentage%."

echo "Indexing $output_bam..."
samtools index "$output_dir/$output_bam"
echo "$output_bam indexed and moved to $output_dir."

echo "Downsampling and indexing complete. Output file has been moved to $output_dir."
