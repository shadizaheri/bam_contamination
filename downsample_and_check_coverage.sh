#!/bin/bash
set -eu

# Check if the correct number of arguments is given
if [ "$#" -ne 4 ]; then
    echo "Usage: $0 [Input BAM 1] [Downsample Percentage 1] [Input BAM 2] [Downsample Percentage 2]"
    exit 1
fi

# Assign variables from arguments
input_bam_1="$1"
downsample_percentage_1="$2"
input_bam_2="$3"
downsample_percentage_2="$4"

# Function to be called when an error occurs
error_handling() {
    echo "An error occurred. Exiting the script."
    exit 1
}

# Trap to call error_handling function on error
trap error_handling ERR

# Create a prefix for mosdepth output
prefix_1="${input_bam_1%.bam}_ds${downsample_percentage_1}"
prefix_2="${input_bam_2%.bam}_ds${downsample_percentage_2}"

# Downsampling and calculating coverage for the first BAM file
output_bam_1="${prefix_1}.bam"
echo "Downsampling $input_bam_1..."
samtools view -b -s $downsample_percentage_1 $input_bam_1 > $output_bam_1
echo "$input_bam_1 downsampled to $downsample_percentage_1%."
echo "Calculating coverage for $output_bam_1..."
mosdepth -t 2 -n -x -Q 1 $prefix_1 $output_bam_1

# Downsampling and calculating coverage for the second BAM file
output_bam_2="${prefix_2}.bam"
echo "Downsampling $input_bam_2..."
samtools view -b -s $downsample_percentage_2 $input_bam_2 > $output_bam_2
echo "$input_bam_2 downsampled to $downsample_percentage_2%."
echo "Calculating coverage for $output_bam_2..."
mosdepth -t 2 -n -x -Q 1 $prefix_2 $output_bam_2

echo "Downsampling and coverage calculation complete."
