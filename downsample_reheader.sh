#!/bin/bash
set -euo pipefail

echo "Starting the downsampling and reheadering script..."

# Validate the number of arguments
if [ "$#" -ne 4 ]; then
    echo "Usage: $0 [Desired Final Coverage] [Contaminant Proportion] [BAM File 1] [BAM File 2]"
    exit 1
fi

# Assign variables from arguments
desired_final_coverage="$1"
contaminant_proportion="$2"
bam_file_1="$3"
bam_file_2="$4"

echo "Script parameters:"
echo "Desired final coverage: $desired_final_coverage"
echo "Contaminant proportion: $contaminant_proportion"
echo "BAM File 1: $bam_file_1"
echo "BAM File 2: $bam_file_2"

# Function to be called when an error occurs
error_handling() {
    echo "An error occurred. Exiting the script."
    exit 1
}

# Trap to call error_handling function on error
trap error_handling ERR

# Function to calculate original coverage using mosdepth
calculate_coverage() {
    local bam_file="$1"
    local prefix="${bam_file%.bam}"
    echo "Calculating original coverage for $bam_file using mosdepth..."
    mosdepth -t 2 -n -x -Q 1 "${prefix}" "$bam_file"
    # Extract mean coverage from mosdepth summary
    local mean_coverage=$(awk '{print $4}' "${prefix}.mosdepth.summary.txt")
    echo "Original coverage for $bam_file: $mean_coverage"
    echo "$mean_coverage"
}

# Calculate original coverages for each BAM file
echo "Calculating original coverages..."
original_coverage_1=$(calculate_coverage "$bam_file_1")
original_coverage_2=$(calculate_coverage "$bam_file_2")

# Calculate downsampling percentages
echo "Calculating downsampling percentages..."
proportion_main_sample=$(echo "scale=2; 1 - $contaminant_proportion" | bc)
downsample_percentage_1=$(echo "scale=2; $desired_final_coverage * $proportion_main_sample / $original_coverage_1" | bc)
downsample_percentage_2=$(echo "scale=2; $desired_final_coverage * $contaminant_proportion / $original_coverage_2" | bc)
echo "Downsampling percentage for BAM File 1: $downsample_percentage_1"
echo "Downsampling percentage for BAM File 2: $downsample_percentage_2"

# Downsampling and indexing for the BAM files
echo "Downsampling BAM files..."
output_bam_1="${bam_file_1%.bam}_${downsample_percentage_1}_downsampled.bam"
output_bam_2="${bam_file_2%.bam}_${downsample_percentage_2}_downsampled.bam"

echo "Downsampling $bam_file_1 to $downsample_percentage_1..."
samtools view -b -s $downsample_percentage_1 $bam_file_1 > $output_bam_1
echo "$bam_file_1 downsampled to $downsample_percentage_1%."
samtools index $output_bam_1

echo "Downsampling $bam_file_2 to $downsample_percentage_2..."
samtools view -b -s $downsample_percentage_2 $bam_file_2 > $output_bam_2
echo "$bam_file_2 downsampled to $downsample_percentage_2%."
samtools index $output_bam_2

# Additional checks for final downsampled BAM files
echo "Calculating final coverage for downsampled BAM files..."
final_coverage_1=$(calculate_coverage "$output_bam_1")
final_coverage_2=$(calculate_coverage "$output_bam_2")
echo "Final coverage for $output_bam_1: $final_coverage_1"
echo "Final coverage for $output_bam_2: $final_coverage_2"

echo "Downsampling and coverage calculation complete. Final coverages after downsampling: $final_coverage_1 (for $output_bam_1), $final_coverage_2 (for $output_bam_2)."
