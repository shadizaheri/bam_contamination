#!/bin/bash
set -euo pipefail

echo "Script initiated: Downsampling and Reheadering BAM Files"

# Validate the number of arguments
if [ "$#" -ne 5 ]; then
    echo "Usage: $0 [Desired Final Coverage] [Contaminant Proportion] [Main BAM File] [Contaminant BAM File] [Main Sample Name]"
    exit 1
fi

# Assign variables from arguments
desired_final_coverage="$1"
contaminant_proportion="$2"
main_bam_file="$3"
contaminant_bam_file="$4"
main_sample_name="$5"

echo "Parameters: Desired Coverage=$desired_final_coverage, Contaminant Proportion=$contaminant_proportion, Main BAM=$main_bam_file, Contaminant BAM=$contaminant_bam_file, Main Sample Name=$main_sample_name"

# Define function for error handling
error_handling() {
    echo "An error occurred. Exiting the script."
    exit 1
}

# Set trap for error handling
trap error_handling ERR

# Function to calculate original coverage using mosdepth
calculate_coverage() {
    local bam_file="$1"
    local prefix="${bam_file%.bam}"
    echo "Calculating coverage for $bam_file..."
    ./mosdepth -t 2 -n -x -Q 1 "${prefix}" "$bam_file"
    local mean_coverage=$(awk '{print $4}' "${prefix}.mosdepth.summary.txt")
    echo "Original coverage for $bam_file is $mean_coverage"
    echo "$mean_coverage"
}

echo "Calculating original coverages..."
# Calculate original coverages for each BAM file
original_coverage_main=$(calculate_coverage "$main_bam_file")
original_coverage_contaminant=$(calculate_coverage "$contaminant_bam_file")

echo "Original coverage - Main: $original_coverage_main, Contaminant: $original_coverage_contaminant"

# Calculate downsampling percentages
echo "Calculating downsampling percentages..."
downsample_percentage_main=$(echo "scale=2; $desired_final_coverage * (1 - $contaminant_proportion) / $original_coverage_main" | bc)
downsample_percentage_contaminant=$(echo "scale=2; $desired_final_coverage * $contaminant_proportion / $original_coverage_contaminant" | bc)

echo "Downsampling percentages - Main: $downsample_percentage_main, Contaminant: $downsample_percentage_contaminant"

# Perform downsampling and indexing for the main BAM file
output_main_bam="${main_bam_file%.bam}_${downsample_percentage_main}_downsampled.bam"
echo "Downsampling $main_bam_file to $downsample_percentage_main..."
samtools view -b -s $downsample_percentage_main $main_bam_file > $output_main_bam
samtools index $output_main_bam

# Perform downsampling and indexing for the contaminant BAM file
output_contaminant_bam="${contaminant_bam_file%.bam}_${downsample_percentage_contaminant}_downsampled.bam"
echo "Downsampling $contaminant_bam_file to $downsample_percentage_contaminant..."
samtools view -b -s $downsample_percentage_contaminant $contaminant_bam_file > $output_contaminant_bam
samtools index $output_contaminant_bam

# Modifying the header for the downsampled contaminant BAM file
out_prefix="${main_sample_name}_${contaminant_proportion}_reheadered"
echo "Modifying header for $output_contaminant_bam..."
samtools view --no-PG -H $output_contaminant_bam > header_"$out_prefix".txt
grep "^@RG" header_"$out_prefix".txt > rg_lines_"$out_prefix".txt
if ! grep -qF "SM:" rg_lines_"$out_prefix".txt; then
    sed -i "s/$/SM:tbd/" rg_lines_"$out_prefix".txt
fi
awk -v sm="$main_sample_name" -F '\t' 'BEGIN {OFS="\t"} { for (i=1; i<=NF; ++i) { if ($i ~ /^SM:/) $i="SM:"sm } print }' rg_lines_"$out_prefix".txt > fixed_rg_lines_"$out_prefix".txt
grep -v "^@RG" header_"$out_prefix".txt > otherlines_"$out_prefix".txt
cat otherlines_"$out_prefix".txt fixed_rg_lines_"$out_prefix".txt > fixed_header_"$out_prefix".txt
samtools reheader fixed_header_"$out_prefix".txt $output_contaminant_bam > "${out_prefix}.bam"
samtools index "${out_prefix}.bam"

echo "Header modification complete. New downsampled contaminant BAM file"
