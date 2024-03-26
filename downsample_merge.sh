#!/bin/bash
set -euxo pipefail

# Combined Script: Downsampling, Reheadering, and Merging BAM Files

# Validate the number of arguments for downsampling, reheadering, and merging
if [ "$#" -ne 6 ]; then
    echo "Usage: $0 [Desired Final Coverage] [Contaminant Proportion] [Main BAM File] [Contaminant BAM File] [Main Sample Name] [Merged BAM Output]"
    exit 1
fi

# Assign variables from arguments for downsampling, reheadering, and merged output name
desired_final_coverage="$1"
contaminant_proportion="$2"
main_bam_file="$3"
contaminant_bam_file="$4"
main_sample_name="$5"
merged_bam="$6"

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
    echo "Calculating coverage for $bam_file..." >&2
    mosdepth -t 2 -n -x -Q 1 "${prefix}" "$bam_file"
    local mean_coverage=$(awk '$1 == "total" {print $4}' "${prefix}.mosdepth.summary.txt")
    echo >&2 "Original coverage for $bam_file is $mean_coverage"
    echo "$mean_coverage"
}

echo "Calculating original coverages..."
original_coverage_main=$(calculate_coverage "$main_bam_file")
original_coverage_contaminant=$(calculate_coverage "$contaminant_bam_file")

echo "Original coverage - Main: $original_coverage_main, Contaminant: $original_coverage_contaminant"
if [[ -z "$original_coverage_main" || -z "$original_coverage_contaminant" ]]; then
    echo "Error: Could not calculate original coverages."
    exit 1
fi

# Calculate downsampling percentages
echo "Calculating downsampling percentages..."
downsample_percentage_main=$(echo "scale=6; $desired_final_coverage * (1 - $contaminant_proportion) / $original_coverage_main" | bc)
downsample_percentage_contaminant=$(echo "scale=6; $desired_final_coverage * $contaminant_proportion / $original_coverage_contaminant" | bc)

echo "Downsampling percentages - Main: $downsample_percentage_main, Contaminant: $downsample_percentage_contaminant"

# Perform downsampling and indexing for the main BAM file
# Modification: Use desired_final_coverage in the file name instead of downsample_percentage_main
output_main_bam="${main_bam_file%.bam}_${desired_final_coverage}_downsampled.bam"
echo "Downsampling $main_bam_file to $downsample_percentage_main..."
samtools view -b -s $downsample_percentage_main $main_bam_file > $output_main_bam
samtools index $output_main_bam
echo "Recalculating coverage for $output_main_bam..."
new_coverage_main=$(calculate_coverage "$output_main_bam")
echo "New coverage for $output_main_bam is $new_coverage_main"


# Perform downsampling and indexing for the contaminant BAM file
# Modification: Use contaminant_proportion in the file name instead of downsample_percentage_contaminant
output_contaminant_bam="${contaminant_bam_file%.bam}_${contaminant_proportion}_downsampled.bam"
echo "Downsampling $contaminant_bam_file to $downsample_percentage_contaminant..."
samtools view -b -s $downsample_percentage_contaminant $contaminant_bam_file > $output_contaminant_bam
samtools index $output_contaminant_bam
echo "Recalculating coverage for $output_contaminant_bam..."
new_coverage_contaminant=$(calculate_coverage "$output_contaminant_bam")
echo "New coverage for $output_contaminant_bam is $new_coverage_contaminant"

# Modifying the header for the downsampled contaminant BAM file
out_prefix="${main_sample_name}_${contaminant_proportion}_reheadered"
echo "Modifying header for $output_contaminant_bam..."
samtools view --no-PG -H $output_contaminant_bam > header_"$out_prefix".txt
grep "^@RG" header_"$out_prefix".txt > rg_lines_"$out_prefix".txt
if ! grep -qF "SM:" rg_lines_"$out_prefix".txt; then
    sed -i "s/$/SM:tbd/" rg_lines_"$out_prefix".txt
fi

# Use awk to replace the sample name with the main sample name
awk -v sm="$main_sample_name" -F '\t' 'BEGIN {OFS="\t"} { for (i=1; i<=NF; ++i) { if ($i ~ /^SM:/) $i="SM:"sm } print }' rg_lines_"$out_prefix".txt > fixed_rg_lines_"$out_prefix".txt
grep -v "^@RG" header_"$out_prefix".txt > otherlines_"$out_prefix".txt
cat otherlines_"$out_prefix".txt fixed_rg_lines_"$out_prefix".txt > fixed_header_"$out_prefix".txt
samtools reheader fixed_header_"$out_prefix".txt $output_contaminant_bam > "${out_prefix}.bam"
samtools index "${out_prefix}.bam"

# Report final output names for clarity
echo "Final output files:"
echo "Main BAM: $output_main_bam"
echo "Main BAM Index: ${output_main_bam}.bai"
echo "Contaminant BAM: $output_contaminant_bam"
echo "Contaminant BAM Index: ${output_contaminant_bam}.bai"
echo "Reheadered Contaminant BAM: ${out_prefix}.bam"
echo "Reheadered Contaminant BAM Index: ${out_prefix}.bam.bai"

echo "Downsampling and reheadering process completed successfully."


# Confirming the files to merge
echo "Preparing to merge the following BAM files:"
echo "Main BAM: $output_main_bam"
echo "Reheadered Contaminant BAM: ${out_prefix}.bam"
# check if the files exist
if [ ! -f "$output_main_bam" ]; then
    echo "Error: Main BAM file $output_main_bam does not exist."
    exit 1
fi

if [ ! -f "${out_prefix}.bam" ]; then
    echo "Error: Reheadered Contaminant BAM file ${out_prefix}.bam does not exist."
    exit 1
fi

# Merge the downsampled main and reheadered contaminant BAM files
echo "Merging downsampled BAM files..."
samtools merge "$merged_bam" "$output_main_bam" "${out_prefix}.bam"
echo "Merging completed: $merged_bam"
# Index the merged BAM file
samtools index "$merged_bam"
echo "Indexing completed: ${merged_bam}.bai"

#print success message
echo "BAM files have been successfully merged and indexed."