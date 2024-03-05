#!/bin/bash
set -eu

# Check if the correct number of arguments is given
if [ "$#" -ne 3 ]; then
    echo "Usage: $0 [local_bam] [new_sample_name] [out_prefix]"
    exit 1
fi

# Define variables from arguments
local_bam="$1"
new_sample_name="$2"
out_prefix="$3"

# Clean up the header
samtools view --no-PG -H "$local_bam" > header_"$out_prefix".txt

# Fix SM in the RG lines
grep "^@RG" header_"$out_prefix".txt > rg_lines_"$out_prefix".txt
if ! grep -qF "SM:" rg_lines_"$out_prefix".txt; then
    sed -i "s/$/SM:tbd/" rg_lines_"$out_prefix".txt
fi
awk -v sm="$new_sample_name" -F '\t' 'BEGIN {OFS="\t"} { for (i=1; i<=NF; ++i) { if ($i ~ /^SM:/) $i="SM:"sm } print}' \
    rg_lines_"$out_prefix".txt \
    > fixed_rg_lines_"$out_prefix".txt

# Paste things back
grep -v "^@RG" header_"$out_prefix".txt > otherlines_"$out_prefix".txt
cat otherlines_"$out_prefix".txt fixed_rg_lines_"$out_prefix".txt > fixed_header_"$out_prefix".txt

# Apply the new header to the BAM file
samtools reheader fixed_header_"$out_prefix".txt "$local_bam" > "${out_prefix}.bam"

# Index the reheadered BAM file
samtools index "${out_prefix}.bam"

echo "Header has been modified and applied to new BAM file: ${out_prefix}.bam"
echo "Indexed BAM file created: ${out_prefix}.bam.bai"
