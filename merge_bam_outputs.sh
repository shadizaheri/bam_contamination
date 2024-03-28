#!/bin/bash
# Author: Shadi Zaheri (szaheri@broadinstitute.org)

set -eu

# Check if there are enough arguments
if [ "$#" -lt 3 ]; then
    echo "Usage: $0 [Merged BAM Output] [Input BAM 1] [Input BAM 2] ..."
    echo "Please provide at least two BAM files to merge and the name for the merged output file."
    exit 1
fi

# Define variable for the merged output from arguments
merged_bam="$1"

# Remove first argument so only BAM files remain
shift

# Check if the input BAM files exist
for bam_file in "$@"; do
    if [ ! -f "$bam_file" ]; then
        echo "Error: BAM file $bam_file does not exist."
        exit 1
    fi
done

# Echo merging command for transparency and confirmation
echo "Merging BAM files into $merged_bam..."
echo "Input files to merge: $@"

# Merge the BAM files using all remaining arguments
samtools merge "$merged_bam" "$@"

echo "Merging completed: $merged_bam"

# Index the merged BAM file
echo "Indexing merged BAM file..."
samtools index "$merged_bam"

echo "Indexing completed: ${merged_bam}.bai"

# Print a success message
echo "BAM files have been successfully merged and indexed."
