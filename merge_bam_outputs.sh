#!/bin/bash
set -eu

# Check if there are enough arguments
if [ "$#" -lt 3 ]; then
    echo "Usage: $0 [Merged BAM Output] [Input BAM 1] [Input BAM 2] ..."
    exit 1
fi

# Define variable for the merged output from arguments
merged_bam="$1"

# Remove first argument so only BAM files remain
shift

# Echo merging command for transparency
echo "Merging BAM files into $merged_bam..."

# Merge the BAM files using all remaining arguments
samtools merge "$merged_bam" "$@"

echo "Merging completed: $merged_bam"

# Index the merged BAM file
echo "Indexing merged BAM file..."
samtools index "$merged_bam"

echo "Indexing completed: ${merged_bam}.bai"

# Print a success message
echo "BAM files have been successfully merged and indexed."
