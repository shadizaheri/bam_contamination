version 1.0

# Define a task for downsampling and reheadering BAM files
task DownsampleAndReheaderMergeBams {
  # Define the input parameters for the task 
    input {
    Float desired_final_coverage
    Float contaminant_proportion
    File main_bam_file
    File main_bam_index  # Adding input for the index of the main BAM file
    File contaminant_bam_file
    File contaminant_bam_index  # Adding input for the index of the contaminant BAM file
    String main_sample_name
    String merged_output_filename               # Filename for the merged output BAM
    String docker_image = "us.gcr.io/broad-dsp-lrma/mosdepth:sz_v3252024"
    String? disk_space  # Optional input for disk space
    Int? cpu    # CPU is now an optional input
    String? memory  # Memory is now an optional input
  }

  # Bash commands to execute the script; ensure 'set -euo pipefail' for strict error handling
# Bash commands to execute the script; ensure 'set -euo pipefail' for strict error handling
command <<<
#!/bin/bash
set -euxo pipefail

# Combined Script: Downsampling, Reheadering, and Merging BAM Files
# Assign input variables for clarity and direct use
desired_final_coverage=${desired_final_coverage}
contaminant_proportion=${contaminant_proportion}
main_bam_file=${main_bam_file}
main_bam_index=${main_bam_index}
contaminant_bam_file=${contaminant_bam_file}
contaminant_bam_index=${contaminant_bam_index}
main_sample_name=${main_sample_name}
merged_bam=${merged_output_filename}

# commnad to run the script
# ./downsample_merge.sh desired_final_coverage contaminant_proportion main_bam_file contaminant_bam_file main_sample_name merged_bam
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
    mosdepth -t ${cpu} -n -x -Q 1 "${prefix}" "$bam_file"
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
# check if the reheadered file exists
if [ ! -f "${out_prefix}.bam" ]; then
    echo "Error: Reheadered BAM file ${out_prefix}.bam does not exist."
    exit 1
fi
echo "Reheadered BAM file created: ${out_prefix}.bam"
echo "Indexing reheadered BAM file..."

samtools index "${out_prefix}.bam"
echo "Indexing completed: ${out_prefix}.bam.bai"
# check if the index file exists
if [ ! -f "${out_prefix}.bam.bai" ]; then
    echo "Error: Index file ${out_prefix}.bam.bai does not exist."
    exit 1
fi
# clean up intermediate files
rm header_"$out_prefix".txt rg_lines_"$out_prefix".txt fixed_rg_lines_"$out_prefix".txt otherlines_"$out_prefix".txt fixed_header_"$out_prefix".txt



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
echo "Merging BAM files..."
# Merge the downsampled main and reheadered contaminant BAM files
echo "Merging downsampled BAM files..."
samtools merge "$merged_bam" "$output_main_bam" "${out_prefix}.bam"
echo "Merging completed: $merged_bam"
# Index the merged BAM file
samtools index "$merged_bam"
echo "Indexing completed: ${merged_bam}.bai"

#print success message
echo "BAM files have been successfully merged and indexed."
>>>

  # Define runtime resources needed
  runtime {
    docker: docker_image  # Specify Docker image
    memory: select_first([memory, "8 GB"])  # Use provided memory or default to 8 GB
    cpu: select_first([cpu, 4])  # Use provided CPU or default to 4
    disks: select_first([disk_space, "local-disk 800 HDD"])  # Use provided disk space or default to 500 GB
  }
  # Specify outputs generated by this task
  output {
      File output_main_bam = sub(main_bam_file, "\\.bam$", "") + "_~{desired_final_coverage}_downsampled.bam"
      File output_main_bam_index = sub(main_bam_file, "\\.bam$", "") + "_~{desired_final_coverage}_downsampled.bam.bai"
      File output_contaminant_bam = sub(contaminant_bam_file, "\\.bam$", "") + "_~{contaminant_proportion}_downsampled.bam"
      File output_contaminant_bam_index = sub(contaminant_bam_file, "\\.bam$", "") + "_~{contaminant_proportion}_downsampled.bam.bai"
      File reheadered_contaminant_bam = "~{main_sample_name}_~{contaminant_proportion}_reheadered.bam"
      File reheadered_contaminant_bam_index = "~{main_sample_name}_~{contaminant_proportion}_reheadered.bam.bai"
      File merged_bam = "~{merged_output_filename}.bam"
      File merged_bam_index = "~{merged_output_filename}.bam.bai"
  }
}

# Define a workflow that combines the downsampling/reheadering and merging tasks
workflow ContaminationWorkflow {
  input {
    Float desired_final_coverage
    Float contaminant_proportion
    File main_bam_file
    File main_bam_index  # Include this
    File contaminant_bam_file
    File contaminant_bam_index  # Include this
    String main_sample_name
    String merged_output_filename
    String? disk_space
    Int? cpu
    String? memory  # Memory is now an optional input
  }


  # Call the downsampling and reheadering task
  call DownsampleAndReheaderMergeBams {
    input:
      desired_final_coverage = desired_final_coverage,
      contaminant_proportion = contaminant_proportion,
      main_bam_file = main_bam_file,
      main_bam_index = main_bam_index,  # Pass the main BAM index
      contaminant_bam_file = contaminant_bam_file,
      contaminant_bam_index = contaminant_bam_index,  # Pass the contaminant BAM index
      main_sample_name = main_sample_name,
      merged_output_filename = merged_output_filename,
      docker_image = "us.gcr.io/broad-dsp-lrma/mosdepth:sz_v3252024",
      disk_space = disk_space,
      cpu = cpu,  # Set CPU to 4
      memory = memory  # Set memory to 8 GB
  }

  # Define the outputs of the entire workflow
  output {
    # make the next output optional
    #File downsampled_main_bam = DownsampleAndReheaderMergeBams.output_main_bam
    #File downsampled_main_bam_index = DownsampleAndReheaderMergeBams.output_main_bam_index  # Adding index output
    #File downsampled_contaminant_bam = DownsampleAndReheaderMergeBams.output_contaminant_bam
    #File downsampled_contaminant_bam_index = DownsampleAndReheaderMergeBams.output_contaminant_bam_index  # Adding index output
    #File reheadered_bam = DownsampleAndReheaderMergeBams.reheadered_contaminant_bam
    #File reheadered_bam_index = DownsampleAndReheaderMergeBams.reheadered_contaminant_bam_index  # Adding index output
    File merged_bam = DownsampleAndReheaderMergeBams.merged_bam
    File merged_bam_index = DownsampleAndReheaderMergeBams.merged_bam_index
  }
}
