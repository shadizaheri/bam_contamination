version 1.0

# Define a task for downsampling and reheadering BAM files
task DownsampleAndReheader {
  input {
    Float desired_final_coverage       # Desired coverage after downsampling
    Float contaminant_proportion       # Proportion of contaminant
    File main_bam_file                 # Main BAM file to be downsampled
    File contaminant_bam_file          # Contaminant BAM file to be downsampled and reheadered
    String main_sample_name            # Sample name for reheadering
    String docker_image = "us.gcr.io/broad-dsp-lrma/mosdepth:sz_v3192024"  # Docker image containing necessary tools
  }

  # Bash commands to execute the script; ensure 'set -euo pipefail' for strict error handling
  command <<<
    set -euo pipefail
    
    # Run the downsampling and reheadering script with specified inputs
    downsample.sh "~{desired_final_coverage}" "~{contaminant_proportion}" "~{main_bam_file}" "~{contaminant_bam_file}" "~{main_sample_name}"
  >>>

  # Define runtime resources needed
  runtime {
    docker: docker_image  # Specify Docker image
    memory: "8 GB"        # Allocate memory; adjust based on needs and environment
    cpu: "4"              # Allocate CPU; adjust as necessary
  }

  # Specify outputs generated by this task
  output {
    File output_main_bam = sub(main_bam_file, "\\.bam$", "") + "_" + "~{desired_final_coverage}" + "_downsampled.bam"
    File output_contaminant_bam = sub(contaminant_bam_file, "\\.bam$", "") + "_" + "~{contaminant_proportion}" + "_downsampled.bam"
    File reheadered_contaminant_bam = "~{main_sample_name}" + "_" + "~{contaminant_proportion}" + "_reheadered.bam"
  }
}

# Define a task for merging BAM files
task MergeBams {
  input {
    File main_bam                        # Main downsampled BAM file
    File contaminant_bam                 # Reheadered downsampled contaminant BAM file
    String output_filename               # Filename for the merged output BAM
    String docker_image = "us.gcr.io/broad-dsp-lrma/mosdepth:sz_v3192024"  # Docker image with samtools
  }

  # Bash commands to execute the merging script
  command <<<
    set -eu
    
    # Merge the specified BAM files and index the merged BAM
    merge.sh "~{output_filename}" "~{main_bam}" "~{contaminant_bam}"
  >>>

  # Define runtime resources needed
  runtime {
    docker: docker_image
    memory: "8 GB"
    cpu: "2"
  }

  # Specify outputs generated by this task
  output {
    File merged_bam = "~{output_filename}"
    File merged_bam_index = "~{output_filename}.bai"
  }
}

# Define a workflow that combines the downsampling/reheadering and merging tasks
workflow ContaminationWorkflow {
  input {
    Float desired_final_coverage
    Float contaminant_proportion
    File main_bam_file
    File contaminant_bam_file
    String main_sample_name
    String merged_output_filename        # Desired name for the final merged BAM file
  }

  # Call the downsampling and reheadering task
  call DownsampleAndReheader {
    input:
      desired_final_coverage = desired_final_coverage,
      contaminant_proportion = contaminant_proportion,
      main_bam_file = main_bam_file,
      contaminant_bam_file = contaminant_bam_file,
      main_sample_name = main_sample_name,
      docker_image = "us.gcr.io/broad-dsp-lrma/mosdepth:sz_v3192024"
  }

  # Call the merging task using outputs from the previous task
  call MergeBams {
    input:
      main_bam = DownsampleAndReheader.output_main_bam,
      contaminant_bam = DownsampleAndReheader.reheadered_contaminant_bam,
      output_filename = merged_output_filename,
      docker_image = "us.gcr.io/broad-dsp-lrma/mosdepth:sz_v3192024"
  }

  # Define the outputs of the entire workflow
  output {
    File downsampled_main_bam = DownsampleAndReheader.output_main_bam
    File downsampled_contaminant_bam = DownsampleAndReheader.output_contaminant_bam
    File reheadered_bam = DownsampleAndReheader.reheadered_contaminant_bam
    File merged_b = MergeBams.merged_bam
    File merged_bam_index = MergeBams.merged_bam_index
  }
}