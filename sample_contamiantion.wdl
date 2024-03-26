version 1.0

# Define a task for downsampling and reheadering BAM files
task DownsampleAndReheader {
  # Define the input parameters for the task 
    input {
    Float desired_final_coverage
    Float contaminant_proportion
    File main_bam_file
    File main_bam_index  # Adding input for the index of the main BAM file
    File contaminant_bam_file
    File contaminant_bam_index  # Adding input for the index of the contaminant BAM file
    String main_sample_name
    String contaminant_sample_name
    String docker_image = "us.gcr.io/broad-dsp-lrma/mosdepth:sz_v3252024"
    String? disk_space  # Optional input for disk space
    Int? cpu    # CPU is now an optional input
    String? memory  # Memory is now an optional input
  }

  # Bash commands to execute the script; ensure 'set -euo pipefail' for strict error handling
  command <<<
    set -euo pipefail
    
    # Run the downsampling and reheadering script with specified inputs
    downsample_reheader.sh "~{desired_final_coverage}" "~{contaminant_proportion}" "~{main_bam_file}" "~{contaminant_bam_file}" "~{main_sample_name}"
  >>>

  # Define runtime resources needed
  runtime {
    docker: docker_image  # Specify Docker image
    memory: select_first([memory, "8 GB"])  # Use provided memory or default to 8 GB
    cpu: select_first([cpu, 4])  # Use provided CPU or default to 4
    disks: select_first([disk_space, "local-disk 1200 HDD"])  # Use provided disk space or default to 500 GB
  }

  # Specify outputs generated by this task
  output {
    # output_main_bam and output_contaminant_bam are the downsampled BAM files
      File output_main_bam =  "~{main_sample_name}_~{desired_final_coverage}_downsampled.bam"
      File output_main_bam_index = "~{main_sample_name}_~{desired_final_coverage}_downsampled.bam.bai"
      File output_contaminant_bam = "~{contaminant_sample_name}_~{contaminant_proportion}_downsampled.bam"
      File output_contaminant_bam_index = "~{contaminant_sample_name}_~{contaminant_proportion}_downsampled.bam.bai"
      File reheadered_contaminant_bam = "~{main_sample_name}_~{contaminant_proportion}_reheadered.bam"
      File reheadered_contaminant_bam_index = "~{main_sample_name}_~{contaminant_proportion}_reheadered.bam.bai"
  }
}

# Define a task for merging BAM files
task MergeBams {
  input {
    File main_bam                        # Main downsampled BAM file
    File contaminant_bam                 # Reheadered downsampled contaminant BAM file
    String output_filename               # Filename for the merged output BAM
    String docker_image = "us.gcr.io/broad-dsp-lrma/mosdepth:sz_v3252024"  # Docker image with samtools
    String? disk_space                   # Optional input for disk space
    Int? cpu                             # CPU is now an optional input
    String? memory                       # Memory is now an optional input
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
    memory: select_first([memory, "8 GB"])  # Use provided memory or default to 8 GB
    cpu: select_first([cpu, 4])  # Use provided CPU or default to 4
    disks: select_first([disk_space, "local-disk 12000 HDD"])  # Use provided disk space or default to 500 GB
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
    File main_bam_index  # Include this
    File contaminant_bam_file
    File contaminant_bam_index  # Include this
    String main_sample_name
    String contaminant_sample_name
    String merged_output_filename
    String? disk_space
    Int? cpu
    String? memory  # Memory is now an optional input
  }

    # Default filename if 'merged_output_filename' is not provided (work on this later)
  #  String default_merged_output_filename = "~{main_sample_name}_mix_~{contaminant_proportion}_~{contaminant_sample_name}.bam"
  #  String final_output_filename = select_first([merged_output_filename, default_merged_output_filename])

  # Call the downsampling and reheadering task
  call DownsampleAndReheader {
    input:
      desired_final_coverage = desired_final_coverage,
      contaminant_proportion = contaminant_proportion,
      main_bam_file = main_bam_file,
      main_bam_index = main_bam_index,  # Pass the main BAM index
      contaminant_bam_file = contaminant_bam_file,
      contaminant_bam_index = contaminant_bam_index,  # Pass the contaminant BAM index
      main_sample_name = main_sample_name,
      contaminant_sample_name = contaminant_sample_name,
      docker_image = "us.gcr.io/broad-dsp-lrma/mosdepth:sz_v3252024",
      disk_space = disk_space,
      cpu = cpu,  # Set CPU to 4
      memory = memory  # Set memory to 8 GB
  }

  # Call the merging task using outputs from the previous task
  call MergeBams {
    input:
      main_bam = DownsampleAndReheader.output_main_bam,
      contaminant_bam = DownsampleAndReheader.reheadered_contaminant_bam,
      output_filename = merged_output_filename,
      docker_image = "us.gcr.io/broad-dsp-lrma/mosdepth:sz_v3252024"
  }

  # Define the outputs of the entire workflow
  output {
    File downsampled_main_bam = DownsampleAndReheader.output_main_bam
    File downsampled_main_bam_index = DownsampleAndReheader.output_main_bam_index  # Adding index output
    File downsampled_contaminant_bam = DownsampleAndReheader.output_contaminant_bam
    File downsampled_contaminant_bam_index = DownsampleAndReheader.output_contaminant_bam_index  # Adding index output
    File reheadered_bam = DownsampleAndReheader.reheadered_contaminant_bam
    File reheadered_bam_index = DownsampleAndReheader.reheadered_contaminant_bam_index  # Adding index output
    File merged_bam = MergeBams.merged_bam
    File merged_bam_index = MergeBams.merged_bam_index
  }
}
