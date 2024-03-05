# sz_downsampling
To use: Please run the following scripts in the order indicated.

1- `downsampling.sh` : `bash downsampling.sh HG00512.bam 0.8137 HG00731.bam 0.0114 &> downsample.sh.log &`


2- `modify_header.sh` : `bash modify_header.sh HG00731_0227_downsampled.bam HG00512 HG00512_0227_downsampled &> modify_header_HG00512_0227_downsampled_.sh.log &`


3- `merge_bam_outputs.sh` : `bash merge_bam_outputs.sh merged_output.bam HG00512.bam HG00731.bam HG00514.bam &> merge_bam_outputs.sh.log &`

