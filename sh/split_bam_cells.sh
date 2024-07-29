#!/bin/bash
# to run: source path/to/split_bam_cells.sh <path_to_bam> <name_of_output_directory> <path_to_barcodes> <nthreads>
# this script will not generate a list of cells
# the file "split_script3.py" must be located in the same folder as this bash source file
# nthreads refers to the number of threads which will be used during sorting
# requires python3 with pysam (and other dependencies) and samtools
declare -a path_to_bam=$1
declare -a output_dir=$2
declare -a path_to_barcodes=$3
declare -a nthreads=$4
declare -a base_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
echo $base_dir
mkdir ${output_dir}
samtools sort -l 0 -t "CB" -O sam -o ${output_dir}/sorted.sam -@ "${nthreads}" ${path_to_bam}
echo "done sorting"
mkdir ${output_dir}/split_bams
python3 ${base_dir}/split_script4.py -b "${output_dir}/sorted.sam" -o "${output_dir}/split_bams" -c "${path_to_barcodes}"
echo "done splitting"
rm ${output_dir}/sorted.sam
# question: do we need to index?
for inf in ${output_dir}/split_bams/*.bam; do samtools index -b $inf; done
# add a link to the original bam file in the , as well as to the barcode file
echo "done indexing"

