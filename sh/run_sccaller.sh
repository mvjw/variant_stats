#!/bin/bash
declare -a sc_path="${1}"
declare -a working="${2}"
declare -a reference="${3}"
declare -a path_to_file="${4}"
declare -a nparjobs="${5}"
declare -a bulk_vcf="${working}/bulk_variants/variants.vcf"
declare -a outf=$(basename "$path_to_file" | cut -d "." -f 1)
mkdir ${working}/sc_variants2/${outf}
python2.7 $sc_path si --bam $path_to_file --fasta $reference --output ${working}/sc_variants2/${outf}/variants.vcf --snp_type hsnp --snp_in $bulk_vcf --cpu_num $nparjobs --min_depth 5 --minvar 1 --engine samtools > /dev/null

