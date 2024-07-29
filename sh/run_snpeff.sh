#!/bin/bash
# run with ./path/to/run_snpeff.sh /path/to/snpeff_dir <path_to_vcf_file>
# where working dir is location of folder containing variants folder
# snpeff should not include forward slash on last character

set -e

declare -a snpEff_path=$1
declare -a infile=$2
declare -a reference="GRCm38.86"
declare outpath="$(dirname "${infile}")"
java -jar ${snpEff_path}/snpEff.jar -c ${snpEff_path}/snpEff.config -s ${snpEff_path}/snpEff_summary.html $reference $infile > "${outpath}/annotated.vcf"
