#!/bin/bash
# run with source /path/to/post_process.sh <working_dir> /path/to/snpeff <njobs>
# where working dir is location of folder containing variants folder
# working should not include forward slash on last character
# snpeff should not include forward slash on last character
# running on 50 cores for 9200 cells will take approximately 5 hours
run_snpEff () {
  # $1 snpEff_path
  # $2 reference
  # $3 vcf_dir
  java -jar ${1}/snpEff.jar -c ${1}/snpEff.config -s ${3}/snpEff_summary.html $2 ${3}/variants.vcf > ${3}/annotated.vcf
}

# first step run snpEff
export -f run_snpEff
declare -a working=$1
declare -a reference="GRCm38.86"
declare -a snpEff_path=$2
declare -a njobs=$3
ls -d ${working}/variants/*/ | tr "\t" "\n" | parallel --gnu --jobs $njobs run_snpEff $snpEff_path $reference {}
ls -d ${working}/variants/*/ | tr "\t" "\n" > ${working}/cell_list.txt

# TODO - consider implementing --nostats and -t options to improve speedup, see manual http://snpeff.sourceforge.net/SnpEff_manual.html
