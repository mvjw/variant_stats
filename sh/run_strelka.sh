#!/bin/bash
# options can be specified and passed to strelka, such as --exome or --rna or --min-qscore=20
declare -a options=$1
declare -a inf=$2
declare -a working=$3
declare -a temp=$4
declare -a reference=$5
declare -a njobs=$6
declare -a outf=$(basename "$inf" | cut -d "." -f 1)
mkdir "${temp}/strelka_temp/${outf}"
configureStrelkaGermlineWorkflow.py $options --bam $inf --referenceFasta $reference --runDir ${temp}/strelka_temp/${outf} > /dev/null
${temp}/strelka_temp/${outf}/runWorkflow.py -m local -j $njobs > /dev/null 2> /dev/null
mkdir ${working}/variants/${outf}
mv ${temp}/strelka_temp/${outf}/results/variants/* ${working}/variants/${outf}/
rm -r ${temp}/strelka_temp/${outf}/
