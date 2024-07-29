#!/bin/bash
# run with /path/to/run_vartrix.sh /path/to/vcf.vcf /path/to/metadata.tsv /path/to/bam.bam /path/to/barcodes.tsv /path/to/ref.fa /path/to/outdir <atac or rna> nthreads
declare -a vcf="${1}"
declare -a meta="${2}"
declare -a bam="${3}"
declare -a barcode="${4}"
declare -a ref="${5}"
declare -a outdir="${6}"
declare -a rna_opt="${7}"
declare -a nthreads="${8}"

export PATH=/<redacted>/bin/vartrix-1.1.16:$PATH

mkdir $outdir
cp $barcode ${outdir}/
if [ "$rna_opt" == "rna" ]; then
  vartrix_linux -v $vcf -b $bam -f $ref -s coverage --out-variants ${outdir}/variants.tsv --ref-matrix ${outdir}/ref_matrix.mtx --out-matrix ${outdir}/alt_matrix.mtx -c $barcode --threads $nthreads --mapq 10 --umi
elif [ "$rna_opt" == "atac" ]; then
  vartrix_linux -v $vcf -b $bam -f $ref -s coverage --out-variants ${outdir}/variants.tsv --ref-matrix ${outdir}/ref_matrix.mtx --out-matrix ${outdir}/alt_matrix.mtx -c $barcode --threads $nthreads --mapq 10
else
  echo "Specify rna or atac" 1>&2
  exit 64
fi
cp "${meta}" "${outdir}/variant_meta_data.tsv"
ln -s $(dirname $barcode) ${outdir}/levels

