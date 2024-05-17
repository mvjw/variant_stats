import os
import argparse
import csv
import allel
import numpy as np
import pandas as pd
import warnings

def parse_args():
    parser = argparse.ArgumentParser(description='Count vartrix entries.')
    parser.add_argument('vartrixdir', help='path to directory containing vartrix files (alt_matrix.mtx, ref_matrix.mtx, variants.tsv)')
    parser.add_argument('vcf', help='VCF file corresponding to the vartrix call')
    parser.add_argument('out', help='Output file')
    args = parser.parse_args()
    return args

def count_vartrix(vardir, vcfdir):
    variants = os.path.join(vardir, "variants.tsv")
    ref_mtx = os.path.join(vardir, "ref_matrix.mtx")
    alt_mtx = os.path.join(vardir, "alt_matrix.mtx")
    variant_list = []
    chromosomes_vartrix = set()
    with open(variants) as f:
        tsvreader = csv.reader(f)
        for row in tsvreader:
            variant_name = row[0]
            if 'chr' in variant_name:
              variant_name = variant_name.split('chr')[1]
            variant_list.append(variant_name)
            chromosomes_vartrix.add('_'.join(variant_name.split('_')[0:-1]))
    cell_list = [set() for _ in range(len(variant_list))]
    cells_w_alt_list = [set() for _ in range(len(variant_list))]
    ref_list = np.zeros(len(variant_list))
    with open(ref_mtx) as f:
        for i, line in enumerate(f.readlines()):
            if i < 3:
                continue
            entry = line.split(' ')
            ref_list[int(entry[0])-1] = ref_list[int(entry[0])-1] + int(entry[2])
            if int(entry[2]) > 0:
              cell_list[int(entry[0])-1].add(int(entry[1]))
    alt_list = np.zeros(len(variant_list))
    with open(alt_mtx) as f:
        for i, line in enumerate(f.readlines()):
            if i < 3:
                continue
            entry = line.split(' ')
            alt_list[int(entry[0])-1] = alt_list[int(entry[0])-1] + int(entry[2])
            if int(entry[2]) > 0:
                cell_list[int(entry[0])-1].add(int(entry[1]))
                cells_w_alt_list[int(entry[0])-1].add(int(entry[1]))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        gt_variants = allel.vcf_to_dataframe(vcfdir, fields='*', alt_number=1)
    chromosomes = list(gt_variants['CHROM'])
    # make sure chromosome namings are compatible
    for i, chromosome in enumerate(chromosomes):
        if 'chr' in chromosome:
            chromosome_split = chromosome.split('chr')[1]
            chromosomes[i] = chromosome_split
    missing_chr = []
    for chrom in list(chromosomes_vartrix):
      if not chrom in chromosomes:
        missing_chr.append(chrom)
    if len(missing_chr) > 0:
      print("### WARNING: Chromosomes {} not found in vcf file, will ignore".format(missing_chr))
    positions = list(gt_variants['POS'])
    ref_alleles = list(gt_variants['REF'])
    alt_alleles = list(gt_variants['ALT'])
    map_ind = {str(chromosome) + '_' + str(int(position)-1):i for i, (chromosome, position) in enumerate(zip(chromosomes, positions))}
    entries = []
    for variant in variant_list:
        try:
          index = map_ind[variant]
        except KeyError:
          continue
        loci = str(chromosomes[index]) + ':' + str(positions[index]) + ':' + str(ref_alleles[index]) + '->' + str(alt_alleles[index])
        allele = str(ref_alleles[index]) + '->' + str(alt_alleles[index])
        entry = {'map':loci, 'chromosome':str(chromosomes[index]), 'position':str(positions[index]), 'allele':allele, 'ref_count':ref_list[index], 'alt_count':alt_list[index], 'cell_count':len(cell_list[index]), 'cells_w_alt_count':len(cells_w_alt_list[index])}
        entries.append(entry)
    counts = pd.DataFrame(entries)
    return counts

def main():
    args = parse_args()
    vardir = args.vartrixdir
    vcfdir = args.vcf
    outfile = args.out
    counts = count_vartrix(vardir, vcfdir)
    counts.to_csv(outfile, sep='\t', index=False)

if __name__ == "__main__":
    main()