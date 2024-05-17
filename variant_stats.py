import os
import sys
import pandas as pd
import json
import argparse
import copy

pd.set_option('display.max_columns', None)
from src.count_vartrix import count_vartrix
from src.count_ase_reader import count_ase_reader
#pd.options.mode.chained_assignment = None  # default='warn'

from src.parsed_vcf import ParsedVCF
from src.cell_list import CellList, get_cell_lists
from src.parsed_mutation_list import ParsedMutationList, get_mutations, save_mutation_array, load_mutation_array, generate_vcf
from src.filtering import filter_mutations
from src.overlap_mutation_list import OverlapMutationList, common_variants, redo_common, rank_overlap

### TEST SCRIPTS ###
    
def filter_single(infile, outfile, outvcf = None, save_intermediary = None, min_qual = -1.0, exonic=False, splice=False):
    print("Parsing vcf...")
    vcf_parsed = ParsedVCF(infile)
    vcf_parsed_as_list = [vcf_parsed]
    # TODO change cell list so that it may process single vcf without list hassle
    print("Compiling variants...")
    vcf_as_cell_list = CellList(vcf_parsed_as_list, "single_vcf")
    vcf_as_mutation_list = ParsedMutationList(vcf_as_cell_list)
    if save_intermediary:
      save_mutation_array(vcf_as_mutation_list, save_intermediary)
    # May wish to add additional options here such as Pass or scVAF
    print("Filtering variants...")
    filtered = filter_mutations(vcf_as_mutation_list, 0.0, 0, 0.0, 0, min_qual, exonic=exonic, splice=splice)[0]
    print("Saving variants...")
    filtered.parsed_mutations.to_csv(outfile, sep='\t', index=False)
    if outvcf:
      generate_vcf(filtered.parsed_mutations, filtered[0], outvcf)

def filter_and_overlap_multiple(sample_list, outfile, outvcf = None, save_intermediary = None, min_qual = -1.0):
    mutation_lists = []
    for sample in sample_list:
      vcf_parsed = ParsedVCF(sample['file'])
      vcf_as_cell_list = CellList([vcf_parsed], sample['name'])
      vcf_as_mutation_list = ParsedMutationList(vcf_as_cell_list)
      mutation_lists.append(vcf_as_mutation_list)
    if save_intermediary:
      save_mutation_array(mutation_lists, save_intermediary)
    filtered = filter_mutations(mutation_lists, 0.0, 0, 0.0, 0, min_qual, exonic=True, splice=False)
    overlap_list = OverlapMutationList(filtered)
    common = common_variants(overlap_list, filtered)
    common.to_csv(outfile, sep='\t', index=False)
    if outvcf:
      generate_vcf(common, filtered[0], outvcf)

def extract_scu_variants(lists_WGS, names_WGS, lists_CTRL, names_CTRL, lists_sc, names_sc, out_path, resume=False, min_cells_w_alt_count=0.0, min_alt_read_count=0, min_cell_vaf_factor=0.0, min_sc_qual=0.0, max_alt_read_count=None, max_cells_w_alt_count=None, min_vaf_for_sc_FN_class=0.01):
    if not resume:
        cells_WGS = get_cell_lists(lists_WGS, names_WGS)
        cell_CTRL = get_cell_lists(lists_CTRL, names_CTRL)
        cells_sc = get_cell_lists(lists_sc, names_sc)
        mutations_WGS = get_mutations(cells_WGS)
        mutations_CTRL = get_mutations(cell_CTRL)
        mutations_sc = get_mutations(cells_sc)
        if not os.path.exists(out_path):
            os.mkdir(out_path)
        save_mutation_array(mutations_WGS, out_path + "/WGS_intermediaries")
        save_mutation_array(mutations_CTRL, out_path + "/CTRL_intermediaries")
        save_mutation_array(mutations_sc, out_path + "/sc_intermediaries")
    else:
        mutations_WGS = load_mutation_array(out_path + "/WGS_intermediaries", names_WGS)
        mutations_CTRL = load_mutation_array(out_path + "/CTRL_intermediaries", names_CTRL)
        mutations_sc = load_mutation_array(out_path + "sc_intermediaries", names_sc)
    filtered_WGS = filter_mutations(mutations_WGS, 0.0, 0, 0.0, 0, 0.0, exonic=False, limexonic=True)
    filtered_CTRL = filter_mutations(mutations_CTRL, 0.0, 0, 0.0, 0, 10.0, exonic=False, limexonic=True)
    filtered_sc = filter_mutations(mutations_sc, min_cells_w_alt_count, min_alt_read_count, min_cell_vaf_factor, 0, min_sc_qual, max_cell_filter_=max_cells_w_alt_count, max_alt_read_filter_=max_alt_read_count, exonic=False, limexonic=True)
    filtered = filtered_sc + filtered_WGS + filtered_CTRL
    overlap_list = OverlapMutationList(filtered)
    common = common_variants(overlap_list, filtered)
    in_WGS = common[names_WGS[0]]
    for name in names_WGS[1:]:
      in_WGS = in_WGS | common[name]
    not_in_WGS = common[~(in_WGS)]
    in_CTRL = not_in_WGS[names_CTRL[0]]
    for name in names_CTRL[1:]:
      in_CTRL = in_CTRL | not_in_WGS[name]
    scu_variants = not_in_WGS[~(in_CTRL)]
    CTRL_variants = not_in_WGS[in_CTRL]
    generate_vcf(scu_variants, filtered[0], out_path + "/variants_scu.vcf")
    generate_vcf(CTRL_variants, filtered[0], out_path + "/variants_CTRL.vcf")
    scu_variants.to_csv(out_path + "/variants_scu.tsv", sep='\t', index=False)
    CTRL_variants.to_csv(out_path + "/variants_CTRL.tsv", sep='\t', index=False)
    stats = {}
    stats['scu_count'] = scu_variants.shape[0]
    return stats
    
def compare_CTRL_variants(names_WGS, names_CTRL_rna, names_CTRL_atac, path_rna, path_atac, out_path, min_alt_read_count=0, min_sc_qual=10.0):
    # note: requires process_scu to be run first to generate temp files
    if not os.path.exists(out_path):
        os.mkdir(out_path)
    mutations_WGS = load_mutation_array(path_rna + "/WGS_intermediaries", names_WGS)
    mutations_CTRL_rna = load_mutation_array(path_rna + "/CTRL_intermediaries", names_CTRL_rna)
    mutations_CTRL_atac = load_mutation_array(path_atac + "/CTRL_intermediaries", names_CTRL_atac)
    filtered_WGS = filter_mutations(mutations_WGS, 0.0, 0, 0.0, 0, 0.0, exonic=False)
    filtered_CTRL_rna = filter_mutations(mutations_CTRL_rna, 0.0, min_alt_read_count, 0.0, 0, min_sc_qual, exonic=False)
    filtered_CTRL_atac = filter_mutations(mutations_CTRL_atac, 0.0, min_alt_read_count, 0.0, 0, min_sc_qual, exonic=False)
    filtered = filtered_WGS + filtered_CTRL_rna + filtered_CTRL_atac
    overlap_list = OverlapMutationList(filtered)
    common = common_variants(overlap_list, filtered)
    in_WGS = common[names_WGS[0]]
    for name in names_WGS[1:]:
      in_WGS = in_WGS | common[name]
    yes_in_WGS = common[in_WGS]
    in_rna = yes_in_WGS[names_CTRL_rna[0]]
    for name in names_CTRL_rna[1:]:
      in_rna = in_rna | yes_in_WGS[name]
    in_atac = yes_in_WGS[names_CTRL_atac[0]]
    for name in names_CTRL_atac[1:]:
      in_atac = in_atac | yes_in_WGS[name]
    in_rna_and_atac_and_WGS = yes_in_WGS[in_rna & in_atac]
    in_rna_only_and_WGS = yes_in_WGS[in_rna & ~(in_atac)]
    in_atac_only_and_WGS = yes_in_WGS[~(in_rna) & in_atac]
    not_in_WGS = common[~(in_WGS)]
    in_rna = not_in_WGS[names_CTRL_rna[0]]
    for name in names_CTRL_rna[1:]:
      in_rna = in_rna | not_in_WGS[name]
    in_atac = not_in_WGS[names_CTRL_atac[0]]
    for name in names_CTRL_atac[1:]:
      in_atac = in_atac | not_in_WGS[name]
    in_rna_and_atac = not_in_WGS[in_rna & in_atac]
    in_rna_only = not_in_WGS[in_rna & ~(in_atac)]
    in_atac_only = not_in_WGS[~(in_rna) & in_atac]
    generate_vcf(in_rna_and_atac, filtered[0], out_path + "/variants_CTRL_both.vcf")
    generate_vcf(in_rna_only, filtered[0], out_path + "/variants_CTRL_rna.vcf")
    generate_vcf(in_atac_only, filtered[0], out_path + "/variants_CTRL_atac.vcf")
    generate_vcf(in_rna_and_atac_and_WGS, filtered[0], out_path + "/variants_WGS_both.vcf")
    generate_vcf(in_rna_only_and_WGS, filtered[0], out_path + "/variants_WGS_rna.vcf")
    generate_vcf(in_atac_only_and_WGS, filtered[0], out_path + "/variants_WGS_atac.vcf")
    in_rna_and_atac.to_csv(out_path + "/variants_CTRL_both.tsv", sep='\t', index=False)
    in_rna_only.to_csv(out_path + "/variants_CTRL_rna.tsv", sep='\t', index=False)
    in_atac_only.to_csv(out_path + "/variants_CTRL_atac.tsv", sep='\t', index=False)
    stats = {}
    stats['atac_rna_overlap'] = in_rna_and_atac.shape[0]
    stats['rna_only'] = in_rna_only.shape[0]
    stats['atac_only'] = in_atac_only.shape[0]
    stats['atac_rna_overlap_WGS'] = in_rna_and_atac_and_WGS.shape[0]
    stats['rna_only_WGS'] = in_rna_only_and_WGS.shape[0]
    stats['atac_only_WGS'] = in_atac_only_and_WGS.shape[0]
    return stats

def process_WGS_and_sc_pair(name_WGS, name_sc, out_path, _mutations_sc, _mutations_WGS, _sc_counts, bulk_sc=False, exonic=False, save_variants=False, min_cells_w_alt_count=0.0, min_alt_read_count=0, min_cell_vaf_factor=0.0, min_sc_qual=0.0, min_vaf=0.01, min_coverage=1):
    names_WGS = [name_WGS]
    names_sc = [name_sc]
    stats = {}
    mutations_sc = copy.deepcopy(_mutations_sc)
    mutations_WGS = copy.deepcopy(_mutations_WGS)
    sc_counts = copy.deepcopy(_sc_counts)
    if bulk_sc:
      filtered_sc = filter_mutations(mutations_sc, 0.0, min_alt_read_count, min_cell_vaf_factor, 0, min_sc_qual, exonic=exonic)
    else:
      filtered_sc = filter_mutations(mutations_sc, min_cells_w_alt_count, min_alt_read_count, min_cell_vaf_factor, 0, min_sc_qual, exonic=exonic)
    filtered_WGS = filter_mutations(mutations_WGS, 0.0, 0, 0.0, 0, 0.0, exonic=exonic)
    # NOTE: realized it did not make sense to keep variants that did not pass filters on sc but were found in WGS and classify these as TP, with commented out sections below, now the FN candidates include these variants
    #mutations = mutations_sc + mutations_WGS
    filtered = filtered_sc + filtered_WGS
    #overlap_list = OverlapMutationList(filtered, unfiltered=mutations).array
    overlap_obj = OverlapMutationList(filtered)
    overlap_list = overlap_obj.array    
    stats['WGS_final_count'] = overlap_obj.get_overlap(name_WGS, name_WGS)
    stats['sc_final_count'] = overlap_obj.get_overlap(name_sc, name_sc)
    stats['mutual_overlap'] = overlap_obj.get_overlap(name_WGS, name_sc)
    common = common_variants(overlap_list, filtered)
    #rank_overlap(common, out_path + "/common_variants_all.png")
    if save_variants:
      common.to_csv(out_path + "/common_variants.tsv", sep='\t', index=False)
    in_WGS = common[names_WGS[0]]
    single_cell = common[~(in_WGS)]
    in_sc = common[names_sc[0]]
    whole_genome = common[~(in_sc)]
    if save_variants:
      single_cell.to_csv(out_path + "/single_cell_variants.tsv", sep='\t', index=False)
      whole_genome.to_csv(out_path + "/whole_genome_variants.tsv", sep='\t', index=False)
      generate_vcf(common, filtered[0], out_path + "/variants_common.vcf")
      generate_vcf(single_cell, filtered[0], out_path + "/variants_single_cell.vcf")
      generate_vcf(whole_genome, filtered[0], out_path + "/variants_whole_genome.vcf")
    WGS_only_sc_counts = sc_counts.loc[sc_counts['map'].isin(whole_genome['map'])].copy()
    exonic_sc = filter_mutations(mutations_sc, 0.0, 0, 0.0, 0, 0.0, exonic=exonic)
    sc_WGS_overlap_sc_counts = sc_counts.loc[sc_counts['map'].isin(exonic_sc[0].parsed_mutations['map'])].copy()
    WGS_sc_counts = sc_counts.loc[sc_counts['map'].isin(filtered_WGS[0].parsed_mutations['map'])].copy()
    print(f"n_cells: {mutations_sc[0].n_cells}, min_cells: {min_cells_w_alt_count}")
    missed_wgs_variants, _, _, _, missing_stats = count_alleles_from_coverage(WGS_only_sc_counts, min_coverage, min_alt_read_count, min_cells_w_alt_count, min_vaf)
    mutual_variants, _, _, _, mutual_stats = count_alleles_from_coverage(sc_WGS_overlap_sc_counts, min_coverage, min_alt_read_count, min_cells_w_alt_count, min_vaf)
    wgs_variants, _, _, _, wgs_stats = count_alleles_from_coverage(WGS_sc_counts, min_coverage, min_alt_read_count, min_cells_w_alt_count, min_vaf)
    missing_stats = {'miss_' + str(key):val for key,val in missing_stats.items()}
    mutual_stats = {'tp_' + str(key):val for key,val in mutual_stats.items()}
    wgs_stats = {'ap_' + str(key):val for key,val in wgs_stats.items()}
    stats.update(mutual_stats)
    stats.update(missing_stats)
    stats.update(wgs_stats)
    if save_variants:
      WGS_only_sc_counts.to_csv(out_path + "/WGS_only_sc_counts.tsv", sep='\t', index=False)
      missed_wgs_variants.to_csv(out_path + "/WGS_missed_by_sc.tsv", sep='\t', index=False)
    #whole_genome_missed_by_sc = whole_genome[whole_genome['map'].isin(missed_wgs_variants['map'])]
    #generate_vcf(whole_genome_missed_by_sc, filtered[0], out_path + "/variants_whole_genome_missed_by_sc.vcf")
    # TODO make a new category counting the number of variants which were rejected by min alt read count filter (or VAF filter)
    #variant_stats(missed_wgs_variants, out_path + "variants_whole_genome_missed_by_sc_stats.png")
    return stats

def process_scu_stats(vartrix_scu_dir, vartrix_CTRL_dir, vartrix_WGS_dir, full_vcf_scu, full_vcf_CTRL, full_vcf_WGS, shortlist_vcf, name, out_path, resume=False, min_coverage=1, min_alt_read_count=0, min_cells_w_alt_count=0, min_vaf=0.0):
    stats = {}
    if not resume:
        sc_counts = count_vartrix(vartrix_scu_dir, full_vcf_scu)
        CTRL_counts = count_vartrix(vartrix_CTRL_dir, full_vcf_CTRL)
        if vartrix_WGS_dir:
            WGS_counts = count_vartrix(vartrix_WGS_dir, full_vcf_WGS)
        if not os.path.exists(out_path):
            os.mkdir(out_path)
        sc_counts.to_csv(out_path + "/vartrix_scu_" + name + "_counted.csv")
        CTRL_counts.to_csv(out_path + "/vartrix_CTRL_" + name + "_counted.csv")
        if vartrix_WGS_dir:
            WGS_counts.to_csv(out_path + "/vartrix_WGS_" + name + "_counted.csv")
    else:
        sc_counts = pd.read_csv(out_path + "/vartrix_scu_" + name + "_counted.csv")
        CTRL_counts = pd.read_csv(out_path + "/vartrix_CTRL_" + name + "_counted.csv")
        if vartrix_WGS_dir:
            WGS_counts = pd.read_csv(out_path + "/vartrix_WGS_" + name + "_counted.csv")
    if shortlist_vcf:
            shortlist = pd.read_csv(shortlist_vcf, sep="\t", comment="#", names=["chrom", "pos", "id", "ref", "alt", "qual", "filter", "info"])
            shortlist['map'] = shortlist['chrom'].astype(str) + ':' + shortlist['pos'].astype(int).astype(str) + ':' + shortlist['ref'] + '->' + shortlist['alt']
    stats['scu_total'] = sc_counts.shape[0]
    stats['CTRL_total'] = CTRL_counts.shape[0]
    if vartrix_WGS_dir:
        stats['WGS_total'] = WGS_counts.shape[0]
    if shortlist_vcf:
        scu_shortlisted = sc_counts[sc_counts['map'].isin(shortlist['map'])]
    else:
        scu_shortlisted = sc_counts
    _, _, _, _, scu_stats = count_alleles_from_coverage(scu_shortlisted, min_coverage, min_alt_read_count, min_cells_w_alt_count, min_vaf)
    _, _, _, _, CTRL_stats = count_alleles_from_coverage(CTRL_counts, min_coverage, min_alt_read_count, min_cells_w_alt_count, min_vaf)
    if vartrix_WGS_dir:
        _, _, _, _, WGS_stats = count_alleles_from_coverage(WGS_counts, min_coverage, min_alt_read_count, min_cells_w_alt_count, min_vaf)
    scu_stats = {'scu_' + str(key):val for key,val in scu_stats.items()}
    CTRL_stats = {'CTRL_' + str(key):val for key,val in CTRL_stats.items()}
    if vartrix_WGS_dir:
        WGS_stats = {'WGS_' + str(key):val for key,val in WGS_stats.items()}
    stats.update(scu_stats)
    stats.update(CTRL_stats)
    if vartrix_WGS_dir:
        stats.update(WGS_stats)
    return stats

def count_alleles_from_coverage(allele_counts, min_coverage, min_alt_coverage, min_cells_w_alt_count, min_vaf):
    stats = {}
    
    # obtain stats
    allele_counts['vaf'] = allele_counts['alt_count'] / (allele_counts['alt_count'] + allele_counts['ref_count'])
    allele_counts['coverage'] = allele_counts['alt_count'] + allele_counts['ref_count']
    allele_counts['depth'] = allele_counts['coverage'] / allele_counts['cell_count']
    
    # filter by min coverage
    total_variants = len(allele_counts)
    with_coverage = allele_counts[allele_counts['coverage'] >= min_coverage]
    total_no_coverage = total_variants - len(with_coverage)
    
    # filter by vaf
    above_vaf_cutoff = with_coverage[with_coverage['vaf'] >= min_vaf]
    total_above_vaf = len(above_vaf_cutoff)
    total_below_vaf = len(with_coverage) - len(above_vaf_cutoff)
    
    # filter by min alt coverage
    with_min_alt_coverage = above_vaf_cutoff[above_vaf_cutoff['alt_count'] >= min_alt_coverage]
    total_below_min_alt_coverage = len(above_vaf_cutoff) - len(with_min_alt_coverage)
    total_above_min_alt_coverage = len(with_min_alt_coverage)
    
    # filter by cell count
    with_min_cell_count = above_vaf_cutoff[above_vaf_cutoff['cells_w_alt_count'] >= min_cells_w_alt_count]
    total_below_min_cell_count = len(above_vaf_cutoff) - len(with_min_cell_count)
    total_above_min_cell_count = len(with_min_cell_count)
    
    # filter by systematic dropout
    systematic_dropout = with_coverage[with_coverage['alt_count'] == 0]
    total_systematic_above_min_coverage = len(systematic_dropout)
    print(f"Total variants: {total_variants}")
    print(f"Variants with no coverage: {total_no_coverage} ({total_no_coverage/total_variants*100}%)")
    print(f"Variants not passing vaf cutoff: {total_below_vaf} ({total_below_vaf/total_variants*100}%)")
    print(f"Variants with allele support: {total_above_vaf} ({total_above_vaf/total_variants*100}%)")
    print(f"Variants without allele support (according to min_alt_coverage): {total_below_min_alt_coverage} ({total_below_min_alt_coverage/(len(above_vaf_cutoff)+sys.float_info.epsilon)*100}%)")
    print(f"Variants with allele support (according to min_alt_coverage): {total_above_min_alt_coverage} ({total_above_min_alt_coverage/(len(above_vaf_cutoff)+sys.float_info.epsilon)*100}%)")
    print(f"Variants without allele support (according to min_cell_count): {total_below_min_cell_count} ({total_below_min_cell_count/(len(above_vaf_cutoff)+sys.float_info.epsilon)*100}%)")
    print(f"Variants with allele support (according to min_cell_count): {total_above_min_cell_count} ({total_above_min_cell_count/(len(above_vaf_cutoff)+sys.float_info.epsilon)*100}%)")
    print(f"Variants with Systematic Dropout (according to min_coverage): {total_systematic_above_min_coverage}")

    stats['below_coverage'] = total_no_coverage
    stats['above_coverage'] = total_variants
    stats['below_vaf'] = total_below_vaf
    stats['above_vaf'] = total_above_vaf
    stats['below_min_alt_cov'] = total_below_min_alt_coverage
    stats['above_min_alt_cov'] = total_above_min_alt_coverage
    stats['below_min_cell_count'] = total_below_min_cell_count
    stats['above_min_cell_count'] = total_above_min_cell_count
    return above_vaf_cutoff, with_min_alt_coverage, with_min_cell_count, systematic_dropout, stats
    
def count_overlap_vartrix_pair(vardir_rna, vcfdir_rna, vardir_atac, vcfdir_atac, vardir_unique_rna, vcfdir_unique_rna, vardir_unique_atac, vcfdir_unique_atac, asedir_WGS, out_path, min_vaf_filter = 0.0, max_vaf_filter = 1.0, alt_filter = 0, min_alt_count_WGS = 0, min_vaf_WGS = 0.0, resume=False):
    stats = {}
    if not resume:
        rna_count = count_vartrix(vardir_rna, vcfdir_rna)
        atac_count = count_vartrix(vardir_atac, vcfdir_atac)
        rna_unique_count = count_vartrix(vardir_unique_rna, vcfdir_unique_rna)
        atac_unique_count = count_vartrix(vardir_unique_atac, vcfdir_unique_atac)
        WGS_count = count_ase_reader(asedir_WGS)
        if not os.path.exists(out_path):
            os.mkdir(out_path)
        rna_count.to_csv(out_path + "/vartrix_rna_counted.csv")
        atac_count.to_csv(out_path + "/vartrix_atac_counted.csv")
        WGS_count.to_csv(out_path + "/asereader_WGS_counted.csv")
        rna_unique_count.to_csv(out_path + "/vartrix_unique_rna_counted.csv")
        atac_unique_count.to_csv(out_path + "/vartrix_unique_atac_counted.csv")
    else:
        rna_count = pd.read_csv(out_path + "/vartrix_rna_counted.csv")
        atac_count = pd.read_csv(out_path + "/vartrix_atac_counted.csv")
        rna_unique_count = pd.read_csv(out_path + "/vartrix_unique_rna_counted.csv")
        atac_unique_count = pd.read_csv(out_path + "/vartrix_unique_atac_counted.csv")
        WGS_count = pd.read_csv(out_path + "/asereader_WGS_counted.csv")
    rna_count['vaf'] = rna_count['alt_count'] / (rna_count['alt_count'] + rna_count['ref_count'])
    atac_count['vaf'] = atac_count['alt_count'] / (atac_count['alt_count'] + atac_count['ref_count'])
    rna_unique_count['vaf'] = rna_unique_count['alt_count'] / (rna_unique_count['alt_count'] + rna_unique_count['ref_count'])
    atac_unique_count['vaf'] = atac_unique_count['alt_count'] / (atac_unique_count['alt_count'] + atac_unique_count['ref_count'])
    WGS_count['vaf'] = WGS_count['alt_count'] / (WGS_count['alt_count'] + WGS_count['ref_count'])
    rna_count = rna_count[rna_count['vaf'] > min_vaf_filter]
    atac_count = atac_count[atac_count['vaf'] > min_vaf_filter]
    rna_unique_count = rna_unique_count[rna_unique_count['vaf'] > min_vaf_filter]
    atac_unique_count = atac_unique_count[atac_unique_count['vaf'] > min_vaf_filter]
    rna_count = rna_count[rna_count['vaf'] <= max_vaf_filter]
    atac_count = atac_count[atac_count['vaf'] <= max_vaf_filter]
    rna_unique_count = rna_unique_count[rna_unique_count['vaf'] <= max_vaf_filter]
    atac_unique_count = atac_unique_count[atac_unique_count['vaf'] <= max_vaf_filter]
    rna_count = rna_count[rna_count['alt_count'] > alt_filter]
    atac_count = atac_count[atac_count['alt_count'] > alt_filter]
    rna_unique_count = rna_unique_count[rna_unique_count['alt_count'] > alt_filter]
    atac_unique_count = atac_unique_count[atac_unique_count['alt_count'] > alt_filter]
    WGS_count = WGS_count[WGS_count['alt_count'] >= min_alt_count_WGS]
    WGS_count = WGS_count[WGS_count['vaf'] >= min_vaf_WGS]
    rna_count_NWGS = rna_count[~(rna_count['map'].isin(WGS_count['map']))]
    atac_count_NWGS = atac_count[~(atac_count['map'].isin(WGS_count['map']))]
    rna_unique_count_NWGS = rna_unique_count[~(rna_unique_count['map'].isin(WGS_count['map']))]
    atac_unique_count_NWGS = atac_unique_count[~(atac_unique_count['map'].isin(WGS_count['map']))]
    #rna_count_NWGS = rna_count
    #atac_count_NWGS = atac_count
    #rna_unique_count_NWGS = rna_unique_count
    #atac_unique_count_NWGS = atac_unique_count
    rna_overlap = rna_count_NWGS[rna_count_NWGS['map'].isin(atac_count_NWGS['map'])]
    atac_overlap = atac_count_NWGS[atac_count_NWGS['map'].isin(rna_count_NWGS['map'])]
    stats['rna_count'] = len(rna_count_NWGS)
    stats['atac_count'] = len(atac_count_NWGS)
    stats['overlap_count'] = len(rna_overlap)
    stats['rna_unique_count'] = len(rna_unique_count_NWGS)
    stats['atac_unique_count'] = len(atac_unique_count_NWGS)
    generate_vcf(rna_overlap, "/<path to>/header.txt", out_path + "/variants_overlap.vcf")
    return rna_overlap, atac_overlap, stats

def iterate_extract_scu_over_filter(scu_sample_dict_list, filter_dict, outfolder, filter_study_name, resume=False):
    stats_lists = {scu_sample_dict['name']:[] for scu_sample_dict in scu_sample_dict_list}
    stats = {}
    for filters in filter_dict:
      for scu_sample_dict in scu_sample_dict_list:
        stats_lists[scu_sample_dict['name']].append(extract_scu_variants(
          scu_sample_dict['lists_WGS'], 
          scu_sample_dict['names_WGS'], 
          scu_sample_dict['lists_CTRL'], 
          scu_sample_dict['names_CTRL'], 
          scu_sample_dict['lists_sc'], 
          scu_sample_dict['names_sc'], 
          outfolder + "/" + scu_sample_dict['name'] + "/", 
          resume=resume,
          min_cells_w_alt_count=filters['min_cell_count'], 
          min_alt_read_count=filters['min_alt_reads'],
          min_cell_vaf_factor=filters['min_cell_vaf_factor'], 
          min_sc_qual=filters['min_sc_qual'],
          max_alt_read_count=filters['max_alt_reads'],
          max_cells_w_alt_count=filters['max_cell_count']))
      resume = True
    for scu_sample_dict in scu_sample_dict_list:
      stats[scu_sample_dict['name']] = pd.DataFrame(stats_lists[scu_sample_dict['name']])
      stats[scu_sample_dict['name']]['min_cell_count'] = [filters['min_cell_count'] for filters in filter_dict]
      stats[scu_sample_dict['name']]['min_alt_reads'] = [filters['min_alt_reads'] for filters in filter_dict]
      stats[scu_sample_dict['name']]['min_cell_vaf_factor'] = [filters['min_cell_vaf_factor'] for filters in filter_dict]
      stats[scu_sample_dict['name']]['min_sc_qual'] = [filters['min_sc_qual'] for filters in filter_dict]
      stats[scu_sample_dict['name']]['max_alt_reads'] = [filters['max_alt_reads'] for filters in filter_dict]
      stats[scu_sample_dict['name']]['max_cell_count'] = [filters['max_cell_count'] for filters in filter_dict]
      stats[scu_sample_dict['name']].to_csv(outfolder + "/" + scu_sample_dict['name'] + "_" + filter_study_name + ".csv")

def iterate_compare_CTRL_variants_over_filter(sample_dict_list, filter_dict, outfolder, filter_study_name, resume=False):
    stats_lists = {sample_dict['name']:[] for sample_dict in sample_dict_list}
    stats = {}
    for filters in filter_dict:
      for sample_dict in sample_dict_list:
        stats_lists[sample_dict['name']].append(compare_CTRL_variants(
          sample_dict['names_WGS'], 
          sample_dict['names_CTRL_rna'], 
          sample_dict['names_CTRL_atac'], 
          sample_dict['path_rna'] + "/" + sample_dict['name'] + "/", 
          sample_dict['path_atac'] + "/" + sample_dict['name'] + "/", 
          outfolder + "/" + sample_dict['name'] + "/",
          min_alt_read_count=filters['min_alt_reads'],
          min_sc_qual=filters['min_sc_qual']))
    for sample_dict in sample_dict_list:
      stats[sample_dict['name']] = pd.DataFrame(stats_lists[sample_dict['name']])
      stats[sample_dict['name']]['min_alt_reads'] = [filters['min_alt_reads'] for filters in filter_dict]
      stats[sample_dict['name']]['min_sc_qual'] = [filters['min_sc_qual'] for filters in filter_dict]
      stats[sample_dict['name']].to_csv(outfolder + "/" + sample_dict['name'] + "_" + filter_study_name + ".csv")

def iterate_process_scu_stats(sample_dict_list, filter_dict, outfolder, filter_study_name, resume=False):
    stats_lists = {sample_dict['name']:[] for sample_dict in sample_dict_list}
    stats = {}
    for filters in filter_dict:
      for sample_dict in sample_dict_list:
        stats_lists[sample_dict['name']].append(process_scu_stats(
          sample_dict['vartrix_scu'], 
          sample_dict['vartrix_CTRL'], 
          sample_dict['vartrix_WGS'], 
          sample_dict['full_vcf_scu'], 
          sample_dict['full_vcf_CTRL'], 
          sample_dict['full_vcf_WGS'],
          sample_dict['shortlist_vcf'],
          sample_dict['name'], 
          outfolder + "/" + sample_dict['name'] + "/",
          resume=resume,
          min_coverage=filters['min_coverage'], 
          min_alt_read_count=filters['min_alt_reads'], 
          min_cells_w_alt_count=filters['min_cell_count'], 
          min_vaf=filters['min_cell_vaf_factor']))
      resume = True
    for sample_dict in sample_dict_list:
      stats[sample_dict['name']] = pd.DataFrame(stats_lists[sample_dict['name']])
      stats[sample_dict['name']]['min_cell_count'] = [filters['min_cell_count'] for filters in filter_dict]
      stats[sample_dict['name']]['min_alt_reads'] = [filters['min_alt_reads'] for filters in filter_dict]
      stats[sample_dict['name']]['min_cell_vaf_factor'] = [filters['min_cell_vaf_factor'] for filters in filter_dict]
      stats[sample_dict['name']]['min_coverage'] = [filters['min_coverage'] for filters in filter_dict]
      stats[sample_dict['name']].to_csv(outfolder + "/" + sample_dict['name'] + "_" + filter_study_name + ".csv")

def iterate_process_pair_over_filter(sample_pair_dict_list, filter_dict, outfolder, filter_study_name, resume=False, bulk_sc=False, exonic=False, save_variants=False):
    stats_lists = {sample_pair_dict['name']:[] for sample_pair_dict in sample_pair_dict_list}
    stats = {}
    for sample_pair_dict in sample_pair_dict_list:
      out_path = outfolder + "/" + sample_pair_dict['name'] + "/"
      lists_WGS = [sample_pair_dict['list_WGS']]
      names_WGS = [sample_pair_dict['name_WGS']]
      lists_sc = [sample_pair_dict['list_sc']]
      names_sc = [sample_pair_dict['name_sc']]
      stats = {}
      if not resume:
          cells_WGS = get_cell_lists(lists_WGS, names_WGS)
          cells_sc = get_cell_lists(lists_sc, names_sc)
          mutations_WGS = get_mutations(cells_WGS)
          mutations_sc = get_mutations(cells_sc)
          sc_counts = count_vartrix(sample_pair_dict['vartrix'], sample_pair_dict['full_vcf'])
          if not os.path.exists(out_path):
              os.mkdir(out_path)
          save_mutation_array(mutations_WGS, out_path + "/WGS_intermediaries")
          save_mutation_array(mutations_sc, out_path + "/sc_intermediaries")
          sc_counts.to_csv(out_path + "/vartrix_WGS_" + sample_pair_dict['name_sc'] + "_counted.csv")
      else:
          mutations_WGS = load_mutation_array(out_path + "/WGS_intermediaries", names_WGS)
          mutations_sc = load_mutation_array(out_path + "sc_intermediaries", names_sc)
          sc_counts = pd.read_csv(out_path + "/vartrix_WGS_" + sample_pair_dict['name_sc'] + "_counted.csv")
      for filters in filter_dict:
        stats_lists[sample_pair_dict['name']].append(process_WGS_and_sc_pair(
          sample_pair_dict['name_WGS'], 
          sample_pair_dict['name_sc'], 
          out_path, 
          mutations_sc,
          mutations_WGS,
          sc_counts,
          bulk_sc=bulk_sc,
          exonic=exonic,
          save_variants=save_variants,
          min_cells_w_alt_count=filters['min_cell_count'], 
          min_alt_read_count=filters['min_alt_reads'],
          min_cell_vaf_factor=filters['min_cell_vaf_factor'], 
          min_sc_qual=filters['min_sc_qual'],
          min_vaf=filters['min_vaf'],
          min_coverage=filters['min_coverage']))
      resume = True
    for sample_pair_dict in sample_pair_dict_list:
      stats[sample_pair_dict['name']] = pd.DataFrame(stats_lists[sample_pair_dict['name']])
      stats[sample_pair_dict['name']]['min_cell_count'] = [filters['min_cell_count'] for filters in filter_dict]
      stats[sample_pair_dict['name']]['min_alt_reads'] = [filters['min_alt_reads'] for filters in filter_dict]
      stats[sample_pair_dict['name']]['min_cell_vaf_factor'] = [filters['min_cell_vaf_factor'] for filters in filter_dict]
      stats[sample_pair_dict['name']]['min_sc_qual'] = [filters['min_sc_qual'] for filters in filter_dict]
      stats[sample_pair_dict['name']]['min_vaf'] = [filters['min_vaf'] for filters in filter_dict]
      stats[sample_pair_dict['name']].to_csv(outfolder + "/" + sample_pair_dict['name'] + "_" + filter_study_name + ".csv")

def iterate_count_overlap_vartrix(sample_pair_dict_list, filter_dict, outfolder, filter_study_name, resume=False):
    stats_lists = {sample_pair_dict['name']:[] for sample_pair_dict in sample_pair_dict_list}
    stats = {}
    for filters in filter_dict:
      for sample_pair_dict in sample_pair_dict_list:
        _, _, stats = count_overlap_vartrix_pair(
          sample_pair_dict['vardir_rna'], 
          sample_pair_dict['vcfdir_rna'], 
          sample_pair_dict['vardir_atac'], 
          sample_pair_dict['vcfdir_atac'],
          sample_pair_dict['vardir_unique_rna'], 
          sample_pair_dict['vcfdir_unique_rna'], 
          sample_pair_dict['vardir_unique_atac'], 
          sample_pair_dict['vcfdir_unique_atac'],
          sample_pair_dict['asedir_WGS'],
          outfolder + "/" + sample_pair_dict['name'] + "/", 
          min_vaf_filter = filters['min_vaf'],
          max_vaf_filter = filters['max_vaf'],
          alt_filter = filters['min_alt_count'],
          min_alt_count_WGS = filters['min_alt_count_WGS'], 
          min_vaf_WGS = filters['min_vaf_WGS'],
          resume=resume)
        stats_lists[sample_pair_dict['name']].append(stats)
      resume = True
    for sample_pair_dict in sample_pair_dict_list:
      stats[sample_pair_dict['name']] = pd.DataFrame(stats_lists[sample_pair_dict['name']])
      stats[sample_pair_dict['name']]['min_vaf'] = [filters['min_vaf'] for filters in filter_dict]
      stats[sample_pair_dict['name']]['max_vaf'] = [filters['max_vaf'] for filters in filter_dict]
      stats[sample_pair_dict['name']]['min_alt_count'] = [filters['min_alt_count'] for filters in filter_dict]
      stats[sample_pair_dict['name']].to_csv(outfolder + "/" + sample_pair_dict['name'] + "_" + filter_study_name + ".csv")

def intersect_scu(tsv1, tsv2, outfile, vcfheader=None):
    list1 = pd.read_csv(tsv1, sep='\t')
    list2 = pd.read_csv(tsv2, sep='\t')
    isec = list1[list1['map'].isin(list2['map'])]
    isec.to_csv(outfile, sep='\t', index=False)
    if vcfheader:
      generate_vcf(isec, vcfheader, outfile.split('.tsv')[0] + '.vcf')

def process_scu():
    resume = True
    outfolder = "/<path to>/hapcaller_sc_rna_max_all_exonic/"
    min_cell_count_base = 0.0
    min_cell_vaf_factor_base = 0.0 
    min_sc_qual_base = 0.0
    min_alt_reads_base = 0
    max_cell_count_base = None
    max_alt_reads_base = None
    min_cell_counts = [0.0, 1.0, 5.0, 10.0, 50.0, 100.0, 200.0]
    min_sc_quals = [0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 20.0, 50.0, 100.0, 200.0]
    #min_cell_vaf_factors = [0.0, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 0.5, 0.9]
    if not os.path.exists(outfolder):
        os.mkdir(outfolder)
    
    list_WGS_EPSRC1 = "/<path to>/WGS/EPSRC1dna/cell_list_consensus.txt"
    list_WGS_EPSRC2 = "/<path to>/WGS/EPSRC2dna/cell_list_consensus.txt"
    list_WGS_EPSRC4 = "/<path to>/WGS/EPSRC4dna/cell_list_consensus.txt"
    list_WGS_CTRL = [list_WGS_EPSRC1, list_WGS_EPSRC2, list_WGS_EPSRC4]
    name_WGS_EPSRC1 = 'EPSRC1dna'
    name_WGS_EPSRC2 = 'EPSRC2dna'
    name_WGS_EPSRC4 = 'EPSRC4dna'
    name_WGS_CTRL = [name_WGS_EPSRC1, name_WGS_EPSRC2, name_WGS_EPSRC4]
    
    list_str_EPSRC1 = "/<path to>/EPSRC_splits/EpSRCrna1/cell_list_bulk_strelka.txt"
    list_str_EPSRC2 = "/<path to>/EPSRC_splits/EpSRCrna2/cell_list_bulk_strelka.txt"
    list_str_EPSRC4 = "/<path to>/EPSRC_splits/EpSRCrnac4/cell_list_bulk_strelka.txt"
    list_str_CTRL = "/<path to>/EPSRC_splits/EpSRCrnaCTRL/cell_list_bulk_strelka.txt"
    name_str_EPSRC1 = 'EPSRC1rna_strelka'
    name_str_EPSRC2 = 'EPSRC2rna_strelka'
    name_str_EPSRC4 = 'EPSRC4rna_strelka'
    name_str_CTRL = 'CTRLrna_strelka'
    
    list_sam_EPSRC1 = "/<path to>/EPSRC_splits/EpSRCrna1/cell_list_bulk_samtools.txt"
    list_sam_EPSRC2 = "/<path to>/EPSRC_splits/EpSRCrna2/cell_list_bulk_samtools.txt"
    list_sam_EPSRC4 = "/<path to>/EPSRC_splits/EpSRCrna4/cell_list_bulk_samtools.txt"
    list_sam_CTRL = "/<path to>/EPSRC_splits/EpSRCrnaCTRL/cell_list_bulk_samtools.txt"
    name_sam_EPSRC1 = 'EPSRC1rna_samtools'
    name_sam_EPSRC2 = 'EPSRC2rna_samtools'
    name_sam_EPSRC4 = 'EPSRC4rna_samtools'
    name_sam_CTRL = 'CTRLrna_samtools'
    
    list_hap_EPSRC1 = "/<path to>/EPSRC_splits/EpSRCrna1/cell_list_bulk_hapcaller.txt"
    list_hap_EPSRC2 = "/<path to>/EPSRC_splits/EpSRCrna2/cell_list_bulk_hapcaller.txt"
    list_hap_EPSRC4 = "/<path to>/EPSRC_splits/EpSRCrna4/cell_list_bulk_hapcaller.txt"
    list_hap_CTRL = "/<path to>/EPSRC_splits/EpSRCrnaCTRL/cell_list_bulk_hapcaller.txt"
    name_hap_EPSRC1 = 'EPSRC1rna_hapcaller'
    name_hap_EPSRC2 = 'EPSRC2rna_hapcaller'
    name_hap_EPSRC4 = 'EPSRC4rna_hapcaller'
    name_hap_CTRL = 'CTRLrna_hapcaller'
    
    list_sc_EPSRC1 = "/<path to>/EPSRC_splits/EpSRCrna1/cell_list_hapcaller.txt"
    list_sc_EPSRC2 = "/<path to>/EPSRC_splits/EpSRCrna2/cell_list_hapcaller.txt"
    list_sc_EPSRC4 = "/<path to>/EPSRC_splits/EpSRCrna4/cell_list_hapcaller.txt"
    list_sc_CTRL = "/<path to>/EPSRC_splits/EpSRCrnaCTRL/cell_list_hapcaller.txt"
    name_sc_EPSRC1 = 'EPSRC1rna'
    name_sc_EPSRC2 = 'EPSRC2rna'
    name_sc_EPSRC4 = 'EPSRC4rna'
    name_sc_CTRL = 'CTRLrna'
    
    EPSRC1_scu_sample_dict = {'name':'EPSRC1','lists_WGS':[list_WGS_EPSRC1],'names_WGS':[name_WGS_EPSRC1],'lists_CTRL':[list_str_EPSRC1,list_sam_EPSRC1,list_hap_EPSRC1],'names_CTRL':[name_str_EPSRC1,name_sam_EPSRC1,name_hap_EPSRC1],'lists_sc':[list_sc_EPSRC1],'names_sc':[name_sc_EPSRC1]}
    EPSRC2_scu_sample_dict = {'name':'EPSRC2','lists_WGS':[list_WGS_EPSRC2],'names_WGS':[name_WGS_EPSRC2],'lists_CTRL':[list_str_EPSRC2,list_sam_EPSRC2,list_hap_EPSRC2],'names_CTRL':[name_str_EPSRC2,name_sam_EPSRC2,name_hap_EPSRC2],'lists_sc':[list_sc_EPSRC2],'names_sc':[name_sc_EPSRC2]}
    EPSRC4_scu_sample_dict = {'name':'EPSRC4','lists_WGS':[list_WGS_EPSRC4],'names_WGS':[name_WGS_EPSRC4],'lists_CTRL':[list_str_EPSRC4,list_sam_EPSRC4,list_hap_EPSRC4],'names_CTRL':[name_str_EPSRC4,name_sam_EPSRC4,name_hap_EPSRC4],'lists_sc':[list_sc_EPSRC4],'names_sc':[name_sc_EPSRC4]}
    CTRL_scu_sample_dict = {'name':'CTRL','lists_WGS':list_WGS_CTRL,'names_WGS':name_WGS_CTRL,'lists_CTRL':[list_str_CTRL,list_sam_CTRL,list_hap_CTRL],'names_CTRL':[name_str_CTRL,name_sam_CTRL,name_hap_CTRL],'lists_sc':[list_sc_CTRL],'names_sc':[name_sc_CTRL]}
    scu_sample_dict_list = [EPSRC1_scu_sample_dict, EPSRC2_scu_sample_dict, EPSRC4_scu_sample_dict, CTRL_scu_sample_dict]
    
    min_cell_count_filter_dict = [{'min_cell_count':i,'min_alt_reads':min_alt_reads_base,'min_cell_vaf_factor':min_cell_vaf_factor_base,'min_sc_qual':min_sc_qual_base,'max_alt_reads':max_alt_reads_base,'max_cell_count':max_cell_count_base} for i in min_cell_counts]
    standard_output_dict = [{'min_cell_count':0.0,'min_alt_reads':0,'min_cell_vaf_factor':0.0,'min_sc_qual':100.0,'max_alt_reads':None,'max_cell_count':None}]
    min_sc_qual_filter_dict = [{'min_cell_count':min_cell_count_base,'min_alt_reads':min_alt_reads_base,'min_cell_vaf_factor':min_cell_vaf_factor_base,'min_sc_qual':i,'max_alt_reads':max_alt_reads_base,'max_cell_count':max_cell_count_base} for i in min_sc_quals]
    #min_cell_vaf_factor_dict = [{'min_cell_count':min_cell_count_base,'min_alt_reads':min_alt_reads_base,'min_cell_vaf_factor':i,'min_sc_qual':100.0,'max_alt_reads':max_alt_reads_base,'max_cell_count':max_cell_count_base} for i in min_cell_vaf_factors]
    
    #iterate_extract_scu_over_filter(scu_sample_dict_list, min_cell_count_filter_dict, outfolder, 'cell_count_filter_stats', resume=resume)
    #iterate_extract_scu_over_filter(scu_sample_dict_list, min_sc_qual_filter_dict, outfolder, 'min_sc_qual_filter_stats', resume=True)
    #iterate_extract_scu_over_filter(scu_sample_dict_list, min_cell_vaf_factor_dict, outfolder, 'min_cell_vaf_factors_stats_100', resume=True)
    iterate_extract_scu_over_filter(scu_sample_dict_list, standard_output_dict, outfolder, 'standard_stats', resume=True)

def scu_stats():
    resume = False
    outfolder = "/<path to>/scu_results/scu_atac_40_100_stats"
    min_cell_count_base = 0.0
    min_cell_vaf_factor_base = 0.0 # changed this, probably should be 0
    min_alt_reads_base = 0
    min_coverage_base = 1
    min_coverages = [1,2,4,10,20,40,100,200,400,1000,2000,4000,10000]
    min_cell_vaf_factors = [0.0, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 0.5, 0.9]
    if not os.path.exists(outfolder):
        os.mkdir(outfolder)
    
    vartrix_WGS_EPSRC1 = "/<path to>/vartrix_WGS_sc_counts_consensus/EPSRC1atac_all_WGS_in_sc_counts/"
    full_vcf_WGS_EPSRC1 = "/<path to>/WGS/EPSRC1_variants/EPSRC1_consensus/consensus3_atac.vcf.gz"
    vartrix_WGS_EPSRC2 = "/<path to>/vartrix_WGS_sc_counts_consensus/EPSRC2atac_all_WGS_in_sc_counts/"
    full_vcf_WGS_EPSRC2 = "/<path to>/EPSRC2_variants/EPSRC2_consensus/consensus3_atac.vcf.gz"
    vartrix_WGS_EPSRC4 = "/<path to>/vartrix_WGS_sc_counts_consensus/EPSRC4atac_all_WGS_in_sc_counts/"
    full_vcf_WGS_EPSRC4 = "/<path to>/WGS/EPSRC4_variants/EPSRC4_consensus/consensus3_atac.vcf.gz"
    vartrix_WGS_CTRL = None
    full_vcf_WGS_CTRL = None
    
    vartrix_CTRL_EPSRC1 = "/<path to>/scu_results/vartrix_counts_CTRL/EPSRC1atac/"
    full_vcf_CTRL_EPSRC1 = "/<path to>/scu_results/CTRL_variants/EPSRC1atac_CTRL.vcf"
    vartrix_CTRL_EPSRC2 = "/<path to>/scu_results/vartrix_counts_CTRL/EPSRC2atac/"
    full_vcf_CTRL_EPSRC2 = "/<path to>/scu_results/CTRL_variants/EPSRC2atac_CTRL.vcf"
    vartrix_CTRL_EPSRC4 = "/<path to>/scu_results/vartrix_counts_CTRL/EPSRC4atac/"
    full_vcf_CTRL_EPSRC4 = "/<path to>/scu_results/CTRL_variants/EPSRC4atac_CTRL.vcf"
    vartrix_CTRL_CTRL = "/<path to>/scu_results/vartrix_counts_CTRL/CTRLatac/"
    full_vcf_CTRL_CTRL = "/<path to>/scu_results/CTRL_variants/CTRLatac_CTRL.vcf"
    
    vartrix_scu_EPSRC1 = "/<path to>/scu_results/vartrix_counts_scu_40_100/EPSRC1atac/"
    full_vcf_scu_EPSRC1 = "/<path to>/scu_results/EPSRC1atac_scu_isec_40_100.vcf"
    vartrix_scu_EPSRC2 = "/<path to>/scu_results/vartrix_counts_scu_40_100/EPSRC2atac/"
    full_vcf_scu_EPSRC2 = "/<path to>/scu_results/EPSRC2atac_scu_isec_40_100.vcf"
    vartrix_scu_EPSRC4 = "/<path to>/scu_results/vartrix_counts_scu_40_100/EPSRC4atac/"
    full_vcf_scu_EPSRC4 = "/<path to>/scu_results/EPSRC4atac_scu_isec_40_100.vcf"
    vartrix_scu_CTRL = "/<path to>/scu_results/vartrix_counts_scu_40_100/CTRLatac/"
    full_vcf_scu_CTRL = "/<path to>/scu_results/CTRLatac_scu_isec_40_100.vcf"
    
    shortlist_vcf_EPSRC1 = None
    shortlist_vcf_EPSRC2 = None
    shortlist_vcf_EPSRC4 = None
    shortlist_vcf_CTRL = None
    
    
    EPSRC1_dict = {'name':'EPSRC1','vartrix_WGS':vartrix_WGS_EPSRC1,'full_vcf_WGS':full_vcf_WGS_EPSRC1,'vartrix_CTRL':vartrix_CTRL_EPSRC1,'full_vcf_CTRL':full_vcf_CTRL_EPSRC1,'vartrix_scu':vartrix_scu_EPSRC1,'full_vcf_scu':full_vcf_scu_EPSRC1,'shortlist_vcf':shortlist_vcf_EPSRC1}
    EPSRC2_dict = {'name':'EPSRC2','vartrix_WGS':vartrix_WGS_EPSRC2,'full_vcf_WGS':full_vcf_WGS_EPSRC2,'vartrix_CTRL':vartrix_CTRL_EPSRC2,'full_vcf_CTRL':full_vcf_CTRL_EPSRC2,'vartrix_scu':vartrix_scu_EPSRC2,'full_vcf_scu':full_vcf_scu_EPSRC2,'shortlist_vcf':shortlist_vcf_EPSRC2}
    EPSRC4_dict = {'name':'EPSRC4','vartrix_WGS':vartrix_WGS_EPSRC4,'full_vcf_WGS':full_vcf_WGS_EPSRC4,'vartrix_CTRL':vartrix_CTRL_EPSRC4,'full_vcf_CTRL':full_vcf_CTRL_EPSRC4,'vartrix_scu':vartrix_scu_EPSRC4,'full_vcf_scu':full_vcf_scu_EPSRC4,'shortlist_vcf':shortlist_vcf_EPSRC4}
    CTRL_dict = {'name':'CTRL','vartrix_WGS':vartrix_WGS_CTRL,'full_vcf_WGS':full_vcf_WGS_CTRL,'vartrix_CTRL':vartrix_CTRL_CTRL,'full_vcf_CTRL':full_vcf_CTRL_CTRL,'vartrix_scu':vartrix_scu_CTRL,'full_vcf_scu':full_vcf_scu_CTRL,'shortlist_vcf':shortlist_vcf_CTRL}
    sample_dict_list = [EPSRC1_dict, EPSRC2_dict, EPSRC4_dict, CTRL_dict]
    
    min_coverage_filter_dict = [{'min_cell_count':min_cell_count_base,'min_alt_reads':min_alt_reads_base,'min_cell_vaf_factor':min_cell_vaf_factor_base,'min_coverage':i} for i in min_coverages]
    min_cell_vaf_factor_filter_dict = [{'min_cell_count':min_cell_count_base,'min_alt_reads':min_alt_reads_base,'min_cell_vaf_factor':i,'min_coverage':min_coverage_base} for i in min_cell_vaf_factors]
    
    iterate_process_scu_stats(sample_dict_list, min_coverage_filter_dict, outfolder, 'min_coverage_stats', resume=resume)
    iterate_process_scu_stats(sample_dict_list, min_cell_vaf_factor_filter_dict, outfolder, 'min_vaf_stats', resume=True)

"""
sample_dict: a list of samples containing according WGS and single-cell variant calling information, which include:
  "name": overall sample name i.e. "EPSRC1"
  "list_WGS": cell list for WGS variants
  "name_WGS": name for WGS dataset, i.e. "EPSRC1dna"
  "list_sc": cell list for single-cell variants
  "name_sc": name for single-cell dataset, i.e. "EPSRC1rna"
  "vartrix":"/<path to>/vartrix_WGS_sc_counts_consensus/EPSRC4rna_all_WGS_in_sc_counts/",
  "full_vcf":"/<path to>/WGS/EPSRC4dna/variants_consensus/pseudo_cell/consensus3.vcf.gz"
"""
def process_WGS_v_sc(sample_dict_path: str, filter_dict_path: str, test_name: str, resume: bool):
    print(f"process_WGS_v_sc: {locals()}")

    with open(sample_dict_path) as f:
      sample_json = json.load(f)

    if not os.path.exists(sample_json["working"]):
      os.mkdir(sample_json["working"])

    with open(filter_dict_path) as f:
      filter_json = json.load(f)
    
    iterables = filter_json["iterables"]
    
    for it_name, it_list in iterables.items():
      valid_filters = ["min_cell_count", "min_alt_reads", "min_cell_vaf_factor", "min_sc_qual", "min_vaf", "min_coverage"]
      if not it_name in valid_filters:
        raise ValueError(f"Invalid filter iterable: {it_name}")

      if not set(filter_json["base"].keys()) == set(valid_filters):
        raise ValueError(f"Invalid base filters: {filter_json['base'].keys()}; should be {valid_filters}")

      filter_dict_list = []
      for it in it_list:
        filter_dict = filter_json["base"].copy()
        filter_dict[it_name] = it
        filter_dict_list.append(filter_dict)
        print(f"running filter: {filter_dict}")
      
      iterate_process_pair_over_filter(sample_json["samples"], filter_dict_list, sample_json["working"], test_name + "_varying_" + it_name, resume=resume, bulk_sc=filter_json["is_bulk"], exonic=filter_json["exonic"], save_variants=filter_json["save_variants"])

def process_compare_CTRL_variants():
    outfolder = "/<path to>/compare_rna_atac_results/CTRL_overlap_all/"
    min_sc_qual_base = 0.0
    min_alt_reads_base = 0
    min_alt_reads = [1,5,10,20,30,50]
    
    if not os.path.exists(outfolder):
        os.mkdir(outfolder)
    
    name_WGS_EPSRC1 = ['EPSRC1dna']
    name_WGS_EPSRC2 = ['EPSRC2dna']
    name_WGS_EPSRC4 = ['EPSRC4dna']
    names_CTRL_EPSRC1_rna = ['EPSRC1rna_hapcaller','EPSRC1rna_samtools','EPSRC1rna_strelka']
    names_CTRL_EPSRC2_rna = ['EPSRC2rna_hapcaller','EPSRC2rna_samtools','EPSRC2rna_strelka']
    names_CTRL_EPSRC4_rna = ['EPSRC4rna_hapcaller','EPSRC4rna_samtools','EPSRC4rna_strelka']
    path_rna = "/<path to>/scu_results/strelka_sc_max_all/"
    names_CTRL_EPSRC1_atac = ['EPSRC1atac_hapcaller','EPSRC1atac_samtools','EPSRC1atac_strelka']
    names_CTRL_EPSRC2_atac = ['EPSRC2atac_hapcaller','EPSRC2atac_samtools','EPSRC2atac_strelka']
    names_CTRL_EPSRC4_atac = ['EPSRC4atac_hapcaller','EPSRC4atac_samtools','EPSRC4atac_strelka']
    path_atac = "/<path to>/scu_results/strelka_sc_atac_max_all/"
    
    EPSRC1_dict = {'name':'EPSRC1',
                   'names_WGS':name_WGS_EPSRC1,
                   'names_CTRL_rna':names_CTRL_EPSRC1_rna,
                   'names_CTRL_atac':names_CTRL_EPSRC1_atac,
                   'path_rna':path_rna,
                   'path_atac':path_atac}
    EPSRC2_dict = {'name':'EPSRC2',
                   'names_WGS':name_WGS_EPSRC2,
                   'names_CTRL_rna':names_CTRL_EPSRC2_rna,
                   'names_CTRL_atac':names_CTRL_EPSRC2_atac,
                   'path_rna':path_rna,
                   'path_atac':path_atac}
    EPSRC4_dict = {'name':'EPSRC4',
                   'names_WGS':name_WGS_EPSRC4,
                   'names_CTRL_rna':names_CTRL_EPSRC4_rna,
                   'names_CTRL_atac':names_CTRL_EPSRC4_atac,
                   'path_rna':path_rna,
                   'path_atac':path_atac}
    sample_dict_list = [EPSRC1_dict, EPSRC2_dict, EPSRC4_dict]
    
    min_alt_read_count_filter_dict = [{'min_alt_reads':i,'min_sc_qual':min_sc_qual_base} for i in min_alt_reads]
    standard_filter_dict = [{'min_alt_reads':min_alt_reads_base,'min_sc_qual':min_sc_qual_base}]
    
    iterate_compare_CTRL_variants_over_filter(sample_dict_list, min_alt_read_count_filter_dict, outfolder, 'min_alt_reads_filter_stats')
    iterate_compare_CTRL_variants_over_filter(sample_dict_list, standard_filter_dict, outfolder, 'standard_filter_stats')

def count_overlap_CTRL():
    outfolder = "/<path to>/compare_rna_atac_results/test"
    min_vaf_base = 0.005
    max_vaf_base = 1.0
    min_vaf_WGS_base = 0.0
    min_alt_count_WGS_base = 1
    min_vafs = [0.0, 0.0005, 0.001, 0.005, 0.01, 0.05]
    max_vafs = [1.0, 0.5, 0.4, 0.3, 0.2, 0.1, 0.05]
    min_alt_count_base = 5
    
    if not os.path.exists(outfolder):
        os.mkdir(outfolder)
    
    vardir_rna_EPSRC1 = "/<path to>/compare_rna_atac_results/CTRL_overlap_all/vartrix_rna/EPSRC1/"
    vardir_rna_EPSRC2 = "/<path to>/compare_rna_atac_results/CTRL_overlap_all/vartrix_rna/EPSRC2/"
    vardir_rna_EPSRC4 = "/<path to>/compare_rna_atac_results/CTRL_overlap_all/vartrix_rna/EPSRC4/"
    vcfdir_rna_EPSRC1 = "/<path to>/compare_rna_atac_results/CTRL_overlap_all/vartrix_rna/EPSRC1/variants_CTRL_both.vcf"
    vcfdir_rna_EPSRC2 = "/<path to>/compare_rna_atac_results/CTRL_overlap_all/vartrix_rna/EPSRC2/variants_CTRL_both.vcf"
    vcfdir_rna_EPSRC4 = "/<path to>/compare_rna_atac_results/CTRL_overlap_all/vartrix_rna/EPSRC4/variants_CTRL_both.vcf"
    
    vardir_atac_EPSRC1 = "/<path to>/compare_rna_atac_results/CTRL_overlap_all/vartrix_atac/EPSRC1/"
    vardir_atac_EPSRC2 = "/<path to>/compare_rna_atac_results/CTRL_overlap_all/vartrix_atac/EPSRC2/"
    vardir_atac_EPSRC4 = "/<path to>/compare_rna_atac_results/CTRL_overlap_all/vartrix_atac/EPSRC4/"
    vcfdir_atac_EPSRC1 = "/<path to>/compare_rna_atac_results/CTRL_overlap_all/vartrix_atac/EPSRC1/variants_CTRL_both.vcf"
    vcfdir_atac_EPSRC2 = "/<path to>/compare_rna_atac_results/CTRL_overlap_all/vartrix_atac/EPSRC2/variants_CTRL_both.vcf"
    vcfdir_atac_EPSRC4 = "/<path to>/compare_rna_atac_results/CTRL_overlap_all/vartrix_atac/EPSRC4/variants_CTRL_both.vcf"
    
    vardir_unique_rna_EPSRC1 = "/<path to>/compare_rna_atac_results/CTRL_overlap_all/vartrix_unique_rna/EPSRC1/"
    vardir_unique_rna_EPSRC2 = "/<path to>/compare_rna_atac_results/CTRL_overlap_all/vartrix_unique_rna/EPSRC2/"
    vardir_unique_rna_EPSRC4 = "/<path to>/compare_rna_atac_results/CTRL_overlap_all/vartrix_unique_rna/EPSRC4/"
    vcfdir_unique_rna_EPSRC1 = "/<path to>/compare_rna_atac_results/CTRL_overlap_all/vartrix_unique_rna/EPSRC1/variants_CTRL_rna.vcf"
    vcfdir_unique_rna_EPSRC2 = "/<path to>/compare_rna_atac_results/CTRL_overlap_all/vartrix_unique_rna/EPSRC2/variants_CTRL_rna.vcf"
    vcfdir_unique_rna_EPSRC4 = "/<path to>/compare_rna_atac_results/CTRL_overlap_all/vartrix_unique_rna/EPSRC4/variants_CTRL_rna.vcf"
    
    vardir_unique_atac_EPSRC1 = "/<path to>/compare_rna_atac_results/CTRL_overlap_all/vartrix_unique_atac/EPSRC1/"
    vardir_unique_atac_EPSRC2 = "/<path to>/compare_rna_atac_results/CTRL_overlap_all/vartrix_unique_atac/EPSRC2/"
    vardir_unique_atac_EPSRC4 = "/<path to>/compare_rna_atac_results/CTRL_overlap_all/vartrix_unique_atac/EPSRC4/"
    vcfdir_unique_atac_EPSRC1 = "/<path to>/compare_rna_atac_results/CTRL_overlap_all/vartrix_unique_atac/EPSRC1/variants_CTRL_atac.vcf"
    vcfdir_unique_atac_EPSRC2 = "/<path to>/compare_rna_atac_results/CTRL_overlap_all/vartrix_unique_atac/EPSRC2/variants_CTRL_atac.vcf"
    vcfdir_unique_atac_EPSRC4 = "/<path to>/compare_rna_atac_results/CTRL_overlap_all/vartrix_unique_atac/EPSRC4/variants_CTRL_atac.vcf"
    
    asedir_WGS_EPSRC1 = "/<path to>/compare_rna_atac_results/CTRL_overlap_all/vartrix_WGS/EPSRC1/counts.table"
    asedir_WGS_EPSRC2 = "/<path to>/compare_rna_atac_results/CTRL_overlap_all/vartrix_WGS/EPSRC2/counts.table"
    asedir_WGS_EPSRC4 = "/<path to>/compare_rna_atac_results/CTRL_overlap_all/vartrix_WGS/EPSRC4/counts.table"
    
    EPSRC1_dict = {'name':'EPSRC1',
                   'vardir_rna':vardir_rna_EPSRC1,
                   'vcfdir_rna':vcfdir_rna_EPSRC1,
                   'vardir_atac':vardir_atac_EPSRC1,
                   'vcfdir_atac':vcfdir_atac_EPSRC1,
                   'vardir_unique_rna':vardir_unique_rna_EPSRC1,
                   'vcfdir_unique_rna':vcfdir_unique_rna_EPSRC1,
                   'vardir_unique_atac':vardir_unique_atac_EPSRC1,
                   'vcfdir_unique_atac':vcfdir_unique_atac_EPSRC1,
                   'asedir_WGS':asedir_WGS_EPSRC1}
    EPSRC2_dict = {'name':'EPSRC2',
                   'vardir_rna':vardir_rna_EPSRC2,
                   'vcfdir_rna':vcfdir_rna_EPSRC2,
                   'vardir_atac':vardir_atac_EPSRC2,
                   'vcfdir_atac':vcfdir_atac_EPSRC2,
                   'vardir_unique_rna':vardir_unique_rna_EPSRC2,
                   'vcfdir_unique_rna':vcfdir_unique_rna_EPSRC2,
                   'vardir_unique_atac':vardir_unique_atac_EPSRC2,
                   'vcfdir_unique_atac':vcfdir_unique_atac_EPSRC2,
                   'asedir_WGS':asedir_WGS_EPSRC2}
    EPSRC4_dict = {'name':'EPSRC4',
                   'vardir_rna':vardir_rna_EPSRC4,
                   'vcfdir_rna':vcfdir_rna_EPSRC4,
                   'vardir_atac':vardir_atac_EPSRC4,
                   'vcfdir_atac':vcfdir_atac_EPSRC4,
                   'vardir_unique_rna':vardir_unique_rna_EPSRC4,
                   'vcfdir_unique_rna':vcfdir_unique_rna_EPSRC4,
                   'vardir_unique_atac':vardir_unique_atac_EPSRC4,
                   'vcfdir_unique_atac':vcfdir_unique_atac_EPSRC4,
                   'asedir_WGS':asedir_WGS_EPSRC4}
    sample_dict_list = [EPSRC1_dict, EPSRC2_dict, EPSRC4_dict]
    
    min_vaf_filter_dict = [{'min_alt_count':min_alt_count_base,'min_vaf':i, 'max_vaf':max_vaf_base, 'min_alt_count_WGS':min_alt_count_WGS_base, 'min_vaf_WGS':min_vaf_WGS_base} for i in min_vafs]
    max_vaf_filter_dict = [{'min_alt_count':min_alt_count_base,'min_vaf':min_vaf_base, 'max_vaf':i, 'min_alt_count_WGS':min_alt_count_WGS_base, 'min_vaf_WGS':min_vaf_WGS_base} for i in max_vafs]
    standard_filter_dict = [{'min_alt_count':5,'min_vaf':min_vaf_base,'max_vaf':max_vaf_base, 'min_alt_count_WGS':min_alt_count_WGS_base, 'min_vaf_WGS':min_vaf_WGS_base}]
    
    #iterate_count_overlap_vartrix(sample_dict_list, min_vaf_filter_dict, outfolder, 'min_vaf_filter', resume=True)
    #iterate_count_overlap_vartrix(sample_dict_list, max_vaf_filter_dict, outfolder, 'max_vaf_filter', resume=True)
    iterate_count_overlap_vartrix(sample_dict_list, standard_filter_dict, outfolder, 'standard_filters', resume=False)

def parse_args():
    parser = argparse.ArgumentParser(description='Variant statistics analysis')
    subparsers = parser.add_subparsers(dest='func')

    pair_parser = subparsers.add_parser('process_pair')
    pair_parser.add_argument('samples', help='path to sample json file')
    pair_parser.add_argument('filters', help='path to filter json file')
    pair_parser.add_argument('test_name', help='name of test (should be unique)')
    pair_parser.add_argument('-r', '--resume', help='if specified, attempt to resume', action="store_true")

    aggvtx_parser = subparsers.add_parser('aggregate_vartrix')
    aggvtx_parser.add_argument('vartrix_dir', help='path to vartrix output folder')
    aggvtx_parser.add_argument('vcf', help='path to variant vcf file used to generate vartrix output')
    aggvtx_parser.add_argument('outcsv', help='path to resulting .csv file')
    args = parser.parse_args()

    return args

def main(): 
    ### Outdated approaches for processing overlaps

    #process_WGS_v_scRNA()
    #process_WGS_v_scATAC()
    #process_WGS_v_bulk()
    #process_scu()
    #scu_stats()
    #process_compare_CTRL_variants()
    #count_overlap_CTRL()
    
    #intersect_scu("/<path to>/scu_results/strelka_sc_atac_max_all/EPSRC1/variants_scu_50.tsv", "/<path to>/scu_results/hapcaller_sc_atac_max_all/EPSRC1/variants_scu_80.tsv", "/<path to>/scu_results/EPSRC1atac_scu_isec_50_80.tsv", vcfheader="/<path to>/scu_results/strelka_sc_atac_max_all/EPSRC1/sc_intermediaries/EPSRC1atac/header.txt")
    #intersect_scu("/<path to>/scu_results/strelka_sc_atac_max_all/EPSRC2/variants_scu_50.tsv", "/<path to>/scu_results/hapcaller_sc_atac_max_all/EPSRC2/variants_scu_80.tsv", "/<path to>/scu_results/EPSRC2atac_scu_isec_50_80.tsv", vcfheader="/<path to>/scu_results/strelka_sc_atac_max_all/EPSRC2/sc_intermediaries/EPSRC2atac/header.txt")
    #intersect_scu("/<path to>/scu_results/strelka_sc_atac_max_all/EPSRC4/variants_scu_50.tsv", "/<path to>/scu_results/hapcaller_sc_atac_max_all/EPSRC4/variants_scu_80.tsv", "/<path to>/scu_results/EPSRC4atac_scu_isec_50_80.tsv", vcfheader="/<path to>/scu_results/strelka_sc_atac_max_all/EPSRC4/sc_intermediaries/EPSRC4atac/header.txt")
    #intersect_scu("/<path to>/scu_results/strelka_sc_atac_max_all/CTRL/variants_scu_50.tsv", "/<path to>/scu_results/hapcaller_sc_atac_max_all/CTRL/variants_scu_80.tsv", "/<path to>/scu_results/CTRLatac_scu_isec_50_80.tsv", vcfheader="/<path to>/scu_results/strelka_sc_atac_max_all/CTRL/sc_intermediaries/CTRLatac/header.txt")
    
    ### Approaches for Statistics on Single VCF ###
    
    #mylist = [{'name':'EPSRC1','file':"/<path to>/WGS/EPSRC1dna/strelka_somatic/CTR1_somatic/somatic_annotated_pass_no_dbsnp.TUMOR.vcf"},
    #          {'name':'EPSRC2','file':"/<path to>/WGS/EPSRC2dna/strelka_somatic/CTR1_somatic/somatic_annotated_pass_no_dbsnp.TUMOR.vcf"},
    #          {'name':'EPSRC4','file':"/<path to>/WGS/EPSRC4dna/strelka_somatic/CTR1_somatic/somatic_annotated_pass_no_dbsnp.TUMOR.vcf"}]
    #filter_and_overlap_multiple(mylist, "/<path to>/WGS/somatic_overlap_no_dbsnp.tsv")
    
    #filter_single("/<path to>/WGS/EPSRC1dna/strelka_variants_somatic_candidates_annotated_no_dbsnp.TUMOR.vcf", "/<path to>/WGS/EPSRC1dna/strelka_variants_somatic_candidates_annotated_no_dbsnp.tsv")
    
    #table = pd.read_csv("/<path to>/scu_results/scu_atac_40_100_stats/EPSRC1/vartrix_WGS_EPSRC1_counted.csv")
    #table = table[table['alt_count'] >= 1]
    #generate_vcf(table, "/<path to>/scu_results/strelka_sc_atac_max_all/EPSRC1/WGS_intermediaries/EPSRC1dna/header.txt", "/<path to>/scu_results/WGS_atac_sigprofiler/EPSRC1_WGS.vcf", atac=False)
    
    #sc_counts = count_vartrix("/<path to>/scu_results/vartrix_counts_CTRL/CTRLrna/", "/<path to>/PTM_overlap/CTRL/CTRL.vcf")
    #sc_counts.to_csv("/<path to>/PTM_overlap/CTRL/CTRL_counted.csv")

    args = parse_args()
    if args.func == 'process_pair':
      process_WGS_v_sc(args.samples, args.filters, args.test_name, resume=args.resume)
    elif args.func == 'aggregate_vartrix':
      counts = count_vartrix(args.vartrix_dir, args.vcf)
      total = len(counts)
      counts['depth'] = counts['ref_count'] + counts['alt_count']
      counts['vaf'] = counts['alt_count'] / counts['depth']
      counts_with_depth = counts[counts['depth'] > 0]
      print(f"{total} variants total, {total - len(counts_with_depth)} variants removed due to 0 depth, {len(counts_with_depth)} final")
      counts_with_depth.to_csv(args.outcsv)

if __name__ == "__main__":
    main()
