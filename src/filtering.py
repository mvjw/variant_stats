import copy
from typing import List

from src_alt.parsed_mutation_list import ParsedMutationList

def filter_mutations(cell_mut_array: List[ParsedMutationList], 
                     cell_filter_: int, 
                     alt_read_filter_: int, 
                     vaf_filter_: float, 
                     depth_filter_: int, 
                     qual_filter_: float, 
                     exonic=False, 
                     splice=False, 
                     max_cell_filter_: int = None, 
                     max_alt_read_filter_: int = None, 
                     limexonic=False
                     ) -> List[ParsedMutationList]:
    # TODO move to within parsed mutations list class
    if not isinstance(cell_mut_array, list):
      cell_mut_array = [cell_mut_array]
    cell_mut_array_copy = []
    for cell_mut in cell_mut_array:
      assert isinstance(cell_mut, ParsedMutationList)
      cell_mut_array_copy.append(copy.deepcopy(cell_mut))
    cell_parsed_array = []
    for cell_mut in cell_mut_array_copy:
      cell_parsed_array.append(cell_mut.parsed_mutations)
    order = len(cell_parsed_array)
    if(exonic):
      for i in range(order):
        #cell_parsed_array[i] = (cell_parsed_array[i])[((cell_parsed_array[i])['impact'] != 'non_coding') & ((cell_parsed_array[i])['impact'] != 'MODIFIER')]
        cell_parsed_array[i] = (cell_parsed_array[i])[((cell_parsed_array[i])['impact'] != 'non_coding') & ((cell_parsed_array[i])['impact'] != 'MODIFIER')]
        print(cell_mut_array_copy[i].name + " total: " + str(len(cell_parsed_array[i])))
    elif(limexonic):
      for i in range(order):
        cell_parsed_array[i] = (cell_parsed_array[i])[((cell_parsed_array[i])['impact'] != 'non_coding') & ((cell_parsed_array[i])['label'] != 'upstream_gene_variant') & ((cell_parsed_array[i])['label'] != 'downstream_gene_variant')]
        print(cell_mut_array_copy[i].name + " total: " + str(len(cell_parsed_array[i])))
    cell_filtered_array = []
    for i in range(order):
      tmp_fil_1 = (cell_parsed_array[i])[(cell_parsed_array[i])['count'].astype(int) > cell_filter_]
      tmp_fil_2 = tmp_fil_1[tmp_fil_1['weighted avg vaf'].astype(float) > vaf_filter_]
      tmp_fil_3 = tmp_fil_2[tmp_fil_2['coverage'].astype(int) > depth_filter_]
      tmp_fil_4 = tmp_fil_3[tmp_fil_3['quality'].astype(float) > qual_filter_]
      tmp_fil_5 = tmp_fil_4[tmp_fil_4['calls'].astype(int) > alt_read_filter_]
      if(splice):
        tmp_fil_6 = tmp_fil_5
      else:
        tmp_fil_6 = tmp_fil_5[~(tmp_fil_5['label'].str.contains("splice_acceptor_variant", regex=False) | tmp_fil_5['label'].str.contains("splice_donor_variant", regex=False) | tmp_fil_5['label'].str.contains("splice_region_variant", regex=False) | tmp_fil_5['label'].str.contains("intron_variant", regex=False))]
      if(max_cell_filter_):
        tmp_fil_7 = tmp_fil_6[tmp_fil_6['count'] <= max_cell_filter_]
      else:
        tmp_fil_7 = tmp_fil_6
      if(max_alt_read_filter_):
        tmp_fil_8 = tmp_fil_7[tmp_fil_7['calls'] <= max_alt_read_filter_]
      else:
        tmp_fil_8 = tmp_fil_7
      cell_filtered_array.append(tmp_fil_8)
      print(cell_mut_array_copy[i].name + " filtered: " + str(len(cell_filtered_array[i])))
    for i, cell_filtered in enumerate(cell_filtered_array):
      cell_mut_array_copy[i].parsed_mutations = cell_filtered.reset_index(drop=True)
    return cell_mut_array_copy    