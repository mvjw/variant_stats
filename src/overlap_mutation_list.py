import numpy as np
import pandas as pd
import matplotlib.pyplot as plt 
from typing import List

from src_alt.parsed_mutation_list import ParsedMutationList

class OverlapMutationList():
  array = None
  names = None
  overlap_matrix = None

  # TODO fix this, is not actually its own class, it just appends values to a
  # list of mutations, so it should be returned as a list of mutations and
  # then common variants should be its own class since it has a new structure
  def __init__(self, mut_filtered_array: List[ParsedMutationList], unfiltered=None):
    names = []
    mut_overlap_array = []
    overlap_matrix = np.zeros((len(mut_filtered_array), len(mut_filtered_array)))
    for filtered_mut in mut_filtered_array:
      assert isinstance(filtered_mut, ParsedMutationList)
      mut_overlap_array.append(filtered_mut.parsed_mutations)
    mut_comp_array = []
    if (unfiltered is not None) and (len(unfiltered) == len(mut_filtered_array)):
      for unfiltered_mut in unfiltered:
        mut_comp_array.append(unfiltered_mut.parsed_mutations)
    else:
      print("Unfiltered not available: using filtered array for comparisons")
      for filtered_mut in mut_filtered_array:
        mut_comp_array.append(filtered_mut.parsed_mutations)
    order = len(mut_overlap_array)
    overlap = []
    for i in range(order):
      overlap.append([])
      for j in [x for x in range(order) if x != i]:
        temp = (mut_overlap_array[j])[(mut_overlap_array[j])['map'].isin((mut_comp_array[i])['map'])]
        temp = (mut_comp_array[i])[(mut_comp_array[i])['map'].isin((mut_overlap_array[j])['map'])]
        print(mut_filtered_array[j].name + " overlap with " + mut_filtered_array[i].name + ": " + str(len(temp)))
        overlap_matrix[j,i] = len(temp)
        (overlap[i]).append(temp)
    for i in range(order):
      for j in range(order-1):
        (mut_overlap_array[i]) = pd.concat([(mut_overlap_array[i]), (overlap[i])[j]], ignore_index=True)
      (mut_overlap_array[i]).drop_duplicates(subset ='map', keep = 'first', ignore_index=True, inplace = True)
      (mut_overlap_array[i]).sort_values(by=['map'], ignore_index=True, inplace=True)
      print(mut_filtered_array[i].name + " final: " + str(len(mut_overlap_array[i])))
      overlap_matrix[i,i] = len(mut_overlap_array[i])
      names.append(mut_filtered_array[i].name)
    self.array = mut_overlap_array
    self.names = names
    self.overlap_matrix = overlap_matrix
  
  def get(self, name):
    idx = self.names.index(name)
    return self.array[idx]
  
  def get_overlap(self, name1, name2):
    # NOTE: in this function, order can matter!
    idx1 = self.names.index(name1)
    idx2 = self.names.index(name2)
    return self.overlap_matrix[idx1, idx2]
  
  # TODO add processing functions within the class (such as automatic overlap
  # calculation, filtering out variants that appear or dont appear in certain
  # datasets, etc.
  # and move up to be with other classes

### Essential Processing Functions ###

# common variants:
# calculates the number of datasets for which a particular variant is in common
# among datasets
# -> returns a single dataframe containing variants with their overlap with each
# dataset, as well as a common dataframe which counts the number of datasets
# that have a given variant in common
def common_variants(overlap_list: OverlapMutationList, cell_mut_list: List[ParsedMutationList]) -> pd.DataFrame:
    # TODO move inside overap class
    list_array = overlap_list.array
    order = len(list_array)
    common = pd.concat([list_array[0], list_array[1]], ignore_index=True)
    for i in range(2,order):
      common = pd.concat([common, list_array[i]], ignore_index=True)
    common.drop_duplicates(subset ='map', keep = 'first', ignore_index=True, inplace = True)
    for i in range(order):
      common[cell_mut_list[i].name] = common['map'].isin((list_array[i])['map'])
    common['common'] = 0 + common[cell_mut_list[0].name]
    for i in range(1, order):
      common['common'] = common['common'] + common[cell_mut_list[i].name]
    common.set_index('map', inplace=True)
    for i in range(order):
      common = common.join((list_array[i]).set_index('map')['quality'], rsuffix = '_' + cell_mut_list[i].name)
    for i in range(order):
      common = common.join((list_array[i]).set_index('map')['incidence'], rsuffix = '_' + cell_mut_list[i].name)
    return common.sort_values(by='gene').reset_index()

# redo common:
# recalculates the common column for a dataframe output of common_variants which
# has been modified or filtered to remove particular subsets of variants or datasets
# where the common column is the counts the number of datasets that have a given
# variant in common
def redo_common(common: pd.DataFrame, names: List[str]) -> pd.DataFrame:
    # TODO move inside overlap class
    order = len(names)
    common['common'] = 0 + common[names[0]]
    for i in range(1, order):
      common['common'] = common['common'] + common[names[i]]
    return common

### Helper Functions ###

# rank overlap:
# generate a barplot of the number of variants which are common to a given number of
# datasets
# -> common is an output from common_variants or redo_common
# -> outfile is the path to the output location of the boxplot
def rank_overlap(common: pd.DataFrame, outfile: str):
    # TODO move inside overlap class
    ranks = common['common'].value_counts(sort=False).sort_index()
    print(ranks)
    figure, axes = plt.subplots(1, 1, figsize=(16,16))
    ranks.plot(kind='bar', ax=axes)
    figure.savefig(outfile)