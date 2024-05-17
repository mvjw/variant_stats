import os
import pandas as pd
import warnings
from typing import List, Union

from src_alt.cell_list import CellList


class ParsedMutationList():
    name = None
    header = None
    n_cells = None
    parsed_mutations = None
    
    def __init__(self, in_var: CellList):
        if isinstance(in_var, CellList):
            self.from_cell_list(in_var)
        elif isinstance(in_var, str):
            self.load_from_folder(in_var)

    def from_cell_list(self, cell_list: CellList):
        df = cell_list.mutations.copy()
        df['calls'] = (df.loc[:, 'vaf']*df.loc[:, 'depth']).copy()
        # df['pass'] = df['pass'].astype(int)
        self.name = cell_list.name
        self.header = cell_list.header
        self.n_cells = cell_list.n_cells
        group = df.groupby(['map'])
        self.parsed_mutations = group.count().reset_index(
        ).sort_values(by=['map'], ascending=True)
        self.parsed_mutations.rename(
            columns={list(self.parsed_mutations)[1]: 'count'}, inplace=True)
        df_chromosome = group['chromosome'].first().reset_index().sort_values(by=[
            'map'], ascending=True)
        df_position = group['position'].first().reset_index().sort_values(by=[
            'map'], ascending=True)
        df_pass = group['pass'].mean().reset_index(
        ).sort_values(by=['map'], ascending=True)
        df_depth = group['depth'].mean().reset_index(
        ).sort_values(by=['map'], ascending=True)
        df_coverage = group['depth'].sum().reset_index(
        ).sort_values(by=['map'], ascending=True)
        df_vaf_avg = group['vaf'].mean().reset_index(
        ).sort_values(by=['map'], ascending=True)
        df_calls = group['calls'].sum().reset_index(
        ).sort_values(by=['map'], ascending=True)
        df_gene = group['gene'].first().reset_index(
        ).sort_values(by=['map'], ascending=True)
        df_allele = group['allele'].first().reset_index(
        ).sort_values(by=['map'], ascending=True)
        df_label = group['label'].first().reset_index(
        ).sort_values(by=['map'], ascending=True)
        df_impact = group['impact'].first().reset_index(
        ).sort_values(by=['map'], ascending=True)
        # df_quality = group['quality'].mean().reset_index().sort_values(by=['map'], ascending=True)
        # df_quality = group['quality'].median().reset_index().sort_values(by=['map'], ascending=True)
        df_quality = group['quality'].max().reset_index(
        ).sort_values(by=['map'], ascending=True)
        self.parsed_mutations['chromosome'] = df_chromosome[list(df_chromosome)[
            1]]
        self.parsed_mutations['position'] = df_position[list(df_position)[1]]
        self.parsed_mutations['gene'] = df_gene[list(df_gene)[1]]
        self.parsed_mutations['allele'] = df_allele[list(df_allele)[1]]
        self.parsed_mutations['label'] = df_label[list(df_label)[1]]
        self.parsed_mutations['impact'] = df_impact[list(df_impact)[1]]
        self.parsed_mutations['pass'] = df_pass[list(df_pass)[1]]
        self.parsed_mutations['depth'] = df_depth[list(df_depth)[1]]
        self.parsed_mutations['coverage'] = df_coverage[list(df_coverage)[1]]
        self.parsed_mutations['calls'] = df_calls[list(df_calls)[1]]
        self.parsed_mutations['avg vaf'] = df_vaf_avg[list(df_vaf_avg)[1]]
        self.parsed_mutations['weighted avg vaf'] = df_calls[list(
            df_calls)[1]]/self.parsed_mutations['coverage']
        self.parsed_mutations['quality'] = df_quality[list(df_quality)[1]]
        self.parsed_mutations['incidence'] = self.parsed_mutations['count']/self.n_cells
        self.parsed_mutations.sort_values(
            by=['chromosome', 'position'], ascending=True, inplace=True)
        self.parsed_mutations.reset_index(drop=True, inplace=True)

    def load_from_folder(self, path: str):
        if not os.path.exists(path):
            raise Exception(path + " does not exist")
        f = open(path + "/header.txt", "r")
        lines = f.readlines()
        f.close()
        self.name = path.split("/")[-1]
        self.n_cells = int(lines[0].split(":")[1])
        self.header = lines[1:]
        self.parsed_mutations = pd.read_csv(path + "/body.tsv", sep='\t')
        print("loaded " + path)

    def save_to_folder(self, path: str):
        path_to_folder = path + "/" + self.name
        print(path_to_folder)
        if not os.path.exists(path_to_folder):
            try:
                os.mkdir(path_to_folder)
            except OSError:
                print("Could not create directory from save_to_folder()")
                return
        header_file_path = path_to_folder + "/header.txt"
        f = open(header_file_path, "w")
        f.write("n_cells:" + str(self.n_cells) + "\n")
        f.write(self.header)
        f.close
        body_file_path = path_to_folder + "/body.tsv"
        self.parsed_mutations.to_csv(body_file_path, sep='\t', index=False)
        print("saved " + path_to_folder)



### helper functions ###

def save_mutation_array(cell_mut_array: Union[ParsedMutationList, List[ParsedMutationList]], path_to_folder: str):
    if not isinstance(cell_mut_array, list):
      cell_mut_array = [cell_mut_array]
    if not os.path.exists(path_to_folder):
      try:
        os.mkdir(path_to_folder)
      except OSError:
        print("Could not create directory from save_mutations()")
        return
    for i, mut_array in enumerate(cell_mut_array):
      assert isinstance(mut_array, ParsedMutationList)
      mut_array.save_to_folder(path_to_folder)
    
def load_mutation_array(path_to_folder: str, names: str) -> ParsedMutationList:
    mutation_array = []
    for name in names:
      pathname = path_to_folder + "/" + name
      if os.path.isdir(pathname):
        mutation_array.append(ParsedMutationList(pathname))
      else:
        raise Exception("Could not find file " + pathname + "\n")
    return mutation_array

def get_mutations(list_cell_lists: List[CellList]) -> List[ParsedMutationList]:
    if not isinstance(list_cell_lists, list):
      list_cell_lists = [list_cell_lists]
    cell_parsed_array = []
    for i, list_cell in enumerate(list_cell_lists):
      assert isinstance(list_cell, CellList)
      cell_parsed_raw = ParsedMutationList(list_cell)
      cell_parsed_array.append(cell_parsed_raw)
      print("finished parsing " + cell_parsed_raw.name)
    return cell_parsed_array

# generate vcf:
# Turn a mutation list dataframe into a vcf file
# The mutation list should have unique mutations and 'chromosome', 'position', and 'allele' fields
# as generated by a parsed mutation list 'parsed_mutations' dataframe
# Likewise, the header should be taken from the 'header' contained in a parsed mutation list
def generate_vcf(mutation_list: pd.DataFrame, header: str, file_name: str, atac=False):
    # TODO update to use mutation list header
    #header = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
    if len(mutation_list) == 0:
      warnings.warn("Skipped {} as dataframe was empty".format(file_name))
      return
    lines = mutation_list.sort_values(by=['chromosome','position'])['chromosome'].astype(str)
    if(atac):
      lines = "chr" + lines
    lines = lines + "\t" + mutation_list.sort_values(by=['chromosome','position'])['position'].astype(str)
    lines = lines + "\t."
    lines = lines + "\t" + mutation_list.sort_values(by=['chromosome','position'])['allele'].str.split('-', n = 1, expand = True)[0]
    lines = lines + "\t" + mutation_list.sort_values(by=['chromosome','position'])['allele'].str.split('>', n = 1, expand = True)[1]
    lines = lines + "\t0"
    lines = lines + "\tPASS"
    #lines = lines + "\tNA=0\n"
    lines = lines + "\t\n"
    outfile = open(file_name, "w")
    try:
      outfile.writelines(header.header)
    except:
      f = open(header, "r")
      header = f.readlines()
      f.close()
      outfile.writelines(header)
    for line in lines:
      outfile.write(line)
    outfile.close()