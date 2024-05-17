import os
import numpy as np
import pandas as pd
import multiprocessing

from typing import List
from src_alt.parsed_vcf import Cell

class CellList():
    name: str = None
    n_cells: int = 0
    header: str = None
    cells: pd.DataFrame = None
    mutations: pd.DataFrame = None

    def __init__(self, list_of_cells: List[Cell], name: str):
        # TODO add asserts here
        self.name = name
        self.n_cells = len(list_of_cells)
        self.header = list_of_cells[0].head
        cell_ids = [c.cell_id for c in list_of_cells]
        barcodes = [c.cell_barcode for c in list_of_cells]
        mut_per_cell = [c.number_variants for c in list_of_cells]
        pass_per_cell = [c.number_pass for c in list_of_cells]
        insert_per_cell = [c.number_insertions for c in list_of_cells]
        coding_per_cell = [c.number_coding for c in list_of_cells]
        # impact_per_cell=[sum((i != 'MODIFIER') & (i != 'non_coding') for i in c.impacts_) for c in list_of_cells]
        avg_depth_per_cell = [((sum(c.depths_)/c.number_variants)
                               if c.number_variants > 0 else 0) for c in list_of_cells]

        self.cells = pd.DataFrame.from_dict({'id': cell_ids})
        self.cells['barcode'] = barcodes
        self.cells['mut_per_cell'] = mut_per_cell
        self.cells['pass_per_cell'] = pass_per_cell
        self.cells['insert_per_cell'] = insert_per_cell
        self.cells['coding_per_cell'] = coding_per_cell
        # self.cells['impact_per_cell']=impact_per_cell
        self.cells['avg_depth_per_cell'] = avg_depth_per_cell
        self.cells.sort_values(by=['id'])

        id_tag = [c.cell_ids_ for c in list_of_cells]
        chr_tag = [c.chromosomes_ for c in list_of_cells]
        pos_tag = [c.positions_ for c in list_of_cells]
        gene_tag = [c.gene_ids_ for c in list_of_cells]
        allele_tag = [c.alleles_ for c in list_of_cells]
        pass_tag = [c.pass_ for c in list_of_cells]
        depth_tag = [c.depths_ for c in list_of_cells]
        vaf_tag = [c.vafs_ for c in list_of_cells]
        label_tag = [c.labels_ for c in list_of_cells]
        impact_tag = [c.impacts_ for c in list_of_cells]
        multia_tag = [c.multi_a_ for c in list_of_cells]
        indel_tag = [c.indels_ for c in list_of_cells]
        qual_tag = [c.qualities_ for c in list_of_cells]

        t_id_tag = np.hstack(np.array([np.asarray(i)
                             for i in id_tag], dtype=object))
        t_chr_tag = np.hstack(
            np.array([np.asarray(i) for i in chr_tag], dtype=object))
        t_pos_tag = np.hstack(
            np.array([np.asarray(i) for i in pos_tag], dtype=object))
        t_gene_tag = np.hstack(
            np.array([np.asarray(i) for i in gene_tag], dtype=object))
        t_allele_tag = np.hstack(
            np.array([np.asarray(i) for i in allele_tag], dtype=object))
        t_pass_tag = np.hstack(
            np.array([np.asarray(i) for i in pass_tag], dtype=object))
        t_depth_tag = np.hstack(
            np.array([np.asarray(i) for i in depth_tag], dtype=object))
        t_vaf_tag = np.hstack(
            np.array([np.asarray(i) for i in vaf_tag], dtype=object))
        t_label_tag = np.hstack(
            np.array([np.asarray(i) for i in label_tag], dtype=object))
        t_impact_tag = np.hstack(
            np.array([np.asarray(i) for i in impact_tag], dtype=object))
        t_multia_tag = np.hstack(
            np.array([np.asarray(i) for i in multia_tag], dtype=object))
        t_indel_tag = np.hstack(
            np.array([np.asarray(i) for i in indel_tag], dtype=object))
        t_qual_tag = np.hstack(
            np.array([np.asarray(i) for i in qual_tag], dtype=object))

        self.mutations = pd.DataFrame.from_dict({'id': t_id_tag})
        self.mutations['chromosome'] = t_chr_tag
        self.mutations['position'] = t_pos_tag
        self.mutations['gene'] = t_gene_tag
        self.mutations['allele'] = t_allele_tag
        self.mutations['map'] = self.mutations['chromosome'].astype(
            str) + ':' + self.mutations['position'].astype(int).astype(str) + ':' + self.mutations['allele']
        self.mutations['pass'] = t_pass_tag.astype(int)
        self.mutations['depth'] = t_depth_tag.astype(int)
        self.mutations['vaf'] = t_vaf_tag.astype(float)
        self.mutations['label'] = t_label_tag
        self.mutations['impact'] = t_impact_tag
        self.mutations['multi_a'] = t_multia_tag
        self.mutations['indel'] = t_indel_tag
        self.mutations['quality'] = t_qual_tag.astype(float)
        self.mutations.sort_values(by=['id'])

        # not yet tracked
        # self.allele_depths_ = []
        # self.lines_ = []
        # self.bads_ = []

    def load_from_folder(self, path):
        if not os.path.exists(path):
            raise Exception(path + " does not exist")
        f = open(path + "/header.txt", "r")
        lines = f.readlines()
        f.close()
        self.name = path.split("/")[-1]
        self.n_cells = int(lines[0].split(":")[1])
        self.header = lines[1:]
        self.parsed_mutations = pd.read_csv(
            path + "/body.tsv", sep='\t', dtype=str)
        print("loaded " + path)

    def __init__(self, in_var):
        if isinstance(in_var, CellList):
            self.from_cell_list(in_var)
        elif isinstance(in_var, str):
            self.load_from_folder(in_var)

    def save_to_folder(self, path):
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

### Bulk Load functions ###

def get_cells(path: str, name: str, nthreads: int = 1) -> CellList:
    cells = []
    directory_list = open(path)
    directories = directory_list.read().splitlines()
    directory_list.close()
    #directories = directories[:100] # only first 100 lines for testing
    if nthreads >= 1:
        file_paths = [d + 'annotated.vcf' for d in directories]
        path_ind_tuples = zip(file_paths, range(len(file_paths)))
        with multiprocessing.Pool(processes=nthreads) as pool:
            cells = pool.starmap(Cell, path_ind_tuples)
        my_list = CellList(cells, name)
    else:
        for i, direct in enumerate(directories):
            if len(direct) > 1:
                direct = direct + 'annotated.vcf'
                cells.append(Cell(direct, i))
            my_list = CellList(cells, name)
    return my_list

def get_cell_lists(path_list: List[str], name_list: List[str], nthreads: int = 1) -> List[CellList]:
    if not isinstance(path_list, list):
      path_list = [path_list]
    if not isinstance(name_list, list):
      name_list = [name_list]
    assert len(path_list) == len(name_list)
    list_cell_list = []
    for i, path in enumerate(path_list):
      cell_list = get_cells(path, name_list[i], nthreads=nthreads)
      list_cell_list.append(cell_list)
      print("finished loading " + cell_list.name)
    return list_cell_list