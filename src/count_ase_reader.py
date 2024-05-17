import pandas as pd

def count_ase_reader(table_path):
    counts = pd.read_csv(table_path, sep="\t")
    counts['map'] = counts['contig'].astype(str) + ':' + counts['position'].astype(str) + ':' + counts['refAllele'].astype(str) + '->' + counts['altAllele'].astype(str)
    counts['allele'] = counts['refAllele'].astype(str) + '->' + counts['altAllele'].astype(str)
    counts.rename(columns = {'contig':'chromosome', 'refCount':'ref_count', 'altCount':'alt_count'}, inplace = True)
    return counts