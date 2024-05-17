import pandas as pd

def to_bed(variants: pd.DataFrame, outfile) -> None:
    variant_indices = variants[['gene', 'impact', 'label', 'count', 'new_count', 'quality', 'new_qual', 'chromosome', 'position']]
    variant_indices['tail'] = variant_indices['position'] + 1
    variant_indices['ref'] = variants['allele'].str.split('-', n = 1, expand = True)[0]
    variant_indices['alt'] = variants['allele'].str.split('>', n = 1, expand = True)[1]
    variant_indices.to_csv(outfile, sep='\t', index=False)

