# filter known PTM sites
#
def filter_known_ptm_sites(mutation_list, ptm_list, isbed=True):
    # TODO find some way to do this in terms of the overlap class
    variants = mutation_list.copy()
    ptms = ptm_list.copy()
    print(ptms.head(n=5))
    if(('chr' in str(ptms['chromosome'].iloc[0])) & (not 'chr' in str(variants['chromosome'].iloc[0]))):
      variants['pmap'] = 'chr' + variants['chromosome'].astype(str) + ':' + variants['position'].astype(int).astype(str)
    else:
      variants['pmap'] = variants['chromosome'].astype(str) + ':' + variants['position'].astype(int).astype(str)
    if(isbed):
      ptms['pmap'] = ptms['chromosome'].astype(str) + ':' + (ptms['position']+1).astype(int).astype(str)
    else:
      ptms['pmap'] = ptms['chromosome'].astype(str) + ':' + ptms['position'].astype(int).astype(str)
    print(variants['pmap'].head(n=10))
    print(ptms['pmap'].head(n=10))
    len1 = variants.shape[0]
    removed = variants[variants['pmap'].isin(ptms['pmap'])]
    variants = variants[~(variants['pmap'].isin(ptms['pmap']))]
    print(removed)
    len2 = variants.shape[0]
    print('removed ', str(len1-len2), ' known PTM sites')
    return variants