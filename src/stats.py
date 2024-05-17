import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

from src_alt.plot import axplot, violin, scatter, bar, hist, axtitle
from src_alt.cell_list import CellList
from src_alt.parsed_mutation_list import ParsedMutationList

def cell_stats(cell_list_: CellList, out_name: str):
    assert isinstance(cell_list_, CellList)
    my_list_mutations = cell_list_.mutations.copy()
    my_list_cells = cell_list_.cells.copy()
    f = plt.figure(figsize=(16,16))
    gs = gridspec.GridSpec(nrows=4, ncols=3, height_ratios=[1, 1, 1, 1])
    axes = []
    axplot(f, axes, my_list_cells['mut_per_cell'], violin, gs, (0,0), 'Mutations/Cell')
    axplot(f, axes, my_list_cells['pass_per_cell'], violin, gs, (0,1), 'Pass Mutations/Cell')
    axplot(f, axes, (my_list_cells['mut_per_cell'], my_list_cells['pass_per_cell']), scatter, gs, (0,2), 'Pass vs Total Mutations')
    axplot(f, axes, my_list_cells['coding_per_cell'], violin, gs, (1,0), 'Coding Mutations/Cell')
    coding_mutations = my_list_mutations[my_list_mutations['impact']!='non_coding'].copy()
    impacts, counts = np.unique(coding_mutations['impact'], return_counts=True)
    axplot(f, axes, (impacts, counts), bar, gs, (1,1), 'Impact Counts')
    my_list_mutations['is_impact'] = (my_list_mutations['impact']!='MODIFIER') & (my_list_mutations['impact']!='non_coding')
    my_list_cells['impact_per_cell'] = my_list_mutations.groupby(['id'])['is_impact'].sum().reset_index().sort_values(by=['id'])['is_impact']
    axplot(f, axes, my_list_cells['impact_per_cell'], violin, gs, (1,2), 'Non-Modifier Coding Mutations/Cell')
    axplot(f, axes, my_list_cells['avg_depth_per_cell'], violin, gs, (2,0), 'Avg Mutation Coverage/Cell')
    axplot(f, axes, (my_list_mutations['depth'], 10, (1,10)), hist, gs, (2,1), 'Total Coverage/Mutation')
    vafs, counts = np.unique(my_list_mutations['vaf'], return_counts=True)
    axplot(f, axes, (vafs, counts), scatter, gs, (2,2), 'VAF Counts/Mutation')
    impact_mutations = my_list_mutations[my_list_mutations['is_impact']].copy()
    my_list_cells['impact_avg_depth_per_cell'] = impact_mutations.groupby(['id'])['depth'].mean().reset_index().sort_values(by=['id'])['depth']
    axplot(f, axes, my_list_cells['impact_avg_depth_per_cell'], violin, gs, (3,0), 'Avg Impactful Mutation Coverage/Cell')
    axplot(f, axes, (impact_mutations['depth'], 10, (1,10)), hist, gs, (3,1), 'Total Coverage/Impactful Mutation')
    vafs, counts = np.unique(impact_mutations['vaf'], return_counts=True)
    axplot(f, axes, (vafs, counts), scatter, gs, (3,2), 'VAF Counts/Impactful Mutation')
    plt.show()
    f.savefig(out_name)

def variant_stats(mutations: ParsedMutationList, out_name: str):
    # NOTE: do not use a mutation list that has been mixed with another (such as after finding overlap between two mutation lists) as this will probably provide inconsistent statistics
    variants_list = mutations.copy()
    f = plt.figure(figsize=(16,16))
    gs = gridspec.GridSpec(nrows=3, ncols=2, height_ratios=[1, 1, 1])
    axes = []
    axplot(f, axes, (variants_list['count'], 100, (1,100)), hist, gs, (0,0), 'Variant Cell Number')
    axtitle(axes, '# cells with variant', 'log histogram count')
    axplot(f, axes, (variants_list['depth'], 100, (1,1.1)), hist, gs, (0,1), 'Variant Average Depths')
    axtitle(axes, 'average depth per cell of variant (when called)', 'log histogram count')
    axplot(f, axes, (variants_list['count'], variants_list['depth']), scatter, gs, (1,0), 'Average Depth vs Number Cells')
    axes[-1].set_xlim([0,100])
    axes[-1].set_ylim([0,40])
    axtitle(axes, '# cells with variant', 'average depth per cell')
    axplot(f, axes, (variants_list['coverage'], 100, (1,100)), hist, gs, (1,1), 'Total Variant Site Coverage')
    axtitle(axes, '# reads at variant site (ignoring cells in which the variant was not called)', 'log histogram count')
    axplot(f, axes, (variants_list['count'],variants_list['coverage']), scatter, gs, (2,0), 'Total Coverage vs Number Cells') 
    axes[-1].set_xlim([0,100])
    axes[-1].set_ylim([0,100])
    axtitle(axes, '# cells with variant', 'reads at variant site')
    axplot(f, axes, (variants_list['vaf'], 100, (0.0, 1.0)), hist, gs, (2,1), 'Variant overall VAF across Cells')
    axtitle(axes, 'vaf for variant', 'log histogram count')
    f.savefig(out_name)