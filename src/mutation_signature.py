import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# mutation signatures:
# generate a barplot of the number of variants which follow a specific mutation
# signature such as C->T
# -> mutation list is any dataframe containing variants with relevant fields of gene, chromosome, position, impact, label, and allele
# -> outfile is the path to the output location of the boxplot
def mutation_signatures(mutation_list: pd.DataFrame, outfile: str):
    # TODO update to use parsed mutation list class
    variants = mutation_list[['gene', 'impact', 'label', 'chromosome', 'position', 'allele']]
    variants['ref'] = variants['allele'].str.split('-').str[0]
    variants['alt'] = variants['allele'].str.split('>').str[1]
    signature_indices = ["A", "C", "T", "G"]
    titles = []
    list = []
    for i, ref in enumerate(signature_indices):
      for j, alt in enumerate(signature_indices):
        if(i==j):
          continue        
        matches = (variants['ref']==ref) & (variants['alt']==alt)
        list.append(matches.sum())
        titles.append(ref + "->" + alt)
    list.append(len(variants)-sum(list))
    titles.append("N/A")
    figure, axes = plt.subplots(1, 1, figsize=(16,16))
    axes.bar(np.arange(len(list)), list)
    axes.set_xticks(np.arange(len(list)))
    axes.set_xticklabels(titles)
    figure.savefig(outfile)
    
# mutation signatures stranded:
# generate a barplot of the number of variants which follow a specific mutation
# signature such as C->T, but corrected for strand direction
# -> mutation list is any dataframe containing variants with relevant fields of gene, chromosome, position, impact, label, and allele
# -> outfile is the path to the output location of the boxplot
def mutation_signatures_stranded(mutation_list: pd.DataFrame, outfile: str):
    # TODO update to use parsed mutation list class
    variants = mutation_list[['gene', 'impact', 'label', 'chromosome', 'position', 'allele', 'strand']]
    variants['ref'] = variants['allele'].str.split('-').str[0]
    variants['alt'] = variants['allele'].str.split('>').str[1]
    signature_indices = ["A", "C", "T", "G"]
    signature_compliments = {"A":"T", "G":"C", "C":"G", "T":"A"}
    classes = ["+", "-", ".", "?"]
    titles = []
    lists = {}
    for h, cla in enumerate(classes):
      list = []
      for i, ref in enumerate(signature_indices):
        for j, alt in enumerate(signature_indices):
          if(i==j):
            continue  
          if(cla == "-"):
            b1 = signature_compliments[ref]
            b2 = signature_compliments[alt]
          else:
            b1 = ref
            b2 = alt
          matches = (variants['ref']==b1) & (variants['alt']==b2) & (variants['strand']==cla)
          list.append(matches.sum())
          if cla == "+":
            titles.append(ref + "->" + alt)
      list.append(sum(variants['strand']==cla)-sum(list))
      lists[cla] = list
      if cla == "+":
        titles.append("N/A")
    df = pd.DataFrame(lists, index=titles)
    ax = df.plot(kind='bar', stacked=True, figsize=(10, 6))
    ax.set_ylabel('No. Mutations')
    plt.legend(title='labels', bbox_to_anchor=(1.0, 1), loc='upper left')
    plt.savefig(outfile)