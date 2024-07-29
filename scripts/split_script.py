# samtools sort -l 0 -t "CB" -O sam -o test.sam -@ 24 /<redacted>/SRCdata/SRC1rna_outs/possorted_genome_bam.bam

import os
import sys
import getopt
import fileinput
import math
import numpy as np
import pysam as ps
import pandas as pd
from collections import defaultdict

def print_to_file(current_tag, out_dir, samfile, read_list):
    name = current_tag.split('_')[0]
    outfile=out_dir+"/"+name+".bam"
    outbam = ps.AlignmentFile(outfile, "wb", template=samfile)
    for item in read_list:
      outbam.write(item)
    outbam.close()

argv = sys.argv[1:]
n_cells = "0"
sam_file_path = ""
barcode_file_path = ""
out_dir = ""
try:
  opts, args = getopt.getopt(argv,"b:c:o:")
except getopt.GetoptError:
  print("split_script4.py -b <path/to/bamfile.bam> -c <path/to/barcodes.tsv> -o <path/to/outfolder>")
  sys.exit(2)
for opt, arg in opts:
  if opt == "-c":
    barcode_file_path = arg
  elif opt == "-o":
    out_dir = arg
  elif opt == "-b":
    sam_file_path = arg
if sam_file_path == "":
  sys.exit("Missing bam_file, use -b option")
if barcode_file_path == "":
  sys.exit("Missing barcode_file, use -c option")
if out_dir == "":
  sys.exit("Missing out_dir, use -o option")

print(ps.__file__)
print(ps.__version__)
barcode_df = pd.read_csv(barcode_file_path, header=None, sep='\t')
barcodes = barcode_df.iloc[:,0].tolist()
samfile = ps.AlignmentFile(sam_file_path, "rb")
#print(samfile.lengths())
missed = 0
bad = 0
b = "finished " + str(0) + " million entries"
sys.stdout.write(b)
sys.stdout.flush()
current_tag = None
read_count = 0
cell_count = 0
read_list = []

for idx,read in enumerate(samfile.fetch(until_eof=True)):
    try:
      tag = read.get_tag("CB", with_value_type=False)
    except KeyError:
      bad = bad + 1
      continue
    if not (tag == current_tag):
      if (current_tag is not None):
        if (current_tag in barcodes):
          print_to_file(current_tag, out_dir, samfile, read_list)
          cell_count = cell_count + 1
        else:
          missed = missed + read_count
      read_list = []
      current_tag = tag
      read_count = 0
    read_list.append(read)
    read_count = read_count + 1
    if idx%10000000 == 0:
      sys.stdout.write("\b" * len(b))
      b = "finished " + str(idx//1000000) + " million entries"
      sys.stdout.write(b)
      sys.stdout.flush()
samfile.close()
print("\n", str(cell_count), " cells were found")
print(str(missed), " entries were from low cell count cells and were not included")
print(str(bad), " entries did not have a CB tag and were not included")
