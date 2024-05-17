"""
This file contains the parsed VCF and cell objects which parse a single-cell vcf file (or any single vcf file)
and hold the parsed mutational information for that file, the results of which are used to generate a cell list
"""

import warnings

class ParsedVCF():
    number_variants = 0     # number of variants processed
    number_pass = 0         # number of variants which pass filters
    number_insertions = 0   # number of indel variants
    number_coding = 0       # number of variants in a protein-coding region (defined by snpEFF)
    pass_ = []              # list of entries for each variant: whether the variant passes filters
    indels_ = []            # list of entries for each variant: whether the variant is an indel
    multi_a_ = []           # list of entries for each variant: whether the variant is multi-allelic
    depths_ = []            # list of entries for each variant: the measured depth as defined by the total reads in the AD field
    vafs_ = []              # list of entries for each variant: the allele frequency as defined as "# supporting alleles / depth"
    impacts_ = []           # list of entries for each variant: the variant impact as defined by snpEFF
    labels_ = []            # list of entries for each variant: the variant label as defined by snpEFF
    allele_depths_ = []     # list of entries for each variant: the alt-allele count
    chromosomes_ = []       # list of entries for each variant: the chromosome of the variant
    positions_ = []         # list of entries for each variant: the position of the variant
    gene_ids_ = []          # list of entries for each variant: the variant gene as defined by snpEFF
    alleles_ = []           # list of entries for each variant: the unique allele id in format ref->alt
    qualities_ = []         # list of entries for each variant: the quality of the called variant
    bads_ = []              # list of entries for each variant: whether the variant is unusable for some reason
    cell_ids_ = []          # list of entries for each variant: the id of the cell for which the variant was called, used downstream
    n_multi = 0             # number of multiallelic variants
    head = ""               # header of the vcf file
    cell_barcode = ""       # barcode of the cell
    cell_id = 0             # id of the cell
    
    def __init__(self, filename: str):
        # TODO solve issues with cell_ids etc. that should only affect cells, not parsed vcfs
        variant_file = open(filename)
        lines = variant_file.readlines()
        variant_file.close()
        for line_o in lines:
            line = line_o.strip()
            if line[0] == '#':
                # line is header
                self.head = self.head + line_o
                if line[1] != '#':
                    # line is column head
                    line = line[1:]
                    line = line.split('\t')
                    chrom_index = line.index('CHROM')
                    pos_index = line.index('POS')
                    filter_index = line.index('FILTER')
                    info_index = line.index('INFO')
                    format_index = line.index('FORMAT')
                    ref_index = line.index('REF')
                    alt_index = line.index('ALT')
                    qual_index = line.index('QUAL')
                continue
            major_cols = line.split('\t')
            if len(major_cols) > format_index+2:
                raise Exception(
                    "{} contains more than one sample".format(filename))
            alts = major_cols[alt_index]
            alts = alts.split(',')
            for alt_n in range(len(alts)):
                if alts[alt_n] == '*':
                    # skip the allele as it is spanning and therefore documented elsewhere
                    continue
                # format contains quality and depth stats on the variant
                format_header = major_cols[format_index]
                format_header = format_header.split(':')
                format_body = major_cols[format_index+1]
                format_body = format_body.split(':')
                try:
                    # AD is depth for ref and alt alleles
                    AD_index = format_header.index('AD')
                    AD = format_body[AD_index]
                    AD = AD.split(',')
                    if len(AD) != len(alts)+1:  # !!! watch out, this may be incorrect
                        AD_i = float("NaN")
                    else:
                        AD_i = AD[alt_n+1]
                    self.allele_depths_.append(AD_i)
                except ValueError:
                    # for strelka somatic snps
                    try:
                        ref_dp_index = major_cols[ref_index] + 'U'
                        RD_index = format_header.index(ref_dp_index)
                        alt_dp_index = alts[alt_n] + 'U'
                        AD_index = format_header.index(alt_dp_index)
                        AD_alt = format_body[AD_index].split(',')[0]
                        AD_ref = format_body[RD_index].split(',')[0]
                        AD = [AD_ref, AD_alt]
                        self.allele_depths_.append(AD_alt)
                    except ValueError:
                        # for strelka somatic indels
                        try:
                            RD_index = format_header.index('TAR')
                            AD_index = format_header.index('TIR')
                            AD_alt = format_body[AD_index].split(',')[0]
                            AD_ref = format_body[RD_index].split(',')[0]
                            AD = [AD_ref, AD_alt]
                            self.allele_depths_.append(AD_alt)
                        except ValueError:
                            # for samtools (?) or varscan (?)
                            info_body = major_cols[info_index]
                            info_body_split = info_body.split(';')
                            info_body_dict = {i.split('=')[0]: i.split(
                                '=')[-1] for i in info_body_split}
                            DP4 = info_body_dict['DP4'].split(',')
                            AD_alt = int(DP4[2]) + int(DP4[3])
                            AD_ref = int(DP4[0]) + int(DP4[1])
                            AD = [AD_ref, AD_alt]
                            self.allele_depths_.append(AD_alt)
                if len(alts) != 2:
                    self.n_multi = self.n_multi + 1
                    self.multi_a_.append(True)
                else:
                    self.multi_a_.append(False)
                try:
                    self.qualities_.append(float(major_cols[qual_index]))
                except ValueError:
                    if major_cols[qual_index] == ".":
                        self.qualities_.append(float(0.0))
                    else:
                        raise Exception("bad!")
                # self.lines_.append(line_o)
                self.number_variants = self.number_variants + 1
                chromosome = major_cols[chrom_index]
                if 'chr' in chromosome:
                    chromosome = chromosome.split('chr')[1]
                self.chromosomes_.append(chromosome)
                self.positions_.append(int(major_cols[pos_index]))
                allele = major_cols[ref_index] + '->' + alts[alt_n]
                self.alleles_.append(allele)
                if major_cols[filter_index] == 'PASS':
                    self.number_pass = self.number_pass + 1
                    self.pass_.append(True)
                else:
                    self.pass_.append(False)
                # The following uses DP to calculate depth, except that it was realized that AD should work fine for _most_ depth calculations (remember, indels filter out low-quality reads from this count)
                # DP_index = None # DP is total depth for variant
                # try:
                #    DP_index = format_header.index('DP')
                # except ValueError:
                #    DPI_index = format_header.index('DPI')
                # if DP_index is not None:
                #    DP = format_body[DP_index]
                # else:
                #    self.number_insertions = self.number_insertions + 1
                #    DP = format_body[DPI_index]
                # self.depths_.append(int(DP))
                DPI_index = None
                try:
                    DPI_index = format_header.index('DPI')
                except ValueError:
                    pass
                if DPI_index is not None:  # not sure if this is best way to check for indels
                    self.number_insertions = self.number_insertions + 1
                    self.indels_.append(True)
                    DPI = format_body[DPI_index]
                else:
                    self.indels_.append(False)
                depth = sum([int(i) for i in AD])
                if depth == 0:
                    try:
                        if int(DPI) != 0:
                            depth = int(DPI)
                            self.bads_.append(False)
                        else:
                            self.bads_.append(True)
                    except:
                        self.bads_.append(True)
                else:
                    self.bads_.append(False)
                self.depths_.append(depth)
                try:
                    # this is incorrect, figure out how to deal with multiple alt alleles
                    VAF = int(AD[1]) / depth
                except:
                    VAF = -1
                self.vafs_.append(VAF)
                info_body = major_cols[info_index]
                info_body = info_body.split(';')
                try:
                    ann = [i for i in info_body if 'ANN=' in i][0]
                except:
                    self.impacts_.append('None')
                    self.gene_ids_.append('None')
                    self.labels_.append('None')
                    warnings.warn(
                        "Could not find annotation for line:\n{}".format(line_o))
                    continue
                ann = ann.split('=')[1]
                ann = ann.split(',')
                # the following work fine with the current version of snpEFF,
                # but may not work in the future
                allele_index = 0  # look at INFO=<ID=ANN...
                coding_index = 7  # look at INFO=<ID=ANN...
                label_index = 1  # look at INFO=<ID=ANN...
                impact_index = 2  # look at INFO=<ID=ANN...
                gene_index = 3  # look at INFO=<ID=ANN...
                curr_imp_len = len(self.impacts_)
                curr_lab_len = len(self.labels_)
                curr_gid_len = len(self.gene_ids_)
                for ann_i in range(len(ann)):
                    ann_ = ann[ann_i]
                    ann_ = ann_.split('|')
                    if ann_[allele_index] == alts[alt_n]:
                        if ann_[coding_index] == 'protein_coding':
                            self.impacts_.append(ann_[impact_index])
                            self.gene_ids_.append(ann_[gene_index])
                            self.number_coding = self.number_coding + 1
                        else:
                            self.impacts_.append('non_coding')
                            self.gene_ids_.append('NA')
                        self.labels_.append(ann_[label_index])
                        break
                # temporary code to double check that offset of appends is not broken in for loop
                if not ((len(self.impacts_) == curr_imp_len+1) and (len(self.labels_) == curr_lab_len+1) and (len(self.gene_ids_) == curr_gid_len+1)):
                    print(f"{len(self.impacts_)} {curr_imp_len+1}")
                    print(f"{len(self.labels_)} {curr_lab_len+1}")
                    print(f"{len(self.gene_ids_)} {curr_gid_len+1}")
                    print(f"{line_o}")
                    print(f"{alts[alt_n]} {ann_[allele_index]}")
                if not (len(self.allele_depths_) == len(self.multi_a_) == len(self.qualities_) == len(self.chromosomes_) == len(self.positions_) == len(self.alleles_) == len(self.pass_) == len(self.indels_) == len(self.bads_) == len(self.depths_) == len(self.vafs_) == len(self.impacts_) == len(self.labels_) == len(self.impacts_) == len(self.labels_) == len(self.gene_ids_) == self.number_variants):
                    # raise Exception("Error in parse_vcf: list length mismatch during parsing")
                    warnings.warn(
                        "Warning in parse_vcf: list length mismatch during parsing")
                    self.impacts_.append(None)
                    self.labels_.append(None)
                    self.gene_ids_.append(None)
        self.cell_ids_ = [0 for i in range(self.number_variants)]


class Cell(ParsedVCF):
  """ This class is identical to the ParsedVCF except that it contains a unique cell_id and barcode """
  def __init__(self, filename: str, cell_id: int):
    # TODO add asserts and catches here
    super().__init__(filename)
    self.cell_barcode = (filename.split("/annotated.vcf")[0]).split("/")[-1]
    self.cell_id = cell_id
    self.cell_ids_ = [cell_id for i in range(self.number_variants)]
