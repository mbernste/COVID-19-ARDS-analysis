import matplotlib as mpl
#mpl.use('Agg')
import seaborn as sns
from matplotlib import pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd
import scanpy as sc
import sys
import os 
from os.path import join
from optparse import OptionParser
from collections import defaultdict
import gseapy as gp
import json


GENE_SETS = ['GO_Biological_Process_2018', 'GO_Molecular_Function_2018', 'KEGG_2016']
GSEA_THRESH = 0.05

def main():
    usage = "" # TODO
    parser = OptionParser(usage=usage)
    parser.add_option("-o", "--out_file", help="OUTPUT file")
    (options, args) = parser.parse_args()

    de_up_genes_f = args[0]
    de_down_genes_f = args[1]
    out_f = options.out_file

    with open(de_up_genes_f, 'r') as f:
        up_genes = [l for l in f]
    with open(de_down_genes_f, 'r') as f:
        down_genes = [l for l in f]
    de_genes = set(up_genes) | set(down_genes)

    db_to_gene_sets = {}
    for db in GENE_SETS:
        enr = gp.enrichr(
            gene_list=list(de_genes),
            gene_sets=[db],
            no_plot=True,
            cutoff=0.05  # test dataset, use lower value from range(0,1)
        )
        enr.results = enr.results[enr.results["Adjusted P-value"] < GSEA_THRESH]
        print(enr.results)
        sig_terms = sorted(set(enr.results['Term']))
        db_to_gene_sets[db] = sig_terms

    with open(out_f, 'w') as f:
        json.dump(db_to_gene_sets, f, indent=4)


if __name__ == '__main__':
    main()

