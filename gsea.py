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


#GENE_SETS = ['GO_Biological_Process_2018', 'GO_Molecular_Function_2018', 'KEGG_2019_Human']
#GENE_SETS = ['/Users/matthewbernstein/Development//gsvapy/gene_sets/c5.bp.v7.1.symbols.gmt']
GENE_SETS = {
        'GO_molecular_function': '../gsvapy/gene_sets/c5.mf.v7.1.symbols.gmt',
        'GO_biological_process': '../gsvapy/gene_sets/c5.bp.v7.1.symbols.gmt',
        'MSigDB_Canonical_Pathways': '../gsvapy/gene_sets/c2.cp.v7.1.symbols.gmt',
        'MSigDB_Immunological_Signatures': '../gsvapy/gene_sets/c7.all.v7.1.symbols.gmt',
        'MSigDB_Hallmark': '../gsvapy/gene_sets/h.all.v7.1.symbols.gmt'
}
GSEA_THRESH = 0.05

def main():
    usage = "" # TODO
    parser = OptionParser(usage=usage)
    parser.add_option("-o", "--out_file", help="OUTPUT file")
    (options, args) = parser.parse_args()

    de_genes_fs = args[0].split(',')
    out_f = options.out_file

    de_genes = set()
    for de_genes_f in de_genes_fs:
        with open(de_genes_f, 'r') as f:
            genes = [l for l in f]
        de_genes.update(genes)

    db_to_gene_sets = {}
    for db_name, db in GENE_SETS.items():
        enr = gp.enrichr(
            gene_list=[x.strip() for x in de_genes],
            gene_sets=[db],
            background=19463,
            no_plot=True,
            cutoff=0.05  # test dataset, use lower value from range(0,1)
        )
        enr.results = enr.results[enr.results["Adjusted P-value"] < GSEA_THRESH]
        sig_terms = {
            str(row[0]): float(row[1])
            for row_i, row in enr.results[['Term', 'Adjusted P-value']].iterrows()
        }
        print(sig_terms)
        db_to_gene_sets[db_name] = sig_terms

    with open(out_f, 'w') as f:
        json.dump(db_to_gene_sets, f, indent=4)


if __name__ == '__main__':
    main()

