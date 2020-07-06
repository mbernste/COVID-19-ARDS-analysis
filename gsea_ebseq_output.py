import matplotlib as mpl
mpl.use('Agg')
from matplotlib_venn import venn2, venn2_circles, venn2_unweighted
from matplotlib_venn import venn3, venn3_circles
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
        #'MSigDB_Immunological_Signatures': '../gsvapy/gene_sets/c7.all.v7.1.symbols.gmt',
        'MSigDB_Hallmark': '../gsvapy/gene_sets/h.all.v7.1.symbols.gmt'
}
GSEA_THRESH = 0.05
TPM_MAX_THRESH = 10.0
FC_THRESH = 2.0
FILTER_BY_FOLD = True

def main():
    usage = "" # TODO
    parser = OptionParser(usage=usage)
    parser.add_option("-t", "--de_table_f", help="DE table file")
    parser.add_option("-c", "--de_table_col", help="DE table column to filter")
    parser.add_option("-o", "--out_file", help="Output file")
    parser.add_option("-k", "--out_kept", help="Output kept genes file")
    parser.add_option("-v", "--out_venn_f", help="Venn diagram file")
    (options, args) = parser.parse_args()

    de_genes_f = args[0]
    ec_f = args[1]
    out_f = options.out_file
    out_kept_f = options.out_kept
    out_venn_f = options.out_venn_f

    # Get all the DE genes
    with open(de_genes_f, 'r') as f:
        de_genes = [l.strip() for l in f]

    print('{} total DE genes'.format(len(de_genes)))

    # Map each gene to its fold change and filter
    # by fold change
    if FILTER_BY_FOLD:
        ec_df = pd.read_csv(ec_f, sep='\t', index_col=0)
        gene_to_fc = {
            gene: ec_df.loc[gene]['FC']
            for gene in ec_df.index
        }
        filtered_genes = [
            gene
            for gene in gene_to_fc
            if gene_to_fc[gene] > FC_THRESH or gene_to_fc[gene] < (1/FC_THRESH)
        ]
        print('{}/{} remain after filtering by fold-change.'.format(len(filtered_genes), len(gene_to_fc)))
        if len(filtered_genes) == 0:
            with open(out_f, 'w') as f:
                f.write('collection\tgene_set\tadjusted_p_value')
            return

    # Filter batch effect if specified in options
    if options.de_table_f:
        de_table_f = options.de_table_f
        de_table_col = options.de_table_col
        de_table_df = pd.read_csv(de_table_f, sep='\t', index_col=0)
        batch_effect_genes = set(de_table_df.loc[de_table_df[de_table_col] == 1][de_table_col].index) 
        
        # Output the removed genes
        present_batch_effect_genes = batch_effect_genes & set(filtered_genes)
        print('Removed {} genes: {}'.format(len(present_batch_effect_genes), present_batch_effect_genes))

        filtered_genes = sorted(set(filtered_genes) - batch_effect_genes)
        with open(out_kept_f, 'w') as f:
            f.write('\n'.join(sorted(filtered_genes)))

        print('{}/{} remain after filtering by batch-effect:'.format(len(filtered_genes), len(gene_to_fc)))
        print(filtered_genes)


        if len(filtered_genes) == 0:
            with open(out_f, 'w') as f:
                f.write('collection\tgene_set\tadjusted_p_value')
            return

    # Draw Venn diagram
    if 'Up' in de_genes_f:
        title = r'DE genes $\bf{higher}$ in COVID-19 ICU patients'
    else:
        title = r'DE genes $\bf{lower}$ in COVID-19 ICU patients'
    fig, ax = plt.subplots(1,1,figsize=(4.5,3.5))
    venn2(
        subsets = (
            len(filtered_genes),
            len(batch_effect_genes),
            len(present_batch_effect_genes),
        ), 
        set_labels = ('COVID-19 ICU\nvs. sepsis ARDS', 'non-ICU\nvs. sepsis non-ARDS'), set_colors=('b', 'y'), alpha = 0.5,
        ax=ax
    )
    plt.gca().set_title(title, fontsize=14)
    plt.tight_layout()
    plt.savefig(options.out_venn_f, format='pdf')

    # Perform gene set enrichment
    db_to_gene_sets = {}
    for db_name, db in GENE_SETS.items():
        enr = gp.enrichr(
            gene_list=[x.strip() for x in filtered_genes],
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
        db_to_gene_sets[db_name] = sig_terms


    # Create final dataframe
    da = []
    for db, gene_set_to_pval in db_to_gene_sets.items():
        for gene_set, pval in gene_set_to_pval.items():
            da.append((db, gene_set, pval))
    df = pd.DataFrame(
        data=da,
        columns=['collection', 'gene_set', 'adjusted_p_value']
    )
    df = df.sort_values(by='adjusted_p_value', axis=0)
    print('{} total enriched gene sets.'.format(len(df)))

    # Write output
    print('Writing to {}.'.format(out_f))
    df.to_csv(out_f, index=False, sep='\t')
    print('done')


if __name__ == '__main__':
    main()

