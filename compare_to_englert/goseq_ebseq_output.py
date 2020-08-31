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

import run_goseq

#GENE_SETS = ['GO_Biological_Process_2018', 'GO_Molecular_Function_2018', 'KEGG_2019_Human']
#GENE_SETS = ['/Users/matthewbernstein/Development//gsvapy/gene_sets/c5.bp.v7.1.symbols.gmt']
GENE_SETS = {
        'GO_molecular_function': '../gene_sets/c5.mf.v7.1.symbols.gmt',
        'GO_biological_process': '../gene_sets/c5.bp.v7.1.symbols.gmt',
        'MSigDB_Canonical_Pathways': '../gene_sets/c2.cp.v7.1.symbols.gmt',
        'MSigDB_Hallmark': '../gene_sets/h.all.v7.1.symbols.gmt'
}
GSEA_THRESH = 0.05
FC_THRESH = 2.0
FILTER_BY_FOLD = True
NUMBER_OF_TOTAL_GENES = 19462 

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
    raw_f = args[2]
    out_f = options.out_file
    out_kept_f = options.out_kept
    out_venn_f = options.out_venn_f

    # Get all the DE genes
    with open(de_genes_f, 'r') as f:
        de_genes = [l.strip() for l in f]
    print('{} total DE genes'.format(len(de_genes)))

    # Parse all genes
    raw_df = pd.read_csv(raw_f, sep='\t', index_col=0)
    all_genes = list(raw_df.index)
    print('{} total genes'.format(len(all_genes)))

    # Map each gene to its fold change and filter
    # by fold change
    ec_df = pd.read_csv(ec_f, sep='\t', index_col=0)
    gene_to_fc = {
        gene: ec_df.loc[gene]['FC'] 
        for gene in ec_df.index
    }
    fold_filtered_genes = [
        gene
        for gene in gene_to_fc
        if gene_to_fc[gene] > FC_THRESH or gene_to_fc[gene] < (1/FC_THRESH)
    ]
    print('{}/{} remain after filtering by fold-change.'.format(len(fold_filtered_genes), len(gene_to_fc)))
    if len(fold_filtered_genes) == 0:
        with open(out_f, 'w') as f:
            f.write('collection\tgene_set\tadjusted_p_value')
        return

    # Filter batch effect if specified in options
    if options.de_table_f:
        de_table_f = options.de_table_f
        de_table_col = options.de_table_col
        de_table_df = pd.read_csv(de_table_f, sep='\t', index_col=0)
        batch_effect_genes = frozenset(
            de_table_df.loc[de_table_df[de_table_col] == 1].index
        ) 
        
        # Output the removed genes
        removed_batch_effect_genes = batch_effect_genes & set(fold_filtered_genes)
        print('Removed {} genes: {}'.format(
            len(removed_batch_effect_genes), 
            removed_batch_effect_genes
        ))

        fold_batch_filtered_genes = sorted(set(fold_filtered_genes) - batch_effect_genes)

        # Write genes that we keep to a file
        with open(out_kept_f, 'w') as f:
            f.write('\n'.join(sorted(fold_batch_filtered_genes)))
        print('{}/{} remain after filtering by batch-effect:'.format(
            len(fold_batch_filtered_genes), 
            len(gene_to_fc)
        ))
        
        if len(fold_batch_filtered_genes) == 0:
            with open(out_f, 'w') as f:
                f.write('collection\tgene_set\tadjusted_p_value')
            return

    # Draw Venn diagram
    if 'Up' in de_genes_f:
        title = r'DE genes $\bf{higher}$ in COVID-19 ARDS patients'
    else:
        title = r'DE genes $\bf{lower}$ in COVID-19 ARDS patients'
    fig, ax = plt.subplots(1,1,figsize=(5.5,3.5))
    venn2(
        subsets = (
            len(fold_batch_filtered_genes),
            len(batch_effect_genes),
            len(removed_batch_effect_genes),
        ), 
        set_labels = (
            'COVID-19 ARDS\nvs. non-COVID-19 ARDS', 
            'non-COVID-19 non-ARDS\n(this study)\nvs. non-COVID-19 non-ARDS\n(Englert et al.)'
        ), 
        set_colors=('b', 'y'), 
        alpha = 0.5,
        ax=ax
    )
    plt.gca().set_title(title, fontsize=15)
    plt.tight_layout()
    plt.savefig(options.out_venn_f, format='pdf')

    # Perform gene set enrichment
    db_to_gene_sets = {}
    for db_name, db in GENE_SETS.items():
        gene_set_to_genes = run_goseq._parse_gene_sets(db)

        # Map each gene to its gene sets
        gene_to_gene_sets = defaultdict(lambda: [])
        for gene_set, genes in gene_set_to_genes.items():
            for gene in genes:
                gene_to_gene_sets[gene].append(gene_set)
        gene_to_gene_sets = dict(gene_to_gene_sets)

        df_goseq = run_goseq.run_GOseq(
            fold_batch_filtered_genes, 
            all_genes,
            gene_to_gene_sets
        )
        df_goseq = df_goseq.loc[df_goseq['fdr'] < GSEA_THRESH]
        gene_set_to_pval = {
            gene_set: df_goseq.loc[gene_set]['fdr']
            for gene_set in df_goseq.index
        }
        db_to_gene_sets[db_name] = gene_set_to_pval

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

