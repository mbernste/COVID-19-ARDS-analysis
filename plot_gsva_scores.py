import matplotlib as mpl
mpl.use('Agg')
from optparse import OptionParser
import seaborn as sns
from matplotlib import pyplot as plt
import sys
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
import numpy as np
from os.path import join
import json

def main():
    usage = "" # TODO 
    parser = OptionParser(usage=usage)
    parser.add_option("-o", "--out_dir", help="Output directory")
    (options, args) = parser.parse_args()

    expr_f = args[0]
    meta_f = args[1]
    gene_sets_f = args[2]
    out_dir = options.out_dir

    expr_df = pd.read_csv(expr_f, sep='\t', index_col=0)
    meta_df = pd.read_csv(meta_f, sep='\t')
    meta_df = meta_df.set_index('Albany_sampleID')

    # Remove patients for which the 28 days has not elapsed
    meta_df = meta_df.loc[meta_df['Hospital_free_days'].notnull()]
    no_expression_data = set(meta_df.index) - set(expr_df.columns)
    meta_df = meta_df.drop(no_expression_data)

    # Rmoeve non-COVID patients
    meta_df = meta_df.loc[meta_df['COVID'] == 1]

    hospital_free = np.array(meta_df['Hospital_free_days'])
    hospital_free = np.nan_to_num(hospital_free)

    expr_df = expr_df[meta_df.index]
    expr_df = expr_df.transpose()
    expr_df['Hospital_free_days'] = hospital_free
    expr_df = expr_df.sort_values(by='Hospital_free_days')

    print(expr_df)

    # Load the enriched gene-sets
    if gene_sets_f.split('.')[-1] == 'json':
        with open(gene_sets_f, 'r') as f:
            db_to_gene_sets = json.load(f)
        all_gene_sets = set()
        for gene_sets in db_to_gene_sets.values():
            all_gene_sets.update(gene_sets)
    elif gene_sets_f.split('.')[-1] == 'tsv':
        with open(gene_sets_f, 'r') as f:
            all_gene_sets = [l.strip() for l in f]

    # Plot bargraphs for each enriched gene-set
    for gene_set in all_gene_sets:
        if gene_set not in expr_df.columns:
            print('Skipping gene set {}. Not found...'.format(gene_set))
            continue
        print('Plotting gene set {}.'.format(gene_set))
        plot_df = pd.DataFrame(
            data=[
                (str(int(hf)), score)
                for hf, score in zip(expr_df['Hospital_free_days'], expr_df[gene_set])
            ],
            columns=['Hospital_free_days', 'Enrichment score']
        )

        # Bar graph
        fig, ax = plt.subplots(1,1,figsize=(0.1*len(hospital_free),3))
        sns.barplot(x='Hospital_free_days', y=gene_set, color="#0052cc", data=expr_df, ax=ax, ci='sd')
        ax.set_title(_prettify_label(gene_set))
        ax.set_ylabel('Enrichment Score')
        ax.set_xlabel('Hospital-free Days')
        plt.tight_layout()
        fig.savefig(
            join(out_dir, '{}.bar.pdf'.format(gene_set)),
            format='pdf',
            dpi=150
        )

        fig, ax = plt.subplots(1,1,figsize=(3,3))
        sns.regplot(y='Hospital_free_days', x=gene_set, color="#0052cc", data=expr_df, ax=ax)
        ax.set_title(_prettify_label(gene_set))
        ax.set_ylabel('Enrichment Score')
        ax.set_xlabel('Hospital-free Days')
        plt.tight_layout()
        fig.savefig(
            join(out_dir, '{}.reg.pdf'.format(gene_set)),
            format='pdf',
            dpi=150
        )

def _prettify_label(gene_set):
    toks = gene_set.split('_')
    mod_toks = []
    for tok_i, tok in enumerate(toks):
        if tok_i == 0 and tok == 'GO':
            continue
        else:
            tok = tok[0].upper() + tok[1:].lower()
            mod_toks.append(tok)
    return ' '.join(mod_toks)

if __name__ == "__main__":
    main()
