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

def main():
    usage = "" # TODO 
    parser = OptionParser(usage=usage)
    parser.add_option("-s", "--subset", help="One of 'COVID', 'NONCOVID'")
    parser.add_option("-l", "--log", action="store_true", help="Take log1")
    parser.add_option("-o", "--output_prefix", help="Output file")
    (options, args) = parser.parse_args()

    expr_f = args[0]
    meta_f = args[1]
    gene_set = args[2] #'GO_MONOCYTE_AGGREGATION'
    prefix = options.output_prefix

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

    plot_df = pd.DataFrame(
        data=[
            (str(int(hf)), score)
            for hf, score in zip(expr_df['Hospital_free_days'], expr_df[gene_set])
        ],
        columns=['Hospital_free_days', 'Enrichment score']
    )

    fig, ax = plt.subplots(1,1,figsize=(0.1*len(hospital_free),3))
    #sns.barplot(x=plot_df.index, y=plot_df['Enrichment score'], color="#0052cc", ax=ax, ci=None)
    sns.barplot(x='Hospital_free_days', y=gene_set, color="#0052cc", data=expr_df, ax=ax, ci='sd')
    #ax.set_xticklabels(plot_df['Hospital_free_days'])

    plt.tight_layout()
    fig.savefig(
        'test.pdf',
        format='pdf',
        dpi=150
    )




if __name__ == "__main__":
    main()
