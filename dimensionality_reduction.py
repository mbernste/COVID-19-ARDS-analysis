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
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis

def main():
    usage = "" # TODO 
    parser = OptionParser(usage=usage)
    parser.add_option("-s", "--subset", help="One of 'COVID', 'NONCOVID'")
    parser.add_option("-l", "--log", action="store_true", help="Take log1")
    parser.add_option("-o", "--output_prefix", help="Output file")
    (options, args) = parser.parse_args()

    expr_f = args[0]
    meta_f = args[1]
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
    print(meta_df)

    hospital_free = np.array(meta_df['Hospital_free_days'])
    hospital_free = np.nan_to_num(hospital_free)

    # Filter expression matrix according to metadata
    expr_df = expr_df[meta_df.index]
    X = np.array(expr_df)
    X = X.T
    if options.log:
        X = np.log(X+1)

    mod = PCA(n_components=2)
    X_pca = mod.fit_transform(X)
    _plot_scatter(
        X_pca, 
        'PCA', 
        hospital_free, 
        'Hospital Free Days',
        '{}.PCA_hospital_free.pdf'.format(prefix)
    )

    mod_tsne = TSNE(n_components=2, perplexity=6)
    mod_pca_100 = PCA(n_components=min([100, len(X)]))
    X_pca_100 = mod_pca_100.fit_transform(X)
    print(X_pca_100.shape)
    print('Fitting t-SNE...')
    X_tsne = mod_tsne.fit_transform(X_pca_100)
    print('done.')
    _plot_scatter(
        X_tsne,
        't-SNE',
        hospital_free,
        'Hospital Free Days',
        '{}.tSNE_hospital_free.pdf'.format(prefix)
    )

    s_hospital_free = sorted(set(hospital_free))
    one_third = int(len(s_hospital_free)/3)
    first_thresh = s_hospital_free[one_third]
    print(first_thresh)
    second_thresh = s_hospital_free[2*one_third]
    discrete_y = []
    for hf in hospital_free:
        if hf < first_thresh:
            discrete_y.append(1)
        elif hf >= first_thresh and hf < second_thresh:
            discrete_y.append(2)
        elif hf >= second_thresh:
            discrete_y.append(3)
    print(discrete_y)
    lda = LinearDiscriminantAnalysis(n_components=2)
    X_lda = lda.fit(X, discrete_y).transform(X)
    _plot_scatter(
        X_lda,
        'LDA',
        hospital_free,
        'Hospital Free Days',
        '{}.LDA_hospital_free.pdf'.format(prefix)
    )

def _plot_scatter(X, units, color_vals, color_variable, out_f):
    df = pd.DataFrame(
        [
            (d1, d2, col)
            for (d1, d2), col in zip(X, color_vals)
        ],
        columns=[
            '{} 1'.format(units),
            '{} 2'.format(units),
            color_variable
        ]
    )
    fig, ax = plt.subplots(1,1,figsize=(4,4))
    sc = ax.scatter(df['{} 1'.format(units)], df['{} 2'.format(units)], c=df[color_variable], cmap='viridis')
    plt.colorbar(sc)
    ax.set_xlabel('{} 1'.format(units))
    ax.set_ylabel('{} 2'.format(units))
    ax.set_title(color_variable)
    ax.set_xticks([])
    ax.set_yticks([])
    plt.tight_layout()
    fig.savefig(
        out_f,
        format='pdf',
        dpi=150
    )




if __name__ == "__main__":
    main()
