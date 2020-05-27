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
    if options.subset == 'COVID':
        meta_df = meta_df.loc[meta_df['COVID'] == 1]

    # Get metadata
    hospital_free = np.array(meta_df['Hospital_free_days'])
    icu_status = []
    for is_icu in meta_df['ICU_1']:
        if is_icu == 1:
            icu_status.append('True')
        else:
            icu_status.append('False')
    covid_status = []
    for is_covid in meta_df['COVID']:
        if is_covid == 1:
            covid_status.append('True')
        else:
            covid_status.append('False')

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
    _plot_scatter_discrete(
        X_pca,
        'PCA',
        icu_status,
        'ICU Status',
        '{}.PCA_icu_status.pdf'.format(prefix)
    )
    if options.subset is None:
        _plot_scatter_discrete(
            X_pca,
            'PCA',
            covid_status,
            'COVID-19\nStatus',
            '{}.PCA_covid_status.pdf'.format(prefix)
        )

    for perp in [5,6,7,8,9,10]:
        mod_tsne = TSNE(n_components=2, perplexity=perp)
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
            '{}.tSNE_perp_{}_hospital_free.pdf'.format(prefix, perp)
        )
        _plot_scatter_discrete(
            X_tsne,
            't-SNE',
            icu_status,
            'ICU Status',
            '{}.tSNE_perp_{}_icu_status.pdf'.format(prefix, perp)
        )
        if options.subset is None:
            _plot_scatter_discrete(
                X_tsne,
                't-SNE',
                covid_status,
                'COVID-19\nStatus',
                '{}.tSNE_perp_{}_covid_status.pdf'.format(prefix, perp)
            )

    # Linear discriminant analysis for hospital-free days
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
    sc = ax.scatter(
        df['{} 1'.format(units)], 
        df['{} 2'.format(units)], 
        c=df[color_variable], 
        cmap='viridis'
    )
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


def _plot_scatter_discrete(X, units, color_vals, color_variable, out_f):
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

    fig, ax = plt.subplots(1,1,figsize=(4,3))
    sns.scatterplot(
        x='{} 1'.format(units),
        y='{} 2'.format(units),
        hue=color_variable,
        data=df,
        ax=ax
    )
    ax.set_xticks([])
    ax.set_yticks([])
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    plt.tight_layout()
    fig.savefig(
        out_f,
        format='pdf',
        dpi=150
    )


if __name__ == "__main__":
    main()
