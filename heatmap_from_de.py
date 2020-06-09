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
    parser.add_option("-a", "--all_genes", action='store_true', help="Don't filter for DE genes")
    parser.add_option("-u", "--up_genes", action='store_true', help="Filter DE-up genes")
    parser.add_option("-d", "--down_genes", action='store_true', help="Filter DE-down genes")
    parser.add_option("-o", "--out_dir", help="Output directory")
    (options, args) = parser.parse_args()

    gene_set_f = args[0]
    gene_set = args[1]
    up_de_list_f = args[2]
    down_de_list_f = args[3]
    albany_meta_f = args[4]
    englert_meta_f = args[5]
    albany_expr_f = args[6]
    englert_expr_f = args[7]
    out_dir = options.out_dir

    # Load gene set
    gene_set_to_genes = _parse_gmt(gene_set_f)
    genes = gene_set_to_genes[gene_set]
    print('There are {} genes in the gene set {}'.format(len(genes), gene_set))

    # Load up-DE genes
    all_de_genes = set()
    gene_to_fc = {}
    df = pd.read_csv(up_de_list_f, sep='\t', index_col=0)
    all_de_genes.update(df.index)
    up_de_genes = set(df.index)
    for gene in df.index:
        gene_to_fc[gene] = abs(df.loc[gene]['L2FC'])   

    # Load down-DE genes
    df = pd.read_csv(down_de_list_f, sep='\t', index_col=0)
    all_de_genes.update(df.index)
    down_de_genes = set(df.index)
    for gene in df.index:
        gene_to_fc[gene] = abs(df.loc[gene]['L2FC']) 

    # Load the data
    print('Loading expression data...')
    albany_expr_df = pd.read_csv(albany_expr_f, sep='\t', index_col=0)
    albany_meta_df = pd.read_csv(albany_meta_f, sep='\t', index_col=0)
    englert_expr_df = pd.read_csv(englert_expr_f, sep='\t', index_col=0)
    englert_meta_df = pd.read_csv(englert_meta_f, sep='\t', index_col=0)

    # Remove genes not in expression data
    genes = set(genes) & set(englert_expr_df.index) & set(albany_expr_df.index)

    # Filter for DE genes
    if options.up_genes:
        genes = set(genes) & up_de_genes
        genes = sorted(genes, key=lambda x: gene_to_fc[x], reverse=True)
    elif options.down_genes: 
        genes = set(genes) & down_de_genes
        genes = sorted(genes, key=lambda x: gene_to_fc[x], reverse=True)
    elif options.up:
        genes = sorted(genes)

    # Filter only COVID ICU patients
    albany_meta_df = albany_meta_df.set_index('Albany_sampleID')
    albany_meta_df = albany_meta_df.loc[list(albany_expr_df.columns)]
    keep_patients = albany_meta_df.loc[(albany_meta_df['COVID'] == 1) & (albany_meta_df['ICU_1'])].index
    albany_expr_df = albany_expr_df[keep_patients]

    # Filter only ARDS Englert patients
    keep_patients = englert_meta_df.loc[englert_meta_df['condition'] == 'ARDS'].index
    englert_expr_df = englert_expr_df[keep_patients]

    # Keep only gene set genes
    print('Restricting to gene sets...')
    englert_expr_df = englert_expr_df.loc[genes]
    albany_expr_df = albany_expr_df.loc[genes]

    # Align their indices
    albany_expr_df = albany_expr_df.loc[englert_expr_df.index]
    albany_meta_df = albany_meta_df.loc[albany_expr_df.columns]

    # Get hospital free days
    hospital_free = albany_meta_df['Hospital_free_days']

    # Take log1
    albany_X = np.array(albany_expr_df)
    englert_X = np.array(englert_expr_df)
    print('Shape of Englert matrix: ', englert_X.shape)
    print('Shape of Albany matrix: ', albany_X.shape)
    
    # Compute z-scores
    X = np.hstack([albany_X, englert_X])
    print(X.shape)
    means = np.expand_dims(np.mean(X, axis=1), axis=0)
    sd = np.expand_dims(np.std(X, axis=1), axis=0)
    albany_X_zscore = (albany_X - means.T) / sd.T 
    englert_X_zscore = (englert_X - means.T) / sd.T 
    albany_X_zscore = np.nan_to_num(albany_X_zscore)
    englert_X_zscore = np.nan_to_num(englert_X_zscore)
    albany_X_zscore[albany_X_zscore == np.inf] = 0
    englert_X_zscore[englert_X_zscore == np.inf] = 0

    # Take log1
    albany_X_log = np.log(albany_X+1)
    englert_X_log = np.log(englert_X+1)

    # Remove genes with zero-mean
    keep_inds = [
        mean_i
        for mean_i, mean in enumerate(np.squeeze(means))
        if mean > 0
    ]
    albany_X_zscore = albany_X_zscore[keep_inds,:]
    englert_X_zscore = englert_X_zscore[keep_inds,:]
    albany_X_log = albany_X_log[keep_inds,:]
    englert_X_log = englert_X_log[keep_inds,:]
    genes = np.array(genes)[keep_inds]

    print(len(genes))

    # Separated heatmap
    fig, axarr = plt.subplots(1,2, figsize=(9,0.2*englert_X_log.shape[0]), gridspec_kw={'width_ratios': (englert_X_log.shape[1], albany_X_log.shape[1])})
    # Re-order the columns according to Seaborns hierarchical clustering
    cg = sns.clustermap(np.hstack([englert_X_log, albany_X_log]))
    row_inds = cg.dendrogram_row.reordered_ind
    gene_labels = np.array(genes)[row_inds]
    cg = sns.clustermap(englert_X_log, row_cluster=True, col_cluster=True)
    englert_col_inds = cg.dendrogram_col.reordered_ind
    englert_X_log = englert_X_log[:,englert_col_inds]
    englert_X_log = englert_X_log[row_inds,:]
    cg =  sns.clustermap(albany_X_log, row_cluster=True, col_cluster=True)
    albany_col_inds = cg.dendrogram_col.reordered_ind
    albany_X_log = albany_X_log[:,albany_col_inds]
    albany_X_log = albany_X_log[row_inds,:]
    # Generate heatmaps
    sns.heatmap(englert_X_log, ax=axarr[0], cbar=False, cmap='viridis', yticklabels=gene_labels, xticklabels=False, cbar_kws={'label': 'log(TPM+1)'})
    sns.heatmap(albany_X_log, ax=axarr[1], cmap='viridis', yticklabels=False, xticklabels=False, cbar_kws={'label': 'log(TPM+1)'})
    axarr[0].set_title('Non-COVID-19 ARDS')
    axarr[1].set_title('COVID-19 ICU')
    plt.tight_layout()
    fig.savefig(
        join(out_dir, '{}.heatmap.pdf'.format(gene_set)),
        format='pdf',
        dpi=150,
        bbox_inches='tight'
    )

    # Genreate z-score
    fig, axarr = plt.subplots(1,2, figsize=(9,0.2*englert_X_zscore.shape[0]), gridspec_kw={'width_ratios': (englert_X_zscore.shape[1], albany_X_zscore.shape[1])})
    
    # Cluster the genes using all samples
    cg = sns.clustermap(np.hstack([englert_X_zscore, albany_X_zscore]))
    row_inds = cg.dendrogram_row.reordered_ind
    gene_labels = np.array(genes)[row_inds]

    # Cluster the Englert columns
    cg = sns.clustermap(englert_X_zscore, row_cluster=True, col_cluster=True)
    englert_col_inds = cg.dendrogram_col.reordered_ind
    englert_X_zscore = englert_X_zscore[:,englert_col_inds]
    englert_X_zscore = englert_X_zscore[row_inds,:]
    # Prepend an empty row to match up with Albany
    englert_X_zscore = np.vstack([np.zeros(englert_X_zscore.shape[1]), englert_X_zscore])
    
    # Cluster the Albany columns
    cg =  sns.clustermap(albany_X_zscore, row_cluster=True, col_cluster=True)
    albany_col_inds = cg.dendrogram_col.reordered_ind
    albany_X_zscore = albany_X_zscore[:,albany_col_inds]
    albany_X_zscore = albany_X_zscore[row_inds,:]
    max_val = np.max(albany_X_zscore)
    min_val = np.min(albany_X_zscore)
    hospital_free = np.array(hospital_free)[albany_col_inds]
    max_hval = np.max(hospital_free)
    min_hval = 0.0
    print('The minimum hospital free days is: ', np.min(hospital_free))
    print('The maximum hospital free days is: ', np.max(hospital_free))


    # Create hospital free days row
    print(albany_X_zscore.shape)
    albany_X_zscore = np.vstack([hospital_free, albany_X_zscore])
    print(albany_X_zscore.shape)

    hospital_free_mask = []
    zscore_mask = []
    for row_i, row in enumerate(albany_X_zscore):
        if row_i == 0:
            h_mask_row = [False for x in row]
            z_mask_row = [True for x in row]
        else:
            h_mask_row = [True for x in row]
            z_mask_row = [False for x in row]
        hospital_free_mask.append(h_mask_row)
        zscore_mask.append(z_mask_row)

    hospital_free_mask = np.array(hospital_free_mask)
    zscore_mask = np.array(zscore_mask)

    sns.heatmap(
        englert_X_zscore, 
        ax=axarr[0], 
        cbar=False, 
        cmap='bwr', 
        center=0, 
        yticklabels=['Hospital-free'] + list(gene_labels), 
        xticklabels=False
    )
    sns.heatmap(albany_X_zscore, ax=axarr[1], cmap='bwr', center=0, yticklabels=False, xticklabels=False, cbar_kws={'label': 'z-score'}, mask=zscore_mask, vmin=-3.0, vmax=3.0)
    sns.heatmap(albany_X_zscore, ax=axarr[1], cmap='viridis', yticklabels=False, xticklabels=False, cbar=False, mask=hospital_free_mask, vmin=min_hval, vmax=max_hval)
    axarr[0].set_title('Non-COVID-19 ARDS')
    axarr[1].set_title('COVID-19 ICU')

    if options.all_genes:
        for gene, k in zip(gene_labels, axarr[0].get_ymajorticklabels()):
            if gene in up_de_genes:
                k.set_color('#F9AF10')
                #k.set_label('{} HAHAHAH $/uparrow$'.format(gene))
            elif gene in down_de_genes:
                k.set_color('#1088F9')


    plt.tight_layout()
    fig.savefig(
        join(out_dir, '{}.heatmap_zscore.pdf'.format(gene_set)),
        format='pdf',
        dpi=150,
        bbox_inches='tight'
    )


def _parse_gmt(gmt_f):
    gene_set_to_genes = {}
    with open(gmt_f, 'r') as f:
        for l in f:
            toks = l.split()
            gene_set = toks[0]
            genes = toks[2:]
            gene_set_to_genes[gene_set] = genes
    return gene_set_to_genes

if __name__ == '__main__':
    main()
