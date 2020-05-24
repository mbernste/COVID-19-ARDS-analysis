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
from sklearn.linear_model import Lasso
from sklearn.model_selection import GridSearchCV
from sklearn import metrics
import json

PARAMETERS = {
    'alpha': [
        0.0000001,
        0.000001,
        0.00001, 
        0.0001, 
        0.001, 
        0.01, 
        0.1, 
        1.0, 
        10.0, 
        100.0, 
        1000.0, 
        10000.0
    ]
}

def main():
    usage = "" # TODO 
    parser = OptionParser(usage=usage)
    parser.add_option("-s", "--subset", help="One of 'COVID', 'NONCOVID'")
    parser.add_option("-o", "--output_prefix", help="Output file")
    (options, args) = parser.parse_args()

    expr_f = args[0]
    meta_f = args[1]
    prefix = options.output_prefix

    # Load the dta and metadata
    expr_df = pd.read_csv(expr_f, sep='\t', index_col=0)
    meta_df = pd.read_csv(meta_f, sep='\t')
    meta_df = meta_df.set_index('Albany_sampleID')

    # Remove patients for which the 28 days has not elapsed
    meta_df = meta_df.loc[meta_df['Hospital_free_days'].notnull()]
    no_expression_data = set(meta_df.index) - set(expr_df.columns)
    meta_df = meta_df.drop(no_expression_data)

    # Remove non-COVID patients
    if options.subset == 'COVID':
        meta_df = meta_df.loc[meta_df['COVID'] == 1]
    elif options.subset == 'NONCOVID':
        meta_df = meta_df.loc[meta_df['COVID'] == 0]
    print(meta_df)

    # Convert to numpy array
    hospital_free = np.array(meta_df['Hospital_free_days'])

    # Filter expression matrix according to metadata
    expr_df = expr_df[meta_df.index]
    X = np.array(expr_df)
    X = X.T
    #if options.log:
    #    X = np.log(X+1)

    # Cross-fold validation grid-search across lambda parameters
    model = Lasso(max_iter=20000)
    clf = GridSearchCV(
        model, 
        PARAMETERS, 
        scoring='neg_mean_squared_error'
    )
    clf.fit(X, hospital_free)
    print(sorted(clf.cv_results_.keys()))
    print(clf.cv_results_['param_alpha'])
    print(clf.cv_results_['mean_test_score'])

    best_alpha = max(
        zip(clf.cv_results_['param_alpha'], clf.cv_results_['mean_test_score']),
        key=lambda x: x[1]
    )[0]
    print('Max alpha: {}'.format(best_alpha))

    model = Lasso(alpha=best_alpha, max_iter=20000)
    #model = Lasso(alpha=0.01, max_iter=20000)
    model.fit(X, hospital_free)
    coeffs = model.coef_

    non_zero_feats = [
        feat 
        for feat, coef in zip(expr_df.index, coeffs)
        if coef > 0.0
    ]
    print("{} non-zero features.".format(len(non_zero_feats)))
   
   
    with open('{}.kept_features.tsv'.format(prefix), 'w') as f:
        f.write('\n'.join(non_zero_feats))

    preds = model.predict(X)
    mae = metrics.mean_absolute_error(hospital_free, preds)
    with open('{}.model_fit.json'.format(prefix), 'w') as f:
        json.dump(
            {
                'Mean absolute error': mae
            },
            f,
            indent=True
        )


if __name__ == "__main__":
    main()
