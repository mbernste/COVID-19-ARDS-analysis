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
from sklearn.linear_model import ElasticNet, SGDRegressor
from sklearn.model_selection import GridSearchCV
from sklearn import metrics
import json

PARAMETERS = {
    'alpha': [
        0.001, 
        0.01, 
        0.1, 
        1.0,
        10.0,
        100.0
    ],
    'l1_ratio': [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
}

MAX_ITER = 100000

def main():
    usage = "" # TODO 
    parser = OptionParser(usage=usage)
    parser.add_option("-s", "--subset", help="One of 'COVID', 'NONCOVID'")
    parser.add_option("-l", "--log", action="store_true", help="Take log1 of data")
    parser.add_option("-o", "--output_prefix", help="Output file")
    (options, args) = parser.parse_args()

    expr_f = args[0]
    meta_f = args[1]
    prefix = options.output_prefix

    expr_df = pd.read_csv(expr_f, sep='\t', index_col=0)
    meta_df = pd.read_csv(meta_f, sep='\t')
    meta_df = meta_df.set_index('Albany_sampleID')
    print(meta_df)

    # Remove patients for which the 28 days has not elapsed
    meta_df = meta_df.loc[meta_df['Hospital_free_days'].notnull()]
    no_expression_data = set(meta_df.index) - set(expr_df.columns)
    meta_df = meta_df.drop(no_expression_data)

    # Rmoeve non-COVID patients
    if options.subset == 'COVID':
        meta_df = meta_df.loc[meta_df['COVID'] == 1]
    elif options.subset == 'NONCOVID':
        meta_df = meta_df.loc[meta_df['COVID'] == 0]
        
    print(meta_df)

    hospital_free = np.array(meta_df['Hospital_free_days'])
    hospital_free = np.nan_to_num(hospital_free)

    # Filter expression matrix according to metadata
    expr_df = expr_df[meta_df.index]
    X = np.array(expr_df)
    X = X.T
    if options.log:
        X = np.log(X+1)

    model = ElasticNet(max_iter=MAX_ITER)
    #model = SGDRegressor(penalty='elasticnet')
    clf = GridSearchCV(model, PARAMETERS, scoring='neg_mean_squared_error')
    clf.fit(X, hospital_free)
    print(sorted(clf.cv_results_.keys()))
    print(clf.cv_results_['param_alpha'])
    print(clf.cv_results_['mean_test_score'])

    best_params = max(
        zip(clf.cv_results_['param_alpha'], clf.cv_results_['param_l1_ratio'], clf.cv_results_['mean_test_score']),
        key=lambda x: x[2]
    )
    best_alpha = best_params[0]
    best_ratio = best_params[1]
    print('Max (alpha, ratio): ({}, {})'.format(best_alpha, best_ratio))

    model = ElasticNet(alpha=best_alpha, l1_ratio=best_ratio, max_iter=MAX_ITER)
    #model = SGDRegressor(penalty='elasticnet', alpha=best_alpha, l1_ratio=best_ratio)
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
