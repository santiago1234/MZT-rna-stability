import numpy as np
import pandas as pd
from sklearn.linear_model import lasso_path
from lassoloaddata import get_folds

# define the grid of lambda values to explore
alphas = np.logspace(-4, -0.5, 30)


def get_lasso_path(X_train, y_train, alphas=alphas):
    """
    compute the lasso path for the given data
    Args:
        X_train: predictors k-mers and PLS, assumes that predictros
        have being scaled, numpy array
        y_train: response var (mRNA indicator) pd.Series or array
    """
    alphas_lasso, coefs_lasso, _ = lasso_path(X_train, y_train, alphas=alphas,
                                              fit_intercept=False)
    return alphas_lasso, coefs_lasso


def path_to_frame(coefs, colnames, rownames):
    """
    puts the lasso path coefs in a pandas
    data frame with adecuate rownames and
    column names
    Args:
        colnames: list, explored lambda values
        rownames: list, name of predicotr variables
    """
    return pd.DataFrame(coefs, index=rownames, columns=colnames)


def lasssopath(predictors, response, alphas=alphas):
    """
    computes the lasso path on the data using 6 folds
    of the data, therefore returns 6 different paths
    Args:
        predictors: predictors pd.DataFrame
        response: response pd.Series
    """
    folds = get_folds(predictors, response,
                      k=6)  # get folds and scales predictors
    paths = []
    i = 1
    for (X_train, y_train), _ in folds:
        print('runing fold {} of 6'.format(i))
        alphas, lasso_path=get_lasso_path(X_train, y_train, alphas=alphas)
        coefs=path_to_frame(lasso_path, alphas, predictors.columns)
        coefs['kfold']=i
        paths.append(coefs)
        i += 1
    return pd.concat(paths)
