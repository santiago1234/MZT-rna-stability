from os import path
import pandas as pd
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import KFold

# this variables set the path they may change depending on
# the location where the module is run
predictorspathtodata = "results_data/"
mrnastabilitydpath = "../19-03-17-ObatainStabilityMeasureInWTConditions/results_data/wt_mRNA_stability_estimate.csv"


def loadpredictors(predictorspathtodata=predictorspathtodata):
    """loads the train and test data, merges into a single file"""
    datafiles = ("kpreds_train.csv", "kpreds_test.csv")
    predictors = [pd.read_csv(path.join(predictorspathtodata, x))
                  for x in datafiles]
    return pd.concat(predictors)


def loadresponsedata(dpath=mrnastabilitydpath):
    """loads the response variable, the mrna
    stability indicator"""
    response = pd.read_csv(dpath)
    # keep the columns we need and reshape
    response = response[['Gene_ID', 'sample_condition', 'estimate']]
    return response.pivot(index='Gene_ID', columns='sample_condition',
                          values='estimate')


def loadlassodata(condition, predictorspathtodata=predictorspathtodata,
                  dpathresponse=mrnastabilitydpath):
    """gets the data to fit a ml model"""
    predictors = loadpredictors(predictorspathtodata)
    response = loadresponsedata(dpathresponse)
    data = pd.merge(
        response[[condition]],
        predictors,
        left_index=True,
        right_on='Gene_ID',
        how='inner')
    # drop the gene id var
    data = data.drop(['Gene_ID'], axis=1)
    y = data[condition]
    X = data.drop([condition], axis=1)
    return X, y


def preprocess(X):
    """
    preporcess the data with StandardScaler, return the scaler
    """
    scaler = StandardScaler()
    scaler.fit(X)
    return scaler


def get_folds(X, y, k=5):
    """
    yields the data folds, scaling is applied to X
    """
    kf = KFold(n_splits=k)
    for train_index, test_index in kf.split(X):
        X_train, X_test = X.iloc[train_index], X.iloc[test_index]
        y_train, y_test = y.iloc[train_index], y.iloc[test_index]
        # scale the predictors
        scaler = preprocess(X_train)
        X_train = scaler.transform(X_train)
        X_test = scaler.transform(X_test)
        yield (X_train, y_train), (X_test, y_test)
