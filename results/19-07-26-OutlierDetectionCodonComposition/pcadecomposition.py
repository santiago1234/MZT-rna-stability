import numpy as np
import pandas as pd

from sklearn.base import BaseEstimator, TransformerMixin
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline

from optimalcodon.codons import codon_composition
from optimalcodon.projects.rnastability.dataprocessing import CodonCompositionTransformer

train = pd.read_csv("results_data/trainseqs.csv")

class CodonCompositionTransformer(BaseEstimator, TransformerMixin):
    """
    get codon composition of coding data,
    input data should be pd.DataFeame
    and should include column named coding
    for the coding DNA sequence
    """

    def fit(self, X, y=None):
        return self

    def transform(self, X, y=None):
        assert 'coding' in X.columns, "not column named coding in data for cds seqs"
        counts = codon_composition(X.coding).values
        counts = np.log(counts + 1)
        return counts


pipeline = Pipeline([('codon composition', CodonCompositionTransformer()),
                          ('scaler', StandardScaler())])

