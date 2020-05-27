import numpy as np
import pandas as pd


sequences = pd.read_csv("../19-04-30-PredictiveModelDecayAllSpecies/19-04-30-EDA/results_data/training_set.csv")

# get only sequence data

sequences = sequences[['gene_id', 'coding']].drop_duplicates()

# generate a train and test set, test set are positve control

train = sequences.sample(frac=0.8, random_state=100)
test = sequences.drop(train.index)

train.to_csv("results_data/trainseqs.csv", index=False)
test.to_csv("results_data/testseqs.csv", index=False)
