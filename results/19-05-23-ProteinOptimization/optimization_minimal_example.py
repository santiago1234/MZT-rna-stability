# Minimal Example To run optimization program

import pandas as pd
from optimalcodon.projects.rnastability.optimization import protein_optimization

# params
# the path to the model file
mdl_path = '../19-04-30-PredictiveModelDecayAllSpecies/19-05-01-TrainModels/results_data/final_model_gbm.joblib'
# numbe of iteration to run
n_iterations = 50

# Gene name
gene_name = "mScarlet-I"
seq = "atggtgagcaagggcgaggcagtgatcaaggagttcatgcggttcaaggtgcacatggagggctccatgaacggccacgagttcgagatcgagggcgagggcgagggccgcccctacgagggcacccagaccgccaagctgaaggtgaccaagggtggccccctgcccttctcctgggacatcctgtcccctcagttcatgtacggctccagggccttcatcaagcaccccgccgacatccccgactactataagcagtccttccccgagggcttcaagtgggagcgcgtgatgaacttcgaggacggcggcgccgtgaccgtgacccaggacacctccctggaggacggcaccctgatctacaaggtgaagctccgcggcaccaacttccctcctgacggccccgtaatgcagaagaagacaatgggctgggaagcgtccaccgagcggttgtaccccgaggacggcgtgctgaagggcgacattaagatggccctgcgcctgaaggacggcggccgctacctggcggacttcaagaccacctacaaggccaagaagcccgtgcagatgcccggcgcctacaacgtcgaccgcaagttggacatcacctcccacaacgaggactacaccgtggtggaacagtacgaacgctccgagggccgccactccaccggcggcatggacgagctgtacaagtccggactcagatctcgagctcaagcttcgaattctgcagtcgacggtaccgcgggcccgggatccaccggatctagataa"

# make sure the sequence is upper case
seq = seq.upper()


# Parameter specification
mdl_params = {'specie': 'fish',  # the specie: (fish, human, mouse, xenopus)
              'cell_type': 'embryo mzt',  # cell type, this depends on the specie
              'datatype': 'aamanitin polya'}  # data type, depends on specie and cell type

gfp_path = protein_optimization(
    protein_sequence=seq, # Input user sequence
    mdl_path=mdl_path, # path to the model file
    mdl_params, # parameters for predictions
    crossover_args={}, # parameters for genetic algorithm
    niterations=n_iterations # number of iterations
)

# the log file is a pd.DataFrame

gfp_path.to_csv('results_data/mScarlet_optimization_path.csv', index=False)
