#! /usr/bin/env python
import argparse
import pandas as pd
from optimalcodon.projects.rnastability.optimization import protein_optimization


mdl_path = "../191005-EvaluateModelLearningCurve/results-data/final_model.joblib"

mdl_params = {'specie': 'human',  # the specie: (fish, human, mouse, xenopus)
              'cell_type': 'hela',  # cell type, this depends on the specie
              'datatype': 'endogenous'}  # data type, depends on speci

n_iterations = 50


def run_optimization(sequence):
    """helper function to run the optimization"""
    search_path = protein_optimization(
        protein_sequence=sequence,  # Input user sequence
        mdl_path=mdl_path,  # path to the model file
        mdl_args=mdl_params,  # parameters for predictions
        crossover_args={},  # parameters for genetic algorithm
        niterations=n_iterations  # number of iterations
    )

    return search_path


def run(args):
    seq = args.seq
    geneid = args.geneid

    search = run_optimization(seq)
    search['gene_id'] = geneid
    search.to_csv(geneid + '-optimization.csv', index=False)


def main():
    parser = argparse.ArgumentParser(
        description='Optimize a mRNA by sysnonimous mutations')
    parser.add_argument("-s", help="coding input sequence",
                        dest="seq", type=str, required=True)
    parser.add_argument("-id", help="id sequence name",
                        dest="geneid", type=str, required=True)
    parser.set_defaults(func=run)
    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
