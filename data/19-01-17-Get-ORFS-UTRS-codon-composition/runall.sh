#!/bin/bash

mkdir -p sequence-data/
Rscript get_fish_cds_and_utr_data.R
Rscript compute_codon_composition.R
