#!/bin/bash
set -e
set -u
set -o pipefail


# script to get the data from LIMS
# make dir to put the files
mkdir -p rna-seq-profiles

# get ribo profile
wget http://bioinfo/n/core/Bioinformatics/secondary/Bazzini/arb/MOLNG-2541/secundo/RSEM_TPM_table.csv
mv RSEM_TPM_table.csv rna-seq-profiles/ribo-zero-profile.csv


# poly-A profile
wget http://bioinfo/n/core/Bioinformatics/secondary/Bazzini/arb/MOLNG-2540/secundo/RSEM_TPM_table.csv
mv RSEM_TPM_table.csv rna-seq-profiles/polyA-profile.csv

# a-aminitin profile
wget http://bioinfo/n/core/Bioinformatics/secondary/Bazzini/arb/MOLNG-2539/secundo/RSEM_TPM_table.csv
mv RSEM_TPM_table.csv rna-seq-profile/alpha-amanitin-prolife.csv

# spike-ins counts
wget http://bioinfo/n/core/Bioinformatics/secondary/Bazzini/arb/MOLNG-2539/secundo/ercc_count.csv
mv ercc_count.csv alpha-amanitin-prolife.csv
