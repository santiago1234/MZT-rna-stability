# make dir to put the files

mkdir -p rna-seq-profiles

# get ribo profile
wget http://bioinfo/n/core/Bioinformatics/secondary/Bazzini/arb/MOLNG-2541/secundo/RSEM_TPM_table.csv
mv RSEM_TPM_table.csv rna-seq-profiles/ribo-zero-profile.csv


# poly-A profile
wget http://bioinfo/n/core/Bioinformatics/secondary/Bazzini/arb/MOLNG-2540/secundo/RSEM_TPM_table.csv
mv RSEM_TPM_table.csv rna-seq-profiles/polyA-profile.csv
