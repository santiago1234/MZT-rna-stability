## Generate log2FC data

I generated a tidy version of the time course: *data/RNAseq_tidytimecourse.csv*

I computed the log2 fold change of every time point with the earliest time point (i.e. 0/0, 1/0, ..., 8/0): *data/log2FC_earlyVSlate_tidytimecourse.csv*

The data contains a column `is_maternal`, where the genes with true are maternal, although i am not 100% sure. But I can take it as true.

