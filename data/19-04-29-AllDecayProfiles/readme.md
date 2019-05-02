## DataSets

+ **./fishxenopus2016/inline-supplementary-material-2.xls|** Suplementary file from ariel paper with the half-life for fish and xenopus (aamanitin-ribo).
+ **/mouseStemCellsSLAMseqNatureBiotechnology/nmeth.4435-S4.xls|** Supplementary Table 2
High-confidence half-life data for 6665 transcripts in mES cells.Table is attached as Excel file and reports for each transcript the chromosome, start and end of the counting windows identifying transcript 3′ ends, Name of the associated transcript, length of the counting window, half-life (h), decay rate (k; h-1), standard error of half-life and decay rate (error of single exponential fit), and accuracy of fit (rsquare).
+ **Wu et al_decay_rates.xlsx |** -> Decay rate profiles for human. From this file I only
extracted the data for RPE and HELA cells.


## Main OUtput

** decay_profiles_and_seqdata.csv |** all the decay profiles for fish,
human, xenopus and mouse with sequence data.

this table contains the following columns:

+ *gene_id |* ensemble gene id, for the orf-ome data the ORF-barcode Name.
+ *datatype |* the datatype used to estimate the decay rates. (aamanitin polya, aamanitin ribo, endogenous, orf-ome, slam-seq)
+ *cell_type |* cell type/tissue (293t, embryo mzt, k562, mES cells)
+ *specie |* specie (fish, xenopus, mouse, human)
+ *decay_rate |* decay rate estimate. NOTE: estimates may be in different scales (consider scaling by specie, datatype, cell_type) befor analysis
+ *coding |* coding dna sequence
+ *3utr |* 3' utr sequence
