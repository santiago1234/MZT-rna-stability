The crosstalk between codon optimality and 3' UTR cis-elements dictates mRNA stability
================
Santiago Gerardo Medina, Gopal Kushawah, María José Blanco Salazar, Ariel A Bazzini
5/26/2020

Published data
==============

A short description of the generated data accompanying this publication.

------------------------------------------------------------------------

**Table 1**: zebrafish RNA-seq time course during MZT
-----------------------------------------------------

RNA-seq gene level quantifications. Raw sequencing data have been deposited in NCBI Gene Expression Omnibus [GSE148391](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE148391). This table contains 35,117 rows and 31 columns. Each row represents a gene.

Columns description:

-   **Gene\_ID -&gt;** zebrafish ensembl gene id
-   **Treatment\_x|y\_hrs|RNAseq\_z** Here each column represents a single RNA-seq experiment containing the gene expression levels (Transcripts Per Million). The variables x, y, and z are placed holders: `x`: Some embryos were treated with [alpha-amanitin](https://en.wikipedia.org/wiki/Alpha-Amanitin) to inhibit zygotic transcription. This takes only two values: `aamanitin` which denotes alpha-amanitin treated embryos and `wt` which represents untreated embryos. `y`: time post fertilization (hours). `z`: RNA-seq method `ribo` for ribosomal-RNA depletion and `polya` for poly-A selection.

A few rows and columns of the table are shown below:

| Gene\_ID           |  Treatment\_wt|0\_hrs|RNAseq\_ribo|  Treatment\_wt|1\_hrs|RNAseq\_ribo|
|:-------------------|----------------------------------:|----------------------------------:|
| ENSDARG00000000001 |                              12.02|                              12.47|
| ENSDARG00000000068 |                              52.57|                              48.97|
| ENSDARG00000000069 |                             219.59|                             163.14|
| ENSDARG00000000019 |                              49.06|                              39.72|
| ENSDARG00000000002 |                               2.71|                               2.20|

------------------------------------------------------------------------

**Table 2**: mRNA stability during MZT
--------------------------------------

Messenger RNA degradation rates estimated from the alpha-amanitin time course. For additional information see the paper's methods section: *"Estimation of mRNA stability"*.

95% confidence interval lower limit Columns description:

-   **Gene\_ID -&gt;** zebrafish ensembl gene id.
-   **decay\_rate -&gt;** estimated decay rate, negative values indicate unstable genes and positive values stable genes.
-   **std.error -&gt;** standard error of the estimate.
-   **conf.low -&gt;** lower limit 95% confidence interval.
-   **conf.high -&gt;** upper limit 95% confidence interval.

A few rows of the table are shown below:

| Gene\_ID           |  decay\_rate|  std.error|    conf.low|   conf.high|
|:-------------------|------------:|----------:|-----------:|-----------:|
| ENSDARG00000000001 |   -0.4317706|  0.0610166|  -0.5677241|  -0.2958171|
| ENSDARG00000000018 |    0.1329993|  0.0487225|   0.0244389|   0.2415597|
| ENSDARG00000000019 |    0.0667906|  0.0190031|   0.0244490|   0.1091321|
| ENSDARG00000000068 |    0.0057338|  0.0202691|  -0.0394285|   0.0508961|
| ENSDARG00000000069 |   -0.3366622|  0.0318811|  -0.4076978|  -0.2656266|
| ENSDARG00000000086 |    0.3106115|  0.0428094|   0.2152261|   0.4059968|

------------------------------------------------------------------------

**Table 3**: Training/testing data to train machine learning predictor of mRNA stability
----------------------------------------------------------------------------------------

Data that was used to train and evaluate machine learning preditor model. This table contains 75,351 rows and 8 columns.

Columns description:

-   **gene\_id -&gt;** ensembl gene id.
-   **specie -&gt;** specie vertebrate (human = *H. sapiens*, fish = *D. rerio*, mouse = *M. musculus*, and xenopus = *X. tropicalis*).
-   **cell\_type -&gt;** cell type from where mRNA stability measurements were derived.
-   **datatype -&gt;** How where mRNA stability mesurements generated?:
    -   `endogenous` Actinomycin D was used to block transcription
    -   `aamanitin ribo` embryos treated with alpha-amanitin (RNA-seq ribosomal-RNA depletion)
    -   `aamanitin polya` embryos treated with alpha-amanitin (RNA-seq poly-A selection)
    -   `slam-seq` SLAM-seq.
-   **decay\_rate -&gt;** mRNA degradation rates (see the note below).
-   **utrlenlog -&gt;** 3' UTR length log-transformed.
-   **cdslenlog -&gt;** cds length log-transformed.
-   **allocation -&gt;** this variable denotes wether the given observations was used for training or for testing (validation). Note: the same gene-id is never in training or testing simultaneously.

Notes:

-   The mRNA degradation rates in this table were standardized (mean = 0 and standard deviation = 1). The next table below shows the mean and standard deviations of the original data. By applying the inverse transform, the values in the original scale can be recovered.
-   This table is missing the codon composition (codon frequencies). Codons frequencies were computed, for each gene, from the longest coding isoform sequence.

| specie  | cell\_type | datatype        |  mean\_decayrate|  stdeviation\_decayrate|
|:--------|:-----------|:----------------|----------------:|-----------------------:|
| fish    | embryo mzt | aamanitin polya |       -0.0629208|               0.1860839|
| fish    | embryo mzt | aamanitin ribo  |       -0.0072156|               0.0032341|
| human   | 293t       | endogenous      |        0.0051652|               0.0432408|
| human   | hela       | endogenous      |       -0.0009781|               0.0542492|
| human   | k562       | endogenous      |       -0.0176985|               0.0668504|
| human   | k562       | slam-seq        |       -0.1313344|               0.0672748|
| human   | RPE        | endogenous      |       -0.0071829|               0.0597563|
| mouse   | mES cells  | slam-seq        |       -0.1961037|               0.0853733|
| xenopus | embryo mzt | aamanitin ribo  |       -0.0016150|               0.0007771|

Decay rates, for the given specie, were obtained from the following publications:

-   zebrafish: This study; AA Bazzini, F del Viso, MA Moreno‐Mateos… - The EMBO journal, 2016.
-   human: Q Wu, SG Medina, G Kushawah, ML DeVore… - Elife, 2019.
-   xenopus: AA Bazzini, F del Viso, MA Moreno‐Mateos… - The EMBO journal, 2016.
-   mouse: VA Herzog, B Reichholf, T Neumann, P Rescheneder… - Nature …, 2017.

A few rows of the table are shown below excluding the `coding` column:

| gene\_id           | specie | cell\_type | datatype   |  decay\_rate|  utrlenlog|  cdslenlog| allocation |
|:-------------------|:-------|:-----------|:-----------|------------:|----------:|----------:|:-----------|
| ENSG00000013523    | human  | RPE        | endogenous |    0.2978391|   8.023225|   7.607878| training   |
| ENSG00000131943    | human  | k562       | slam-seq   |   -0.8748520|   8.249837|   6.131227| training   |
| ENSG00000104325    | human  | 293t       | endogenous |   -0.1878186|   7.459915|   6.916715| training   |
| ENSG00000227124    | human  | k562       | endogenous |    0.0023865|   6.916715|   7.917901| training   |
| ENSMUSG00000018865 | mouse  | mES cells  | slam-seq   |    0.6431109|   7.225481|   6.752270| training   |

------------------------------------------------------------------------

**Table 4**: Predictions and residual values during MZT for zebrafish and xenopus.
----------------------------------------------------------------------------------

This data is part of **Fig. 2** (see paper). This table contains the codon predicted stability for zebrafish and xenopus during MZT. This table contains 10,668 rows and 4 columns.

Columns description:

-   **gene\_id -&gt;**: ensembl gene id, the genes here are maternally deposited.
-   **specie -&gt;**: either zebrafish or xenopus.
-   **residual -&gt;**: residual values, values close to zero indicate that the model predicts well the stability and large values that the model is far from the observed value. The residual is the difference between observed and predicted mRNA stability.

A few rows of the table are shown below:

| gene\_id           | specie  |   predicted|    residual|
|:-------------------|:--------|-----------:|-----------:|
| ENSDARG00000006031 | fish    |   0.2479914|  -0.0736010|
| ENSXETG00000013410 | xenopus |   0.3176923|  -1.2910357|
| ENSDARG00000070447 | fish    |  -0.0268868|  -0.0312370|
| ENSXETG00000024933 | xenopus |  -0.3469115|   1.4982524|
| ENSXETG00000020045 | xenopus |  -0.6369484|   0.6133757|

**Table 5**: Gene level mesurements of codon optimality.
--------------------------------------------------------

See the methods section “Measuring codon optimality at the gene level” and Additional file Figure S3. This table contains numerical measurements of codon optimality in some endogenous genes, we have generated two such measurements PLS1 and PLS2, positive values are associated with enrichment in optimal codons (stabilizing codon) and negative values with enrichment in non-optimal codons (destabilizing codon).

This table contains 57,627 rows and 4 columns.

Columns description:

-   **gene\_id -&gt;**: ensembl gene id.
-   **PLS1 -&gt;** measurement 1 of codon optimality
-   **PLS2 -&gt;** measurement 2 of codon optimality
-   **specie -&gt;** vertebrate (human, mouse, xenopus, or zebrafish).

| gene\_id           |        PLS1|        PLS2| specie  |
|:-------------------|-----------:|-----------:|:--------|
| ENSXETG00000027421 |   3.5829424|  -1.9113006| xenopus |
| ENSMUSG00000031731 |  -1.0684200|   2.9623095| mouse   |
| ENSMUSG00000022000 |  -2.3651121|   0.9490589| mouse   |
| ENSG00000197647    |  -6.3191219|  -5.4298493| human   |
| ENSMUSG00000030538 |   0.4288678|  -3.3887081| mouse   |

------------------------------------------------------------------------

Reporters sequences
-------------------

Bewlow you can find the coding and 3' UTR sequences that were used to clone the reporter genes used in **Fig. 4e** and **Figure S4e**.

#### Coding

Optimal sequence: *GACATCTTTGGCTTTGAGAACTTTGAGGTCAACCGCTTTGAGCAGTTCAACATTAACTATGCAAACGAGAAGCTTCAGGAGTATTTCAACAAGCACATTTTCTCACTGGAGCAGCTTGAGTTCAGGAAGGTGCAGCATGAGCTGGAGGAGGCTCAGGAGAGAGCTGACATCGCCGAGTCCCAGGTCAACAAGCTCAGAGCTAAAAGCCGTGAATTTGGAAAGGGTAAAGAGGCTGAGGAGGCTGACTCCTTCGACTATAAGAGCTTCTTCGCCAAGGTTGGGCTGTCCGCCAAGACTCCTGATGACATCAAGAAGGCTTTTGCTGTCATTGACCAGGACAAGAGCGGCTTCATTGAGGAGGATGTGGAGGACTCCCTCTGTGAGGCCAAAGAGCTGTTCATCAAGACAGTCAAGCACTTCGGTGAGGACGCTGATAAGATGCAGCCTGATGAGTTCTTTGGGATTTTCGACCAGTTCTTGCGTATCCCCAAGGAGCAGGGCTTCCTGTCGTTCTGGAGAGGAAACTTGGCCAACGTCATCAGATACTTCCCCACACAGGCCCTCAACTTTGCTTTCAAGGACAAGTACAAGAAGGTCTTCGACATCACAGACAAGCTGGAGAACGAGCTGGCCAATAAGGAGGCTTTCCTCAGACAGATGGAGGAGAAGAACAGGCAGTTGCAGGAGCGGCTTGAGTTGGCAGAGCAGAAGCTCCAGCAG*

Non-optimal sequence: *ACATCTTTGGCTTTGAGAACTTTGAGGTCAACCGCTTTGAGCAGTTCAACATTAACTATGCAAACGAGAAGCTTCAGGAGTATTTCAACAAGCACATTTTCTCACTGGAGCAGCTTGAGTTCAGGAAGGTGCAGCATGAGCTGGAGGAGGCTCAGGAGAGAGCTGACATCGCCGAGTCCCAGGTCAACAAGCTCAGAGCTAAAAGCCGTGAATTTGGAAAGGGTAAAGAGGCTGAGGAGGCTGACTCCTTCGACTATAAGAGCTTCTTCGCCAAGGTTGGGCTGTCCGCCAAGACTCCTGATGACATCAAGAAGGCTTTTGCTGTCATTGACCAGGACAAGAGCGGCTTCATTGAGGAGGATGTGGAGGACTCCCTCTGTGAGGCCAAAGAGCTGTTCATCAAGACAGTCAAGCACTTCGGTGAGGACGCTGATAAGATGCAGCCTGATGAGTTCTTTGGGATTTTCGACCAGTTCTTGCGTATCCCCAAGGAGCAGGGCTTCCTGTCGTTCTGGAGAGGAAACTTGGCCAACGTCATCAGATACTTCCCCACACAGGCCCTCAACTTTGCTTTCAAGGACAAGTACAAGAAGGTCTTCGACATCACAGACAAGCTGGAGAACGAGCTGGCCAATAAGGAGGCTTTCCTCAGACAGATGGAGGAGAAGAACAGGCAGTTGCAGGAGCGGCTTGAGTTGGCAGAGCAGAAGCTCCAGCAGG*

#### 3' UTRs

3' UTR weak seeds. This UTR was used in **Figure S4e**. This UTR contains a 6-mer and a 7-mer:

*CTTCCACCAATAGAGGAACTAGGAAACAACCAATGGGCTCTGATAAATCAGAGGAAAGGTGAAGAGAAAAAAAACCTGTGCTGCCACAATACCCCGAGAGTGCATGGAATAGTTACACATTAAACAACTGTGAAGAAGAGAGACACAGTTAGAACAAATACAGTATCCAAAAC***GCACTT***CTTACATTC***AGCACTT***GCAGTTTTCGTTTTACTTCATTATTTTGAAATAATAACAGTAGACACAACTTCAGCCATATGCGCTTACTTTGCAcCCAATAATCCGTTTTCATAAAGAACGACTTACATTTTATGTTTTAATAACTTGTTAATTGTATTATGGGAACAGATTATATGGTCTGTTATTAATGGCTCTTACTAATAAATGTGTCAATTGAATTATTGGCTGTTTTT*

3' UTR strong seeds. This UTR was used in **Fig. 4e**. This UTR contains two 7-mers and an 8-mer:

*TGATCCAACGTGGAGTTCTCGCATCTTGAAATGTCCAAATATATTATCATTATTATTTTCCCTCTTCATTTGCTCTTCAACGAAAGGTGAACAGGTGGAACGATATTTATTTTTCTAGTCTGTAATTTTCTGCTGAAAGTATTTTTTTGTGTGTATGTTAAAACAAATGTGGAGAGGTTTCAAAGACCAAACGATTGGCTGACGAAGACGAATATCTTCCTAATAATCGATGAATGAAACTATAAATCAGGTGAAAATGAACTGGCGGAGTTTCGAGATGAGTTATG***AGCACTT***CTGAGGGTGTTACCAGCGTGAAATGAAGCCAGTGCAGCCCAGCTGTCTTAAGTCCCACACAAAAACT***GCACTTA***ACTG***AGCACTTA***AAAAGAAAAACATATTTGACCACAACAGTCGAAGCAGTTTCCAAAATTTAGGAAAGACTTTTAATTTATTTTTCAGTATTCTGCTAACTTAAAAGGCCATCATGTGTTGTGTGTATGTGACATGTGCCATTACCTCTGCTGTGGCACCTTTTCCTTCTCAGGTAAAATCTGCTTATATTTGTACTT*

Note: The mutated versions of these 3' UTRs were created by introducing nucleotides substitutions in the microRNA seed site. This should disrupt microRNA recognition.
