The crosstalk between codon optimality and 3’ UTR cis-elements dictates
mRNA stability
================
Santiago Gerardo Medina-Muñoz, Gopal Kushawah, Luciana Andrea
Castellano, Michay Diez, Michelle Lynn DeVore, María José Blanco
Salazar, Ariel A Bazzini
5/26/2020

# Published data

A short description of each dataset accompanying this publication.

-----

## **Table 1**: zebrafish RNA-seq time course during MZT

RNA-seq gene level quantifications. Raw sequencing data have been
deposited in NCBI Gene Expression Omnibus
[GSE148391](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE148391).
This table contains 35,117 rows and 31 columns. Each row represents a
gene.

The column names contains all the relevant sample description:

  - **Gene\_ID -\>** zebrafish ensembl gene id
  - **Treatment\_x-y\_hrs-RNAseq\_z** Here each column represents a
    single RNA-seq experiment containing the gene expression levels
    (Transcripts Per Million). The variables x, y, and z are placed
    holders: `x`: Some embryos were treated with alpha-amanitin to
    inhibit zygotic transcription. This takes only two values:
    `aamanitin` which denotes alpha-amanitin treated embryos and `wt`
    which represents untreated embryos. `y`: time post fertilization
    (hours). `z`: RNA-seq method `ribo` for ribosomal-RNA depletion and
    `polya` for poly-A selection.

A few rows and columns of the table are shown below:

| Gene\_ID           | Treatment\_wt-0\_hrs-RNAseq\_ribo | Treatment\_wt-1\_hrs-RNAseq\_ribo |
| :----------------- | --------------------------------: | --------------------------------: |
| ENSDARG00000000002 |                              2.71 |                              2.20 |
| ENSDARG00000000001 |                             12.02 |                             12.47 |
| ENSDARG00000000019 |                             49.06 |                             39.72 |
| ENSDARG00000000018 |                            105.91 |                             87.21 |
| ENSDARG00000000069 |                            219.59 |                            163.14 |

Table 1

-----

## **Table 2**: mRNA stability during MZT

Messenger RNA degradation rates estimated from the alpha-amanitin time
course. For additional information see the paper’s methods section:
*“Estimation of mRNA stability”*.

95% confidence interval lower limit Columns description:

  - **Gene\_ID -\>** zebrafish ensembl gene id.
  - **decay\_rate -\>** estimated decay rate, negative values indicate
    unstable genes and positive values stable genes.
  - **std.error -\>** standard error of the estimate.
  - **conf.low -\>** lower limit 95% confidence interval.
  - **conf.high -\>** upper limit 95% confidence interval.

A few rows of the table are shown below:

| Gene\_ID           | decay\_rate | std.error |    conf.low |   conf.high |
| :----------------- | ----------: | --------: | ----------: | ----------: |
| ENSDARG00000000001 | \-0.4317706 | 0.0610166 | \-0.5677241 | \-0.2958171 |
| ENSDARG00000000018 |   0.1329993 | 0.0487225 |   0.0244389 |   0.2415597 |
| ENSDARG00000000019 |   0.0667906 | 0.0190031 |   0.0244490 |   0.1091321 |
| ENSDARG00000000068 |   0.0057338 | 0.0202691 | \-0.0394285 |   0.0508961 |
| ENSDARG00000000069 | \-0.3366622 | 0.0318811 | \-0.4076978 | \-0.2656266 |
| ENSDARG00000000086 |   0.3106115 | 0.0428094 |   0.2152261 |   0.4059968 |

Table 2

-----

## **Table 3**: Training/testing data to train machine learning predictor of mRNA stability

Data that was used to train and evaluate machine learning preditor
model. This table contains 75,351 rows and 8 columns.

Columns description:

  - **gene\_id -\>** ensembl gene id.
  - **specie -\>** specie vertebrate (human = *H. sapiens*, fish = *D.
    rerio*, mouse = *M. musculus*, and xenopus = *X. tropicalis*).
  - **cell\_type -\>** cell type from where mRNA stability measurements
    were derived.
  - **datatype -\>** How where mRNA stability mesurements generated?:
      - `endogenous` Actinomycin D was used to block transcription
      - `aamanitin ribo` embryos treated with alpha-amanitin (RNA-seq
        ribosomal-RNA depletion)
      - `aamanitin polya` embryos treated with alpha-amanitin (RNA-seq
        poly-A selection)
      - `slam-seq` SLAM-seq.
  - **decay\_rate -\>** mRNA degradation rates (see the note below).
  - **utrlenlog -\>** 3’ UTR length log-transformed.
  - **cdslenlog -\>** cds length log-transformed.
  - **allocation -\>** this variable denotes wether the given
    observations was used for training or for testing (validation).
    Note: the same gene-id is never in training or testing
    simultaneously.

Notes:

  - The mRNA degradation rates in this table were standardized (mean = 0
    and standard deviation = 1). The next table below shows the mean and
    standard deviations of the original data. By applying the inverse
    transform, the values in the original scale can be recovered.
  - This table is missing the codon composition (codon frequencies).
    Codons frequencies were computed, for each gene, from the longest
    coding isoform sequence.

| specie  | cell\_type | datatype        | mean\_decayrate | stdeviation\_decayrate |
| :------ | :--------- | :-------------- | --------------: | ---------------------: |
| fish    | embryo mzt | aamanitin polya |     \-0.0629208 |              0.1860839 |
| fish    | embryo mzt | aamanitin ribo  |     \-0.0072156 |              0.0032341 |
| human   | 293t       | endogenous      |       0.0051652 |              0.0432408 |
| human   | hela       | endogenous      |     \-0.0009781 |              0.0542492 |
| human   | k562       | endogenous      |     \-0.0176985 |              0.0668504 |
| human   | k562       | slam-seq        |     \-0.1313344 |              0.0672748 |
| human   | RPE        | endogenous      |     \-0.0071829 |              0.0597563 |
| mouse   | mES cells  | slam-seq        |     \-0.1961037 |              0.0853733 |
| xenopus | embryo mzt | aamanitin ribo  |     \-0.0016150 |              0.0007771 |

means and strandar deviations

Decay rates, for the given specie, were obtained from the following
publications:

  - zebrafish: This study; AA Bazzini, F del Viso, MA Moreno‐Mateos… -
    The EMBO journal, 2016.
  - human: Q Wu, SG Medina, G Kushawah, ML DeVore… - Elife, 2019.
  - xenopus: AA Bazzini, F del Viso, MA Moreno‐Mateos… - The EMBO
    journal, 2016.
  - mouse: VA Herzog, B Reichholf, T Neumann, P Rescheneder… - Nature …,
    2017.

A few rows of the table are shown below excluding the `coding` column:

| gene\_id           | specie | cell\_type | datatype        | decay\_rate | utrlenlog | cdslenlog | allocation |
| :----------------- | :----- | :--------- | :-------------- | ----------: | --------: | --------: | :--------- |
| ENSG00000149428    | human  | 293t       | endogenous      |   0.9591786 |  7.705713 |  8.007700 | training   |
| ENSG00000183826    | human  | 293t       | endogenous      | \-0.2203851 |  8.784928 |  7.517521 | training   |
| ENSDARG00000007181 | fish   | embryo mzt | aamanitin ribo  |   0.6516244 |  5.236442 |  7.165494 | training   |
| ENSDARG00000025820 | fish   | embryo mzt | aamanitin ribo  | \-2.1045739 |  5.056246 |  7.321850 | training   |
| ENSDARG00000017244 | fish   | embryo mzt | aamanitin polya |   1.5186344 |  8.373323 |  9.088625 | training   |

Table 3

-----

## **Table 4**: Predictions and residual values during MZT for zebrafish and xenopus.

This data is part of **Fig. 2** (see paper). This table contains the
codon predicted stability for zebrafish and xenopus during MZT. This
table contains 10,668 rows and 4 columns.

Columns description:

  - **gene\_id -\>**: ensembl gene id, the genes here are maternally
    deposited.
  - **specie -\>**: either zebrafish or xenopus.
  - **residual -\>**: residual values, values close to zero indicate
    that the model predicts well the stability and large values that the
    model is far from the observed value. The residual is the difference
    between observed and predicted mRNA stability.

A few rows of the table are shown below:

| gene\_id           | specie  |   predicted |    residual |
| :----------------- | :------ | ----------: | ----------: |
| ENSXETG00000017158 | xenopus |   0.0736592 |   0.4100507 |
| ENSDARG00000041113 | fish    | \-0.0636230 |   0.6425244 |
| ENSXETG00000019722 | xenopus | \-0.6945703 | \-0.0959954 |
| ENSXETG00000020444 | xenopus |   0.4641538 |   0.7438621 |
| ENSXETG00000013779 | xenopus |   0.1119481 | \-2.0316289 |

Table 4

-----

## **Table 5**: Gene level mesurements of codon optimality.

See the methods section “Measuring codon optimality at the gene level”
and Additional file Figure S3. This table contains numerical
measurements of codon optimality in some endogenous genes, we have
generated two such measurements PLS1 and PLS2, positive values are
associated with enrichment in optimal codons (stabilizing codon) and
negative values with enrichment in non-optimal codons (destabilizing
codon).

This table contains 57,627 rows and 4 columns.

Columns description:

  - **gene\_id -\>**: ensembl gene id.
  - **PLS1 -\>** measurement 1 of codon optimality
  - **PLS2 -\>** measurement 2 of codon optimality
  - **specie -\>** vertebrate (human, mouse, xenopus, or zebrafish).

| gene\_id           |        PLS1 |        PLS2 | specie |
| :----------------- | ----------: | ----------: | :----- |
| ENSDARG00000099464 |   0.7106535 | \-2.4052864 | fish   |
| ENSG00000181396    |   2.6401594 | \-1.4055193 | human  |
| ENSG00000123576    |   0.6598256 | \-1.5668842 | human  |
| ENSMUSG00000040034 | \-3.7472500 |   2.0859899 | mouse  |
| ENSDARG00000101190 | \-1.1129394 |   0.7257975 | fish   |

Table 5

-----

## **Table 6**: Reporter sequences

This table contains 4 rows and 2 columns (Fig. 4e).

Columns description:

  - **sequence\_id -\>**: The sequence id.
  - **PLS1 -\>** The reporter DNA sequence.

Next, you can find the first column of this table.

| sequence\_id        |
| :------------------ |
| CODING-optimal      |
| CODING-non\_optimal |
| 3UTR-mir17          |
| 3UTR-mir17\_mutant  |

Table 6

-----

## **Table 7**: Codon frequencies of endogenous genes used to train machine learning model

This table contains the codon frequencies of the endogenous genes for
zebrafish, *Xenopus*, mouse, and human. The frequencies were determined
from the longest coding sequence isoform for each transcript.

Together this table and **Table 3** can be used to train the machine
learning model to predict mRNA stability.

Columns description:

  - **gene\_id -\>**: ensembl gene id.
  - **AAA -\>** frequency in transcript for codon AAA.
  - **AAC -\>** frequencies in transcript for codon AAC.
  - etc.

| gene\_id           | AAA | AAC | AAG | AAT |
| :----------------- | --: | --: | --: | --: |
| ENSDARG00000079869 |   5 |   7 |  12 |   3 |
| ENSG00000179604    |   1 |   0 |  17 |   5 |
| ENSG00000169813    |   6 |   8 |  10 |   3 |
| ENSG00000167700    |   1 |   1 |   7 |   2 |
| ENSG00000083454    |   7 |  16 |  22 |   4 |

Table 7
