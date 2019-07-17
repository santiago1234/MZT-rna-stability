## What is the strongest predictor of mRNA stability?

Here, I will compare pathways of mRNA stability.


### m6A pathways

I obtained the m6a target genes for mouse and human from [Batista et al](https://www.sciencedirect.com/science/article/pii/S1934590914004512#app3)

Data sets were obtained from the supplementary information of the papper:


**1-s2.0-S1934590914004512-mmc2.xlsx |** Table S1. All Mouse High-Confidence Peaks. Related to Figure 1 and Figure 4. Coordinates of m6A peaks in mouse genome (mm9), position of the m6A peak in the transcript, type of transcript, and gene symbol are displayed. For the Difference in Mettl3, the ratio between the IP and the Input is represented.

**1-s2.0-S1934590914004512-mmc4.xlsx |** Table S3. All Human High-Confidence Peaks. Related to Figure 5. Coordinates of m6A peaks in human genome (mm9), type of transcript, and gene symbols are shown.


### microRNAs

**TargetScan7.1__miR-291-3p_294-3p_295-3p_302-3p.predicted_targets.txt |** mouse miR-291a predicted targets. (these targets were used in the slam-seq paper).


### additional data

** |** mapping of mouse gene names to ensemble gene ids.
** mart_export.txt|** mappings human gene ids to transcript ids

### generate pathway information and response data (mRNA stability)


** make_predictor_data.R|** Generates the table X, this table contains the mRNA stability in the 
slam seq data for human and mouse and for fish it contains the log2FC at 6 hrs for the poly-a data
with respect to the initial time point (3 hrs). The table has the PLS decomposition of the codon space, the m6a
targets informaion and the microRNAs targets.
