library(tidyverse)

###########################################################################
# get the optimality ------------------------------------------------------
###########################################################################

optimality <-
  read_csv("../19-08-19-OverlapFinal/results_data/pls_species.csv") %>%
  filter(specie %in% c("fish", "xenopus"))

# get the orf and the 3’ UTR ----------------------------------------------

fish_orfs <- read_csv("../../data/19-01-17-Get-ORFS-UTRS-codon-composition/sequence-data/fish_seq_data_cds_3utr.csv")
xen_orfs <- read_csv("../../data/19-06-12-XenopusData/xenopus_seq_data_cds_3utr.csv")

dtafunctional <- bind_rows(fish_orfs, xen_orfs) %>%
  rename(gene_id = ensembl_gene_id) %>%
  right_join(optimality)

## add the % of optimal codons

opt_ratio <- function(orf) {
  # get the percentage of optimal codons according to bazzini embo2016
  codonr::count_codons(orf) %>%
    full_join(codonr::optimality_code_embo2016) %>%
    replace_na(list(n = 0)) %>%
    group_by(optimality) %>%
    summarise(suma = sum(n)) %>%
    spread(key = optimality, value = suma) %>%
    mutate(
      ratio_log = log2(optimal / (`non-optimal` + 1))
    ) %>%
    pull(ratio_log)
}

dtafunctional <-
  dtafunctional %>%
  mutate(
    optimality_ratio = map_dbl(coding, opt_ratio),
    miR430 = str_count(`3utr`, "GCACTT")
  )

write_csv(dtafunctional, "results_data/optimality_information.csv")


###########################################################################
# log 2 fold change -------------------------------------------------------
###########################################################################

xen_rna <- read_csv("../../data/19-06-12-XenopusData/time_course_xenopus.csv")

xen_rna <-
  xen_rna %>%
  rename(Gene_ID = ensembl_gene_id) %>%
  filter(
    time %in% c(2, 9),
    sample_condition == "wt_polya"
  ) %>%
  spread(key = time, value = expression_quantification) %>%
  arrange(-`2`) %>%
  slice(1:6000) %>%
  mutate(
    log2FC = log2(`9` / (`2` + 0.00001)),
    specie = "xenopus"
  ) %>%
  select(Gene_ID, log2FC, specie)


fish_fc <- read_csv("../../data/19-02-05-FoldChangeData/data/log2FC_earlyVSlate_tidytimecourse.csv") %>%
  filter(
    sample_condition == "wt_polya",
    time == 6,
    is_maternal
  ) %>%
  select(Gene_ID, log2FC) %>%
  mutate(specie = "fish")

foldchanges <-
  bind_rows(xen_rna, fish_fc) %>%
  filter(Gene_ID %in% dtafunctional$gene_id) %>%
  rename(gene_id = Gene_ID)

write_csv(foldchanges, "results_data/fold_change_late_vs_early.csv")

###########################################################################
# orthologs ---------------------------------------------------------------
###########################################################################

orto <- read_csv("../19-07-30-OrthologsEvoAnalysis/results_data/fish_xenopus_ortologs.csv") %>%
  select(`Gene stable ID`, `Xenopus gene stable ID`) %>%
  distinct() %>%
  rename(Gene_ID = `Gene stable ID`, xen_ort = `Xenopus gene stable ID`)

## load the maternal genes only

fish_maternal <- read_csv("../19-02-24-OverlapPathwaysFig3/results_data/regulatory_pathways_matrix.csv") %>%
  filter(is_maternal) %>%
  select(Gene_ID) %>%
  distinct()

orto <- inner_join(fish_maternal, orto)

# keep just 1 to 1 mapping

to_keep <-
  orto %>%
  group_by(Gene_ID) %>%
  summarise(n = n()) %>%
  filter(n == 1) %>%
  pull(Gene_ID)

orto <- filter(orto, Gene_ID %in% to_keep) %>%
  rename(ortolog_fish = Gene_ID, ortolog_xen = xen_ort)

write_csv(orto, "results_data/orthologs_fishMaternal_to_Xen.csv")
