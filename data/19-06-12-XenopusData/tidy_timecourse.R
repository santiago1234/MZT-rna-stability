library(tidyverse)
library(biomaRt)

# first deal with the ids -------------------------------------------------
# the ids reporter in the quant files are probabily from JGI
# i will take the gene name from the ones that have and recover the ensembl gene
# id

ensembl <- useMart("ensembl",dataset="xtropicalis_gene_ensembl")

ids_maps <- getBM(attributes=c('external_gene_name', 'ensembl_gene_id'), 
             mart = ensembl) %>% 
  as_tibble() %>% 
  filter(str_length(external_gene_name) != 0)

## first generate the ids mappings

ids_quants <- read_tsv("GSE65785_clutchApolyA_absolute_TPE_gene_isoform.txt", skip = 1) %>% 
  rename(
    gene_or_isoform = X1,
    external_gene_name = X2
  ) %>% 
  filter(gene_or_isoform == "Gene")
  
ids_maps <- 
  ids_quants %>% 
  dplyr::select(external_gene_name) %>% 
  separate(external_gene_name, into = c("jgi_id", "external_gene_name"), sep = "\\|") %>% 
  filter(str_length(external_gene_name) != 0) %>% 
  inner_join(ids_maps)

  

# tidy data sets ----------------------------------------------------------


laod_quants <- function(dpath, sample_condition) {
  # only quantifications at the gene level are maintained
  quants <- read_tsv(dpath, skip = 1) %>% 
    rename(
      gene_or_isoform = X1,
      external_gene_name = X2
    ) %>% 
    filter(gene_or_isoform == "Gene") # I dont neet the isophorms
  
  # add the enmble gene ids
  quants <- quants %>% 
    separate(external_gene_name, into = c("jgi_id", "external_gene_name"), sep = "\\|") %>% 
    filter(str_length(external_gene_name) != 0) %>% 
    dplyr::select(-gene_or_isoform) %>% 
    inner_join(ids_maps)
  
  # tidy the data
  quants %>% 
    gather(
      key = time, value = expression_quantification,
      -jgi_id, -external_gene_name, -ensembl_gene_id
    ) %>% 
    mutate(
      time = as.numeric(time),
      sample_condition = sample_condition
    )
  
}


timecourse_tidy <- 
  bind_rows(
  laod_quants("GSE65785_clutchA_rdRNA_absolute_TPE_gene_isoform.txt", sample_condition = "wt_ribo"),
  laod_quants("GSE65785_clutchApolyA_absolute_TPE_gene_isoform.txt", sample_condition = "wt_polya")
)

# define the top 3500 genes expressed at time 0 hrs
# I take the averge of the two cases



write_csv(timecourse_tidy, "time_course_xenopus.csv")


# - -----------------------------------------------------------------------
# - -----------------------------------------------------------------------
# get cds and 3’ UTR data -------------------------------------------------
# I donwloaded the data from the biomart website

utr <- read_tsv("xenopus_3utr.txt", col_names = c("id", "3utr")) %>% 
  separate(id, into = c("ensembl_gene_id", "transcript_id"), sep = "\\|") %>% 
  dplyr::select(-transcript_id)

coding <- read_tsv("xenopus_coding.txt", col_names = c("id", "coding")) %>% 
  separate(id, into = c("ensembl_gene_id", "transcript_id"), sep = "\\|") %>% 
  dplyr::select(-transcript_id)

coding <- coding %>% 
  filter(
    !str_detect(coding, "Sequence unavailable"),
    str_length(coding) %% 3 == 0, # coding should be a multiple of 3
    !str_detect(coding, "N") # No Ns is better
  ) %>% 
  group_by(ensembl_gene_id) %>% 
  do(
    arrange(., -str_length(coding)) %>% # get only 1 isophorm (keep longest)
      slice(1)
  ) %>% 
  ungroup()

utr <- 
  utr %>% 
  filter(
    !str_detect(`3utr`, "Sequenceunavailable"),
    !str_detect(`3utr`, "N")
  ) %>% 
  group_by(ensembl_gene_id) %>% 
  do(
    arrange(., -str_length(`3utr`)) %>% # get only 1 isophorm (keep longest)
      slice(1)
  ) %>% 
  ungroup()


left_join(coding, utr) %>% 
  write_csv("xenopus_seq_data_cds_3utr.csv")
