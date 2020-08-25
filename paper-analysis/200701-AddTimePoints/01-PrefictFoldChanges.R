library(tidyverse)

# input files -------------------------------------------------------------

predicted_stability_infile <- snakemake@input$model_predicted_stability
precomputed_fish_foldchanges_infile <- snakemake@input$foldchange_data_fish
xen_timecourse_infile <- snakemake@input$expression_profile_xenopus
table_with_genesList_infile <- snakemake@input$gene_list

outfile <- snakemake@output[[1]]
# predicted stability -----------------------------------------------------

predicted_stability <- read_csv(predicted_stability_infile) %>% 
  select(-coding, -cell_type, -datatype, -utrlenlog, -cdslenlog)

## helper function

get_results <- function(df_grp) {
  df_grp
  # fit the model to take the predicted_stability into
  # a better scale
  fit <- lm(log2FC ~ predicted_stability, data = df_grp)
  df_grp$pred_log2FC <- predict(fit)
  df_grp$residual <- resid(fit)
  df_grp
}

# fish --------------------------------------------------------------------

stability_mzt_fish <- read_csv(precomputed_fish_foldchanges_infile) %>% 
  filter(
    sample_condition == "wt_ribo",
    time %in% 4:8,
    is_maternal
  ) %>% 
  filter(
    !is.infinite(log2FC),
    !is.na(log2FC)
  )

stability_mzt_fish <- 
  stability_mzt_fish %>% 
  select(Gene_ID, log2FC, time) %>% 
  rename(gene_id = Gene_ID)

## add the predicted stability

stability_mzt_fish <- inner_join(stability_mzt_fish, predicted_stability)

# fit a linear model to scale and add the predicted during MZT and --------

fish_predictions <- 
  stability_mzt_fish %>% 
  group_by(time) %>% 
  nest() %>% 
  mutate(
    resultado = map(data, get_results)
  ) %>% 
  select(time, resultado) %>% 
  unnest(resultado) %>% 
  ungroup()


# xenopus -----------------------------------------------------------------

xen_rna <- read_csv(xen_timecourse_infile)

# subset the genes that i've been using in the paper
# the top 6K genes expressed at 2hrs
xen_genes <- read_csv(table_with_genesList_infile) %>% 
  filter(specie == "xenopus") %>% 
  pull(gene_id)

# subset timepoints 2-12 hrs and wt_polya data
xen_rna <-
  xen_rna %>%
  rename(Gene_ID = ensembl_gene_id) %>%
  filter(
    time %in% c(2:20),
    sample_condition == "wt_polya",
    Gene_ID %in% xen_genes
  ) %>% 
  rename(gene_id = Gene_ID) %>% 
  select(gene_id, time, expression_quantification)


# compute the log2FC ------------------------------------------------------

xen_rna_expression_at_2 <- xen_rna %>% 
  filter(time == 2) %>% 
  rename(exp_at_2hrs = expression_quantification) %>% 
  select(-time)
  
stability_mzt_xen <- 
  xen_rna %>% 
  inner_join(xen_rna_expression_at_2) %>% 
  mutate(
    log2FC = log2(expression_quantification / (exp_at_2hrs)),
  ) %>% 
  select(gene_id, time, log2FC) %>% 
  filter(time > 2, !is.infinite(log2FC))

xen_predictions <- 
  stability_mzt_xen %>% 
  inner_join(predicted_stability) %>% 
  group_by(time) %>% 
  nest() %>% 
  mutate(
    resultado = map(data, get_results)
  ) %>% 
  select(time, resultado) %>% 
  unnest(resultado) %>% 
  ungroup()

# save results ------------------------------------------------------------

bind_rows(
  xen_predictions,
  fish_predictions
) %>% 
  write_csv(outfile)



