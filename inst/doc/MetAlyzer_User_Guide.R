## ----setup, echo=FALSE, message=FALSE-----------------------------------------
knitr::opts_chunk$set(echo=TRUE, cache=TRUE, collapse=T, comment='#>')
library(MetAlyzer)
library(SummarizedExperiment)
library(ggplot2)
library(dplyr)

## ----install_cran, eval=FALSE-------------------------------------------------
#  install.packages("MetAlyzer")

## ----install_github, eval=FALSE-----------------------------------------------
#  library(devtools)
#  install_github("nilsmechtel/MetAlyzer")

## ----initialize_extraction----------------------------------------------------
fpath <- example_extraction_data()
metalyzer_se <- MetAlyzer_dataset(file_path = fpath)

metalyzer_se

## ----get_meta_data------------------------------------------------------------
meta_data <- colData(metalyzer_se)
head(meta_data)

## ----get_metabolites----------------------------------------------------------
metabolites <- rowData(metalyzer_se)
head(metabolites)

## ----get_concentration_values-------------------------------------------------
concentration_values <- assays(metalyzer_se)$conc_values
head(concentration_values, c(5, 5))

## ----get_quantification_status------------------------------------------------
quantification_status <- assays(metalyzer_se)$quant_status
head(quantification_status, c(5, 5))

## ----get_aggregated_data------------------------------------------------------
aggregated_data <- aggregatedData(metalyzer_se)
head(aggregated_data)

## ----filter_metabolites_extraction--------------------------------------------
metalyzer_se <- filterMetabolites(metalyzer_se, drop_metabolites = "Metabolism Indicators")
metalyzer_se

## ----filter_meta_data---------------------------------------------------------
metalyzer_se <- filterMetaData(metalyzer_se, `Sample Description` %in% 1:6)

## ----renameMetaData-----------------------------------------------------------
metalyzer_se <- renameMetaData(metalyzer_se, "Extraction_Method" = "Sample Description")

meta_data <- colData(metalyzer_se)
head(meta_data)

## ----load_replicates----------------------------------------------------------
replicate_meta_data <- example_meta_data()
head(replicate_meta_data)

## ----updateMetaData-----------------------------------------------------------
metalyzer_se <- updateMetaData(
  metalyzer_se,
  Date = Sys.Date(),
  Replicate = replicate_meta_data$Replicate
)

meta_data <- colData(metalyzer_se)
head(meta_data)

## ----calculateCV--------------------------------------------------------------
metalyzer_se <- calculate_cv(
  metalyzer_se,
  groups = c("Tissue", "Extraction_Method", "Metabolite"),
  cv_thresholds = c(0.1, 0.2, 0.3),
  na.rm = TRUE
)

aggregated_data <- aggregatedData(metalyzer_se) %>%
  select(c(Extraction_Method, Metabolite, Mean, SD, CV, CV_thresh))
head(aggregated_data)

## ----calculateANOVA-----------------------------------------------------------

metalyzer_se <- calculate_anova(
  metalyzer_se,
  categorical = "Extraction_Method",
  groups = c("Tissue", "Metabolite"),
  impute_perc_of_min = 0.2,
  impute_NA = TRUE
)

aggregated_data <- aggregatedData(metalyzer_se) %>%
  select(c(Extraction_Method, Metabolite, imputed_Conc, log2_Conc, ANOVA_n, ANOVA_Group))
head(aggregated_data)

## ----imputation_results-------------------------------------------------------
cat("Number of zero values before imputation:",
    sum(aggregatedData(metalyzer_se)$Concentration == 0, na.rm = TRUE), "\n")

cat("Number of zero values after imputation:",
    sum(aggregatedData(metalyzer_se)$imputed_Conc == 0, na.rm = TRUE), "\n")

## ----initialize_treatment-----------------------------------------------------
fpath <- example_mutation_data_xl()
metalyzer_se <- MetAlyzer_dataset(file_path = fpath)

metalyzer_se

## ----prepare_metabolites_treatment--------------------------------------------
metalyzer_se <- filterMetabolites(metalyzer_se, drop_metabolites = "Metabolism Indicators")
metalyzer_se

## ----show_sample_description--------------------------------------------------
meta_data <- colData(metalyzer_se)
meta_data$`Sample Description`

## ----prepare_control_mutant---------------------------------------------------
control_mutant <- factor(colData(metalyzer_se)$`Sample Description`, levels = c("Control", "Mutant"))
metalyzer_se <- updateMetaData(metalyzer_se, Control_Mutant = control_mutant)

meta_data <- colData(metalyzer_se)
meta_data$Control_Mutant

## ----calculate_log2FC---------------------------------------------------------
metalyzer_se <- calculate_log2FC(
  metalyzer_se,
  categorical = "Control_Mutant",
  impute_perc_of_min = 0.2,
  impute_NA = TRUE
)

## ----get_log2FC---------------------------------------------------------------
log2FC(metalyzer_se)

## ----plot_log2FC_vulcano, fig.width=7, fig.height=4.5-------------------------

log2fc_vulcano <- plot_log2FC(
  metalyzer_se,
  hide_labels_for = rownames(rowData(metalyzer_se)),
  vulcano=TRUE
)

log2fc_vulcano


## ----plot_log2FC_scatter, fig.width=9, fig.height=9---------------------------

log2fc_by_class <- plot_log2FC(
  metalyzer_se,
  hide_labels_for = rownames(rowData(metalyzer_se)),
  vulcano=FALSE
)

log2fc_by_class

## ----plot_network, fig.width=9, fig.height=9----------------------------------

log2fc_network <- plot_network(
  metalyzer_se,
  q_value=0.05,
  metabolite_text_size=2,
  connection_width=0.75,
  pathway_text_size=4,
  pathway_width=4,
  scale_colors = c("green", "black", "magenta")
)

log2fc_network

