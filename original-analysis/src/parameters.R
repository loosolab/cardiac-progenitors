library(scales)
#
# General options
#
parameters <- list(
  "general" = list("project" = 'Cardiac progenitor single cells',
                   "path_rdata" = 'data',
                   "path_rextdata" = "ext_data",
                   "path_rfigures" = 'supplementary_figures',
                   "path_supdata" = 'supplementary_data',
                   "lower_detection_limit" = 10,
                   "biomart_dataset" = "mmusculus_gene_ensembl"),
  "filtering" = list("qc_columns" = c("total_features", "pct_dropout", "totalalignments", "exonic", "intronic", "intergenic", "secondaryalignments","exprs_feature_controls_housekeeping", "pct_counts_feature_controls_mito", "percent_genes"),
                     "log_qc_columns" = c("total_features", "totalalignments", "exonic", "intronic", "intergenic","secondaryalignments", "exprs_feature_controls_housekeeping"),
                     "cell_outlier_columns" = list(pct_counts_feature_controls_mito = c(1.5, "higher", F), total_features = c(2, "both", F), pct_dropout = c(2, "higher", F), exprs_feature_controls_housekeeping = c(2, "lower", F), percent_genes = c(1.5, "both", F)),
                     "cell_outlier_maxfailed" = 1,
                     "EpOverA_A" = 2000,
                     "kOverA_k" = 10,
                     "kOverA_A" = 10),
  "normalisation" = list("cycle_markers" = "mouse_cycle_markers.rds"),
  "heterogeneity" = list("het_zscore" = qnorm(0.99)),
  "diffusionmaps" = list("de_fdr" = 0.01,
                         "de_abs_logfc" = 2,
                         "marker_fdr" = 0.01,
                         "marker_auroc" = 0.8),
  "correlation" = list(),
  "colors" = list("nkx_cluster_palette" = c("#8B0000", "#CDAD00", "#00008B"),
                  "nkx_timepoint_palette" = c("#FF003C", "#248F8D", "#987F69"),
                  "nkx_run_palette" = c("#ECD078", "#C02942", "#53777A", "#542437"),
                  "nkxko_cluster_palette" = c(alpha(c("#00008B", "#8B0000", "#CDAD00"), 1), "red"),
                  "isl_cluster_palette" = c("#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "grey"),
                  "isl_timepoint_palette" = c("#FF003C", "#248F8D", "#987F69"),
                  "islko_cluster_palette" = c(alpha(c("#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E"), 0.8), "#FF0001", "grey"),
                  "run_palette" = c("#ECD078", "#C02942", "#53777A", "#542437", "#0B486B", "#D95B43"))
)

# #
# # Options for chapter 02_filtering.Rmd
# #
# # [cellfilter]
# qc_columns <- c("total_features", "pct_dropout", "totalalignments", "exonic", "intronic", "intergenic", "secondaryalignments","exprs_feature_controls_housekeeping", "pct_counts_feature_controls_mito", "percent_genes")
# log_qc_columns <- c("total_features", "totalalignments", "exonic", "intronic", "intergenic","secondaryalignments", "exprs_feature_controls_housekeeping")
# cell_outlier_columns <- list(pct_counts_feature_controls_mito = c(1.5, "higher", F), total_features = c(2, "both", F), pct_dropout = c(2, "higher", F), exprs_feature_controls_housekeeping = c(2, "lower", F), percent_genes = c(1.5, "both", F))
# cell_outlier_maxfailed <- 1
#
# # [genefilter]
# EpOverA_A <- 2000
# kOverA_k <- 10
# kOverA_A <- lower_detection_limit
#
# # [downsampling]
#
# #
# # Options for chapter 03_normalisation.Rmd
# #
# # [cell cylce]
# cycle_markers <- "mouse_cycle_markers.rds"
#
# # [PCA]
# init_pca_variables = c("Platform", "pct_dropout", "log10_total_counts", "pct_counts_feature_controls_mito", "pct_counts_top_200_endogenous_features", "totalalignments", "total_features", "cellcycle", "Type", "Chip")
#
# #
# # Options for chapter 04_heterogeneity.Rmd
# #
# het_zscore <- qnorm(0.99)
#
# #
# # Options for chapter 05_clustering.Rmd
# #
# timpoint_col_palette <- "Set1"
# cycle_col_palette <- "Set1"
# cluster_col_palette <- "Set2"
# cluster_na_col <- "#969696"
#
# # [differential expression]
# logfc_threshold <- 3
# fdr_threshold <- 0.01
#
# # [violin plots]
# show_max_de_genes <- 20
#
# #
# # Options for 07_diffusion.Rmd
# #
# proj_col <- "#bd0026"



