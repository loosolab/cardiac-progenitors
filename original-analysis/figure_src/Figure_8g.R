library(magrittr)

#
# Load scATAC-seq data
#
load("data/scData_filtered_2.Rda")

#
# Load scRNA-seq data
#
RNAseq <- read.table("data/Isl1_endo_branch_expression.txt", header = T, row.names = 1, sep = "\t")

RNAseq <- RNAseq[c("Sox13", "Sox6", "Sox9", "Tead1", "Gata6", "Zeb1", "Fosl2", "Junb", "Gata2", "Tal1"), ]
RNAseq <- RNAseq[c("Gata4", "Gata5", "Sox6", "Sox9", "Sox13", "Tead1", "Tead2", "Tead4", "Gata6", "Gata2"), ]
split <- c(rep("first", 6), rep("second", 4))

#RNAseq_scaled <- t(scale(t(RNAseq), scale = F))
RNAseq_smoothed <- t(zoo::rollapply(t(RNAseq), width = 15, FUN = mean, fill = "extend", by.column = T))

#
# Get TF deviation data
#
load("data/dev_cluster_specific_2.Rda")

# Extract deviation scores
tfs_limited <- c("SOX13_MA1120.1", "Sox6_MA0515.1", "SOX9_MA0077.1", "TEAD1_MA0090.2",
                 "GATA6_MA1104.1",
                 "ZEB1_MA0103.3", "FOSL2_JUNB_MA1138.1", "FOSL2_JUNB_MA1138.1", "GATA2_MA0036.3", "GATA1_TAL1_MA0140.2")
tfs_limited <- c("Gata4_MA0482.1", "GATA5_MA0766.1", "Sox6_MA0515.1", "SOX9_MA0077.1", "SOX13_MA1120.1", "TEAD1_MA0090.2", "TEAD2_MA1121.1", "TEAD4_MA0809.1", "GATA6_MA1104.1", "GATA2_MA0036.3")

var_data_limited <- chromVAR::deviationScores(dev)[tfs_limited, ]

#
# Get cluster specific cells ordered by their pseudotime
#
heatmaps <- list("reference" = list(branch = "cardiac", cluster = 2),
                 "intermediate" = list(branch = "cl3", cluster = 3),
                 "endo" = list(branch = "endo", cluster = 5))

var_data_list <- lapply(heatmaps, function(x) {
  i <- which(SummarizedExperiment::colData(scData_filtered)$.cluster_5 == x$cluster)
  data <- var_data_limited[, i]
  if(!is.na(x$branch)) {
    data <- data[, order(SummarizedExperiment::colData(scData_filtered)[i, paste0("dpt_", x$branch)])]
  }
  data
})
var_data_limited <- do.call("cbind", var_data_list)

#
# Smooth variability data
#
var_data_limited_smoothed <- t(zoo::rollapply(t(var_data_limited), width = 13, FUN = mean, fill = "extend", by.column = T))

#
# Define color mapping
#
col_fun_rna <- circlize::colorRamp2(seq(from = -4, to = 4, length.out = 11), colors = rev(RColorBrewer::brewer.pal(11, "RdBu")))
col_fun_atac <- circlize::colorRamp2(seq(from = -1, to = 5, length.out = 9), colors = RColorBrewer::brewer.pal(9, "YlGnBu"))

#
# Draw Heatmap
#
h_left <- ComplexHeatmap::Heatmap(RNAseq_smoothed,
                                  col = col_fun_rna,
                                  split = split,
                                  cluster_columns = F,
                                  cluster_rows = F,
                                  show_row_names = T,
                                  row_names_side = "left",
                                  show_row_dend = F,
                                  name = "Expression",
                                  column_title_side = "top",
                                  column_title = "RNA expression of TF",
                                  width = grid::unit(90, "mm"),
                                  heatmap_legend_param = list(legend_direction = "horizontal", at = c(-4, 4), labels = c("low", "high")))

h_right <- ComplexHeatmap::Heatmap(var_data_limited_smoothed,
                                   col = col_fun_atac,
                                   split = split,
                                   cluster_columns = F,
                                   cluster_rows = F,
                                   show_column_names = F,
                                   show_row_names = F,
                                   show_row_dend = F,
                                   name = "Motif accessibility",
                                   column_title_side = "top",
                                   column_title = "Motif accessibility of TF",
                                   width = grid::unit(90, "mm"),
                                   heatmap_legend_param = list(legend_direction = "horizontal", at = c(-1, 5), labels = c("low", "high")))

plot_endo_RNA_ATAC <- grid::grid.grabExpr(ComplexHeatmap::draw(ComplexHeatmap::add_heatmap(h_left, h_right), column_title_gp = grid::gpar(fontface = "bold"), heatmap_legend_side = "bottom"))

figure_endo_RNA_ATAC <- multipanelfigure::multi_panel_figure(width = 230, height = 205, columns = 1, rows = 1)
figure_endo_RNA_ATAC <- multipanelfigure::fill_panel(figure_endo_RNA_ATAC, plot_endo_RNA_ATAC)

ggplot2::ggsave(figure_endo_RNA_ATAC,
                file = "figures_2/Endo_RNA_ATAC.pdf",
                height = 205,
                width = 230,
                unit = "mm")
