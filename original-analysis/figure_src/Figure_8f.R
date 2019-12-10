library(magrittr)

#
# Load scATAC-seq data
#
load("data/scData_filtered_2.Rda")

#
# Load scRNA-seq data
#
RNAseq <- read.table("data/Isl1_cardiac_branch_expression.txt", header = T, row.names = 1, sep = "\t")

#RNAseq <- RNAseq[c("Hoxa9", "Hoxc10", "Hoxa10", "Hoxd8", "Foxc1", "Foxo4", "Tcf4", "Foxp1", "Gata4", "Gata6", "Gata5"), ]
RNAseq <- RNAseq[c("Hoxa9", "Hoxc10", "Hoxd8", "Tead1", "Tead2", "Tead4", "Hand1", "Tbx5", "Gata4", "Gata6"), ]
split <- c(rep("first", 3), rep("second", 3), rep("third", 4))

#RNAseq_scaled <- t(scale(t(RNAseq), scale = F))
RNAseq_smoothed <- t(zoo::rollapply(t(RNAseq), width = 15, FUN = mean, fill = "extend", by.column = T))

#
# Get TF deviation data
#
load("data/dev_cluster_specific_2.Rda")

# Extract deviation scores
tfs_limited <- c("Hoxa9_MA0594.1", "HOXC10_MA0905.1", "HOXA10_MA0899.1", "Hoxd8_MA0910.1",
                 "FOXC1_MA0032.2", "FOXO4_MA0848.1", "TCF4_MA0830.1", "FOXP1_MA0481.2",
                 "Gata4_MA0482.1", "GATA6_MA1104.1", "GATA5_MA0766.1")
tfs_limited <- c("Hoxa9_MA0594.1", "HOXC10_MA0905.1", "Hoxd8_MA0910.1",
                "TEAD1_MA0090.2", "TEAD2_MA1121.1", "TEAD4_MA0809.1",
                "Hand1_Tcf3_MA0092.1", "TBX5_MA0807.1","Gata4_MA0482.1", "GATA6_MA1104.1")

var_data_limited <- chromVAR::deviationScores(dev)[tfs_limited, ]

#
# Get cluster specific cells with their pseudotime
#
cells <- which(SummarizedExperiment::colData(scData_filtered)$.cluster_5 %in% c(1, 2))

var_data_limited <- var_data_limited[, cells]

# Order var_data by dpt
var_data_limited <- var_data_limited[, order(SummarizedExperiment::colData(scData_filtered)$dpt_cardiac[cells])]

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
                                  cluster_rows = T,
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

plot_cardiac_RNA_ATAC <- grid::grid.grabExpr(ComplexHeatmap::draw(ComplexHeatmap::add_heatmap(h_left, h_right), column_title_gp = grid::gpar(fontface = "bold"), heatmap_legend_side = "bottom"))

figure_cardiac_RNA_ATAC <- multipanelfigure::multi_panel_figure(width = 230, height = 205, columns = 1, rows = 1)
figure_cardiac_RNA_ATAC <- multipanelfigure::fill_panel(figure_cardiac_RNA_ATAC, plot_cardiac_RNA_ATAC)

ggplot2::ggsave(figure_cardiac_RNA_ATAC,
       file = "figures_2/Cardiac_RNA_ATAC.pdf",
       height = 205,
       width = 230,
       unit = "mm")
