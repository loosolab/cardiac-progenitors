library(magrittr)

#
# Load QC filtered data
#
load("data/scData_filtered_2.Rda")

#
# Load TF motif variability data
#
load("data/dev_cluster_specific_2.Rda")
load("data/dev_clustered_cluster_specific_2.Rda")

load("data/variability_cluster_specific_2.Rda")
load("data/variability_clustered_cluster_specific_2.Rda")

#
# Get variability data
#
var_cutoff <- 1.5
n_var_cutoff <- which(variability$variability > var_cutoff)
var_data <- chromVAR::deviationScores(dev)[rownames(variability)[n_var_cutoff], ]

rownames(var_data) <- sapply(rownames(var_data), function(s) {
  s <- tolower(gsub("_MA[0-9]+", "", s))
  l <- unlist(strsplit(s, ""))
  paste0(toupper(l[1]), paste0(l[2:length(l)], collapse = ""))
})

#
# Use TF families to split heatmap rows
#
tf_distance <- as.dist(read.table("data/teichmann_similarity.txt", sep = "\t", header = T, row.names = 1))

hc <- hclust(d = tf_distance, method = "average")
tf_clusters <- factor(cutree(hc, h = 0.6))

m <- match(rownames(variability)[n_var_cutoff], names(tf_clusters))
split <- droplevels(tf_clusters[m])

split <- factor(split, levels = c("72", "50", "51", "52", "63", "46", "157", "160", "161", "180", "56", "75", "76", "21", "82", "94", "78", "107", "186", "57"), ordered = T)

#
# Generate Heatmaps
#
heatmaps <- list("reference" = list(branch = "cardiac", cluster = 2),
                 "intermediate" = list(branch = "cl3", cluster = 3),
                 "endo" = list(branch = "endo", cluster = 5),
                 "cardiac" = list(branch = "cardiac", cluster = 1),
                 "cluster4" = list(branch = "cl4", cluster = 4))

col_fun_atac <- circlize::colorRamp2(seq(from = -1, to = 5, length.out = 9), colors = RColorBrewer::brewer.pal(9, "YlGnBu"))

htlist_input <- lapply(heatmaps, function(x) {
  i <- which(SummarizedExperiment::colData(scData_filtered)$.cluster_5 == x$cluster)

  var_data_subset <- var_data[, i]
  if(!is.na(x$branch)) {
    var_data_subset <- var_data_subset[, order(SummarizedExperiment::colData(scData_filtered)[i, paste0("dpt_", x$branch)])]
  }
  var_data_smooth <- t(zoo::rollapply(t(var_data_subset), width = 9, FUN = mean, fill = "extend", by.column = T))

  ComplexHeatmap::Heatmap(var_data_smooth,
                          col = col_fun_atac,
                          cluster_columns = F,
                          split = split,
                          cluster_rows = T,
                          show_column_names = F,
                          show_row_names = ifelse(x$cluster == 2, TRUE, FALSE),
                          row_names_side = "left",
                          show_row_dend = F,
                          name = paste("Cluster", x$cluster),
                          column_title_side = "top",
                          column_title = paste("Cluster", x$cluster),
                          show_heatmap_legend = ifelse(x$cluster == 2, TRUE, FALSE),
                          row_names_gp = grid::gpar(fontsize = 6),
                          heatmap_legend_param = list(legend_direction = "horizontal", at = c(-1, 5), labels = c("low", "high"), title = "Motif accessibility of TF"))
})

htlist <- Reduce(ComplexHeatmap::add_heatmap, htlist_input)
plot_chromvar_dpt <- grid::grid.grabExpr(ComplexHeatmap::draw(htlist, column_title_gp = grid::gpar(fontface = "bold"), heatmap_legend_side = "bottom"))

figure_chromvar_dpt <- multipanelfigure::multi_panel_figure(rows = 4, columns = 1, height = 230, width = 205)
figure_chromvar_dpt <- multipanelfigure::fill_panel(figure_chromvar_dpt, plot_chromvar_dpt, row = 2:4)

ggplot2::ggsave(figure_chromvar_dpt,
                file = "figures_2/chromVar_dpt_heatmap.pdf",
                height = 230,
                width = 205,
                unit = "mm")

