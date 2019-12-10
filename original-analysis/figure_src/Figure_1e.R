library(RColorBrewer)
library(ComplexHeatmap)
library(circlize)
library(rioja)
library(gridExtra)
library(cowplot)
library(multipanelfigure)

#
# Load data
#
source("src/parameters.R")
do_zscore <- T
n_top_genes <- 30
load(file.path(parameters$general$path_supdata, "c1_subset.Rdata"))

#
# Nkx2-5 differential expression heatmap
#
nkx_cells <- which(pData(c1_subset)$Background == "nkx2-5" & pData(c1_subset)$Platform == "C1")
nkx_de_genes <- rownames(fData(c1_subset))[which(fData(c1_subset)$nkx_de & fData(c1_subset)$nkx_marker)]

# Find top 30 genes per cluster
load(file.path(parameters$general$path_supdata, "Nkx2-5-diffExprs.Rdata"))
load(file.path(parameters$general$path_supdata, "Nkx2-5-markers.Rdata"))

nkx_de_gene_order <- lapply(names(nkx_diff_data), function(n) {
  data <- nkx_diff_data[[n]]
  j <- which(data$lfc.hi < -2)
  data <- data[j, ]
  m <- na.omit(match(nkx_de_genes, data$geneID))
  i <- order(abs(data$lfc.hi[m]), decreasing = T)
  data$geneID[m][i[1:min(n_top_genes, length(m))]]
})

nkx_de_selected <- unique(unlist(nkx_de_gene_order))
nkx_subset <- c1_subset[nkx_de_selected, nkx_cells]

#
# Constrained hierarchical clustering
#
nkx_subset <- nkx_subset[, !is.na(pData(nkx_subset)$cluster)]

# Get clustering to order cells
nkx_cl_order <- order(pData(nkx_subset)$cluster)
nkx_expr_data <- get_exprs(nkx_subset, "norm_exprs_sf")[, nkx_cl_order]
nkx_dist <- as.dist(1 - cor(nkx_expr_data, method="pearson"))

nkx_const_clust <- chclust(nkx_dist)
nkx_const_clust_dendro <- as.dendrogram(nkx_const_clust)

rownames(nkx_expr_data) <- fData(nkx_subset)$symbol

# Zscore
if (do_zscore) {
  nkx_expr_data <- t(scale(t(nkx_expr_data), scale = F))
}

#
# Isl1 differential expression heatmap
#
isl1_cells <- which(pData(c1_subset)$Background == "isl1" & pData(c1_subset)$Platform == "C1")
isl1_de_genes <- rownames(fData(c1_subset))[which(fData(c1_subset)$isl1_de & fData(c1_subset)$isl1_marker)]

# Find top 30 genes per cluster
load(file.path(parameters$general$path_supdata, "Isl1-diffExprs.Rdata"))
load(file.path(parameters$general$path_supdata, "Isl1-markers.Rdata"))

isl1_de_gene_order <- lapply(names(isl1_diff_data), function(n) {
  data <- isl1_diff_data[[n]]
  j <- which(data$lfc.hi < -2)
  data <- data[j, ]
  m <- na.omit(match(isl1_de_genes, data$geneID))
  i <- order(abs(data$lfc.hi[m]), decreasing = T)
  na.omit(data$geneID[m][i[1:min(n_top_genes, length(m))]])
})

isl1_de_selected <- na.omit(unique(unlist(isl1_de_gene_order)))
isl1_subset <- c1_subset[isl1_de_selected, isl1_cells]

#
# Constrained hierarchical clustering
#
isl1_subset <- isl1_subset[, !is.na(pData(isl1_subset)$cluster)]

# Get clustering to order cells
isl1_cl_order <- order(pData(isl1_subset)$cluster)
isl1_expr_data <- get_exprs(isl1_subset, "norm_exprs_sf")[, isl1_cl_order]
isl1_dist <- as.dist(1 - cor(isl1_expr_data, method="pearson"))

isl1_const_clust <- chclust(isl1_dist)
isl1_const_clust_dendro <- as.dendrogram(isl1_const_clust)

rownames(isl1_expr_data) <- fData(isl1_subset)$symbol

# Zscore
if (do_zscore) {
  isl1_expr_data <- t(scale(t(isl1_expr_data), scale = F))
}

#
# Color mapping function to bring both heatmaps to same scale
#
e_min <- min(nkx_expr_data, isl1_expr_data)
e_max <- max(nkx_expr_data, isl1_expr_data)

bounds <- max(abs(e_min), abs(e_max))

breaks <- seq(-bounds, bounds, length.out = 99)
colors <- colorRampPalette(rev(brewer.pal(10, "RdBu")))(99)

col_fun <- colorRamp2(breaks, colors)

#
# Nkx2-5 Heatmap annotation and Heatmap
#
nkx_cell_annotation_data <- data.frame(Timepoint = pData(nkx_subset)$Timepoint[nkx_cl_order], Cluster = pData(nkx_subset)$cluster[nkx_cl_order])
nkx_cell_annotation <- HeatmapAnnotation(nkx_cell_annotation_data, col = list(Timepoint = c("e7.5" = "#FF003C", "e8.5" = "#248F8D", "e9.5" = "#987F69"),
                                                                              Cluster = c("3" = "#00008B", "1" = "#8B0000", "2" = "#CDAD00")),
                                         show_legend = F,
                                         show_annotation_name = T,
                                         annotation_name_gp = gpar(fontsize = 4),
                                         annotation_name_offset = unit(0.8, "mm"),
                                         gap = unit(0.1, "mm"),
                                         annotation_height = unit.c(unit(1.5, "mm"), unit(1.5, "mm")))

nkx_de_heatmap <- Heatmap(nkx_expr_data,
                          col = col_fun,
                          top_annotation = nkx_cell_annotation,
                          cluster_columns = nkx_const_clust_dendro,
                          cluster_rows = T,
                          show_column_names = F,
                          show_row_dend = F,
                          show_column_dend = F,
                          name = "nkx",
                          #column_title = "Nkx2-5+ cells",
                          column_title_gp = gpar(fontface = "bold", fontsize = 8),
                          split = rep(paste("Marker genes\ncluster", c(1,2,3)), times = sapply(nkx_de_gene_order, length)),
                          gap = unit(0.3, "mm"),
                          row_title_rot = 90,
                          row_title_gp = gpar(fontsize = 4),
                          row_names_gp = gpar(fontsize = 3),
                          heatmap_legend_param = list("plot" = F, "color_bar" = "continuous", "legend_direction" = "horizontal", "title" = "Scaled expression", "title_position" = "leftcenter", "title_gp" = gpar(fontsize = 4), "labels_gp" = gpar(fontsize = 4), "grid_height" = unit(2, "mm")))

#
# Isl1 Heatmap annotation and Heatmap
#
isl1_cell_annotation_data <- data.frame(Timepoint = pData(isl1_subset)$Timepoint[isl1_cl_order], Cluster = pData(isl1_subset)$cluster[isl1_cl_order])
isl1_cell_annotation <- HeatmapAnnotation(isl1_cell_annotation_data, col = list(Timepoint = c("e7.5" = "#FF003C", "e8.5" = "#248F8D", "e9.5" = "#987F69"),
                                                                                Cluster = c("1" = "#1B9E77", "2" = "#D95F02", "3" = "#7570B3", "4" = "#E7298A", "5" = "#66A61E")),
                                          show_legend = F,
                                          show_annotation_name = T,
                                          annotation_name_gp = gpar(fontsize = 4),
                                          annotation_name_offset = unit(0.8, "mm"),
                                          gap = unit(0.1, "mm"),
                                          annotation_height = unit.c(unit(1.5, "mm"), unit(1.5, "mm")))

isl1_de_heatmap <- Heatmap(isl1_expr_data,
                           col = col_fun,
                           top_annotation = isl1_cell_annotation,
                           cluster_columns = isl1_const_clust_dendro,
                           cluster_rows = T,
                           show_column_names = F,
                           show_row_dend = F,
                           show_column_dend = F,
                           name = "isl1",
                           #column_title = "Isl1+ cells",
                           column_title_gp = gpar(fontface = "bold", fontsize = 8),
                           split = rep(paste("Marker genes\ncluster", c(1,2,3,4,5)), times = sapply(isl1_de_gene_order, length)),
                           gap = unit(0.3, "mm"),
                           row_title_rot = 90,
                           row_title_gp = gpar(fontsize = 4),
                           row_names_gp = gpar(fontsize = 3),
                           heatmap_legend_param = list("plot" = T ,"color_bar" = "continuous", "legend_direction" = "horizontal", "title" = "Scaled expression", "title_position" = "leftcenter", "title_gp" = gpar(fontsize = 4), "labels_gp" = gpar(fontsize = 4), "grid_height" = unit(2, "mm")))

#
# Build heatmap for Figure 1e
#
fig_1e <- grid.arrange(grid.grabExpr(draw(nkx_de_heatmap, heatmap_legend_side = "bottom", padding = unit(c(0, 0, 0, 0), "mm"))), grid.grabExpr(draw(isl1_de_heatmap, heatmap_legend_side = "bottom", padding = unit(c(0, 0, -5, 0), "mm"))), nrow = 2)

m <- multi_panel_figure(width = 205, height = 230, columns = 4, rows = 4)
m <- fill_panel(m, fig_1e, label = "E", column = c(1,2), row = c(2,3,4))

ggsave(
  plot = m,
  filename = file.path(parameters$general$path_rfigures, "Figure_1_e2.pdf"),
  width = 205,
  height = 230,
  units = "mm",
  dpi = 600)

#
# Heatmap for all Nkx2-5 DE genes
#
nkx_de_gene_order <- lapply(names(nkx_diff_data), function(n) {
  data <- nkx_diff_data[[n]]
  j <- which(data$lfc.hi < -1)
  data <- data[j, ]
  m <- na.omit(match(nkx_de_genes, data$geneID))
  i <- order(abs(data$lfc.hi[m]), decreasing = T)
  data$geneID[m]
})

nkx_de_selected <- unique(unlist(nkx_de_gene_order))
nkx_subset <- c1_subset[nkx_de_selected, nkx_cells]

#
# Constrained hierarchical clustering
#
nkx_subset <- nkx_subset[, !is.na(pData(nkx_subset)$cluster)]

# Get clustering to order cells
nkx_cl_order <- order(pData(nkx_subset)$cluster)
nkx_expr_data <- get_exprs(nkx_subset, "norm_exprs_sf")[, nkx_cl_order]
nkx_dist <- as.dist(1 - cor(nkx_expr_data, method="pearson"))

nkx_const_clust <- chclust(nkx_dist)
nkx_const_clust_dendro <- as.dendrogram(nkx_const_clust)

rownames(nkx_expr_data) <- fData(nkx_subset)$symbol

# Zscore
if (do_zscore) {
  nkx_expr_data <- t(scale(t(nkx_expr_data), scale = F))
}

# Heatmap
nkx_all_de_heatmap <- Heatmap(nkx_expr_data,
                          col = colorRampPalette(rev(brewer.pal(10, "RdBu")))(99),
                          top_annotation = nkx_cell_annotation,
                          cluster_columns = nkx_const_clust_dendro,
                          cluster_rows = T,
                          show_column_names = F,
                          show_row_dend = F,
                          name = "nkx",
                          column_title = "Nkx2-5+ cells",
                          column_title_gp = gpar(fontface = "bold", fontsize = 8),
                          split = rep(paste("Marker genes\ncluster", c(1,2,3)), times = sapply(nkx_de_gene_order, length)),
                          gap = unit(0.5, "mm"),
                          row_title_rot = 90,
                          row_title_gp = gpar(fontsize = 4),
                          row_names_gp = gpar(fontsize = 1),
                          heatmap_legend_param = list("color_bar" = "continuous", "legend_direction" = "horizontal", "title" = "Scaled expression", "title_position" = "leftcenter", "title_gp" = gpar(fontsize = 4), "labels_gp" = gpar(fontsize = 4), "grid_height" = unit(2, "mm")))

#
# Heatmap for all Isl1 DE genes
#
isl1_de_gene_order <- lapply(names(isl1_diff_data), function(n) {
  data <- isl1_diff_data[[n]]
  j <- which(data$lfc.hi < -2)
  data <- data[j, ]
  m <- na.omit(match(isl1_de_genes, data$geneID))
  i <- order(abs(data$lfc.hi[m]), decreasing = T)
  na.omit(data$geneID[m])
})
# Remove Etv2 from cluster 5
isl1_de_gene_order[[5]] <- isl1_de_gene_order[[5]][! isl1_de_gene_order[[5]] %in% c("ENSMUSG00000006311")]
# Rescue two missing genes
isl1_de_gene_order[[2]] <- c(isl1_de_gene_order[[2]], "ENSMUSG00000029432")
isl1_de_gene_order[[5]] <- c(isl1_de_gene_order[[5]], "ENSMUSG00000050953")
isl1_de_selected <- na.omit(unique(unlist(isl1_de_gene_order)))
isl1_subset <- c1_subset[isl1_de_selected, isl1_cells]

#
# Constrained hierarchical clustering
#
isl1_subset <- isl1_subset[, !is.na(pData(isl1_subset)$cluster)]

# Get clustering to order cells
isl1_cl_order <- order(pData(isl1_subset)$cluster)
isl1_expr_data <- get_exprs(isl1_subset, "norm_exprs_sf")[, isl1_cl_order]
isl1_dist <- as.dist(1 - cor(isl1_expr_data, method="pearson"))

isl1_const_clust <- chclust(isl1_dist)
isl1_const_clust_dendro <- as.dendrogram(isl1_const_clust)

rownames(isl1_expr_data) <- fData(isl1_subset)$symbol

# Zscore
if (do_zscore) {
  isl1_expr_data <- t(scale(t(isl1_expr_data), scale = F))
}

# Heatmap
isl1_all_de_heatmap <- Heatmap(isl1_expr_data,
                           col = rev(brewer.pal(10, "RdBu")),
                           top_annotation = isl1_cell_annotation,
                           cluster_columns = isl1_const_clust_dendro,
                           cluster_rows = T,
                           show_column_names = F,
                           show_row_dend = F,
                           name = "isl1",
                           column_title = "Isl1+ cells",
                           column_title_gp = gpar(fontface = "bold", fontsize = 8),
                           split = rep(paste("Marker genes\ncluster", c(1,2,3,4,5)), times = sapply(isl1_de_gene_order, length)),
                           gap = unit(0.5, "mm"),
                           row_title_rot = 90,
                           row_title_gp = gpar(fontsize = 4),
                           row_names_gp = gpar(fontsize = 1),
                           heatmap_legend_param = list("color_bar" = "continuous", "legend_direction" = "horizontal", "title" = "Scaled expression", "title_position" = "leftcenter", "title_gp" = gpar(fontsize = 4), "labels_gp" = gpar(fontsize = 4), "grid_height" = unit(2, "mm")))

#
# Supplementary Figure X
#
sup_figX_a <- grid.grabExpr(draw(nkx_all_de_heatmap, heatmap_legend_side = "bottom"))
sup_figX_b <- grid.grabExpr(draw(isl1_all_de_heatmap, heatmap_legend_side = "bottom"))

m <- multi_panel_figure(width = 205, height = 230, columns = 1, rows = 2)
m <- fill_panel(m, sup_figX_a, row = 1)
m <- fill_panel(m, sup_figX_b, row = 2)

ggsave(
  plot = m,
  filename = file.path(parameters$general$path_rfigures, "Supplementary_Figure_DE_genes.pdf"),
  width = 205,
  height = 230,
  units = "mm",
  dpi = 600)
