library(scater)
library(dplyr)
library(tidyr)
library(singlecellutils)
library(RColorBrewer)
library(viridis)
library(multipanelfigure)
library(ComplexHeatmap)
library(circlize)
library(zoo)
library(MPAgenomics)
library(gridExtra)

#
# Load data
#
source("src/parameters.R")
min_cor_cluster <- 0.5
min_cor <- 0.7
load(file.path(parameters$general$path_supdata, "c1_subset.Rdata"))

#
# Nkx2-5 Rank based correlation
#
nkx_cells <- which(pData(c1_subset)$Background == "nkx2-5" & pData(c1_subset)$Platform == "C1")
nkx_de_genes <- rownames(fData(c1_subset))[which(fData(c1_subset)$nkx_de & fData(c1_subset)$nkx_marker)]
nkx_cluster <- factor(pData(c1_subset[, nkx_cells])$cluster)

nkx_data <- get_exprs(c1_subset[nkx_de_genes, nkx_cells], "norm_exprs_sf")

nkx_cor <- cor(t(nkx_data), method = "pearson")
nkx_dpt <- pData(c1_subset[, nkx_cells])$dpt

nkx_cor_list <- lapply(levels(nkx_cluster), function(c) {
  ind <- which(nkx_cluster == c)
  cc <- cor(t(nkx_data[, ind]), cbind(nkx_dpt[ind]), method = "spearman")
})
nkx_cor_list$dpt <- cor(t(nkx_data), cbind(nkx_dpt), method = "spearman")

nkx_correlation <- as.data.frame(do.call("cbind", nkx_cor_list))
colnames(nkx_correlation) <- c(paste0("cluster", levels(nkx_cluster)), "nkx_dpt")

write.table(nkx_correlation, file = file.path(parameters$general$path_supdata,"Nkx2-5-correlation.txt"), quote = F, row.names = T, col.names = T, sep = "\t")

#
# Extract highly correlated genes
#
# From cluster correlation
nkx_i <- unique(unlist(apply(nkx_correlation[, -ncol(nkx_correlation)], 2, function(x) which(abs(x) > min_cor_cluster))))

# From dpt correlation
nkx_i <- unique(c(nkx_i, which(abs(nkx_correlation$nkx_dpt) > min_cor)))

nkx_cor_data <- nkx_cor[nkx_i, nkx_i]
colnames(nkx_cor_data) <- fData(c1_subset[nkx_de_genes, ][nkx_i, ])$symbol
rownames(nkx_cor_data) <- fData(c1_subset[nkx_de_genes, ][nkx_i, ])$symbol

#
# Isl1 Rank based correlation
#
isl1_cells <- which(pData(c1_subset)$Background == "isl1" & pData(c1_subset)$Platform == "C1")
isl1_de_genes <- rownames(fData(c1_subset))[which(fData(c1_subset)$isl1_de & fData(c1_subset)$isl1_marker)]
isl1_cluster <- factor(pData(c1_subset[, isl1_cells])$cluster)

isl1_data <- get_exprs(c1_subset[isl1_de_genes, isl1_cells], "norm_exprs_sf")

isl1_cor <- cor(t(isl1_data), method = "pearson")
isl1_dpt <- pData(c1_subset[, isl1_cells])$dpt

isl1_cor_list <- lapply(levels(isl1_cluster), function(c) {
  ind <- which(isl1_cluster == c)
  cc <- cor(t(isl1_data[, ind]), cbind(isl1_dpt[ind]), method = "spearman")
})
isl1_cor_list$dpt <- cor(t(isl1_data), cbind(isl1_dpt), method = "spearman")

isl1_correlation <- as.data.frame(do.call("cbind", isl1_cor_list))
colnames(isl1_correlation) <- c(paste0("cluster", levels(isl1_cluster)), "isl1_dpt")

write.table(isl1_correlation, file = file.path(parameters$general$path_supdata,"Isl1-correlation.txt"), quote = F, row.names = T, col.names = T, sep = "\t")

#
# Extract highly correlated genes
#
# From cluster correlation
isl1_i <- unique(unlist(apply(isl1_correlation[, -ncol(isl1_correlation)], 2, function(x) which(abs(x) > min_cor_cluster))))

# From dpt correlation
isl1_i <- unique(c(isl1_i, which(abs(isl1_correlation$isl1_dpt) > min_cor)))

isl1_cor_data <- isl1_cor[isl1_i, isl1_i]
colnames(isl1_cor_data) <- fData(c1_subset[isl1_de_genes, ][isl1_i, ])$symbol
rownames(isl1_cor_data) <- fData(c1_subset[isl1_de_genes, ][isl1_i, ])$symbol

#
# Nkx2-5 Heatmaps
#

nkx_cor_heatmap <- Heatmap(nkx_cor_data,
              col = colorRamp2(seq(from = -1, to = 1, by = 0.25), rev(brewer.pal(9, "RdBu"))),
              row_names_gp = gpar(fontsize = 2),
              column_names_gp = gpar(fontsize = 2),
              show_column_dend = F,
              column_title = "Nkx2-5+ gene correlation",
              column_title_gp = gpar(fontface = "bold", fontsize = 8),
              name = "nkx_cor",
              heatmap_legend_param = list("at" = c(-1, 0, 1), "color_bar" = "continuous", "legend_direction" = "horizontal", "title" = "Pearson correlation", "title_position" = "leftcenter", "title_gp" = gpar(fontsize = 4), "labels_gp" = gpar(fontsize = 4), "grid_height" = unit(2, "mm")))

write.table(file = file.path(parameters$general$path_supdata, "Source_Data_SupplementaryFigure_7a.txt"),
            x = nkx_cor_data,
            row.names = T,
            col.names = T,
            sep = "\t",
            quote = F)
write.table(file = file.path(parameters$general$path_supdata, "Source_Data_SupplementaryFigure_7b.txt"),
            x = isl1_cor_data,
            row.names = T,
            col.names = T,
            sep = "\t",
            quote = F)

#
# Smoothed expression heatmap
#
nkx_data <- nkx_data[nkx_i,]
rownames(nkx_data) <- fData(c1_subset[nkx_de_genes, ][nkx_i, ])$symbol

nkx_pseudotime_order <- order(nkx_dpt)
nkx_data_smoothed <- t(apply(nkx_data[, nkx_pseudotime_order], 1, function(r) rollapply(r, 11, mean, fill="extend")))

write.table(file = file.path(parameters$general$path_supdata, "Source_Data_SupplementaryFigure_7c.txt"),
            x = nkx_data_smoothed,
            row.names = T,
            col.names = T,
            sep = "\t",
            quote = F)

# Annotation of columns
nkx_cell_annotation_data <- data.frame(Timepoint = pData(c1_subset[, nkx_cells])$Timepoint[nkx_pseudotime_order], Cluster = pData(c1_subset[, nkx_cells])$cluster[nkx_pseudotime_order], Pseudotime = nkx_dpt[nkx_pseudotime_order])
nkx_cell_annotation <- HeatmapAnnotation(nkx_cell_annotation_data, col = list(Timepoint = c("e7.5" = "#36648B", "e8.5" = "#4F94CD", "e9.5" = "#63B8FF"),
                                                                              Cluster = c("1" = "#00008B", "2" = "#8B0000", "3" = "#CDAD00"),
                                                                              Pseudotime = colorRamp2(c(min(nkx_dpt), mean(nkx_dpt), max(nkx_dpt)), colors = c("#36648B", "#4F94CD", "#63B8FF"))),
                                         show_legend = F,
                                         show_annotation_name = T,
                                         annotation_name_gp = gpar(fontsize = 4),
                                         annotation_name_offset = unit(0.8, "mm"),
                                         gap = unit(0.1, "mm"),
                                         annotation_height = unit.c(unit(1.5, "mm"), unit(1.5, "mm"), unit(1.5, "mm")))

nkx_smooth_heatmap <- Heatmap(nkx_data_smoothed,
                              col = inferno(99),
                              cluster_columns = F,
                              #clustering_distance_rows = 'manhattan',
                              clustering_method_rows = 'ward.D',
                              #split = cl$clustering,
                              row_names_gp = gpar(fontsize = 1.5),
                              column_names_gp = gpar(fontsize = 1.5),
                              show_column_names = F,
                              name = "nkx_smooth",
                              top_annotation = nkx_cell_annotation,
                              heatmap_legend_param = list("color_bar" = "continuous", "legend_direction" = "horizontal", "title" = "Expression", "title_position" = "leftcenter", "title_gp" = gpar(fontsize = 4), "labels_gp" = gpar(fontsize = 4), "grid_height" = unit(2, "mm")))

#
# Gene classification into "upregulated", "downregulated" and "other"
#
nkx_segments <- factor(apply(nkx_data_smoothed, 1, function(g) {
  smntn <- segmentation(as.numeric(g), method = "PELT", plot = F, verbose = F)
  s <- smntn$segment
  s$level <- ifelse(s$means < 1, "off", ifelse(s$means < 5, "low", ifelse(s$means < 8, "moderate", ifelse(s$means < 10, "high", "very high"))))
  s$level <- factor(s$level, levels = c("off", "low", "moderate", "high", "very high"), ordered = T)
  n <- nrow(s)

  direction <- paste(sapply(2:nrow(s), function(i) {
    diff <- as.numeric(s$level[i]) - as.numeric(s$level[i-1])
    if(diff >= 1) return("up")
    if(diff <= -1) return("down")
  }), collapse = "-")

  class = "other"
  if(nrow(s) > 1) {
    direction <- paste(sapply(2:nrow(s), function(i) {
      diff <- as.numeric(s$level[i]) - as.numeric(s$level[i-1])
      if(diff >= 1) return("up")
      if(diff <= -1) return("down")
    }), collapse = "-")
    if (startsWith(direction, "up") & endsWith(direction, "up")) class <- "upregulated"
    if (startsWith(direction, "down") & endsWith(direction, "down")) class <- "downregulated"
    #if (startsWith(direction, "up") & endsWith(direction, "down")) class <- "intermediate"
    #if (startsWith(direction, "up") & endsWith(direction, "up") & s$level[1] > "off") class <- "upregulated - primed"
    #if (startsWith(direction, "up") & endsWith(direction, "up") & s$level[1] == "off") class <- "upregulated - de novo"
    #if (startsWith(direction, "down") & endsWith(direction, "down")) class <- "downregulated"
  }
  return(class)
}), levels = c("downregulated", "intermediate", "upregulated", "upregulated - primed", "upregulated - de novo", "other"), ordered = T)

nkx_smooth_heatmap_class <- Heatmap(nkx_data_smoothed,
                                 col = inferno(99),
                                 cluster_columns = F,
                                 split = nkx_segments,
                                 row_names_gp = gpar(fontsize = 3),
                                 column_names_gp = gpar(fontsize = 1.5),
                                 row_title_gp = gpar(fontsize = 4),
                                 gap = unit(0.5, "mm"),
                                 show_column_names = F,
                                 show_row_dend = F,
                                 name = "nkx_smooth",
                                 top_annotation = nkx_cell_annotation,
                                 heatmap_legend_param = list("color_bar" = "continuous", "legend_direction" = "horizontal", "title" = "Expression", "title_position" = "leftcenter", "title_gp" = gpar(fontsize = 4), "labels_gp" = gpar(fontsize = 4), "grid_height" = unit(2, "mm")))

#
# Isl1 heatmaps
#

isl1_cor_heatmap <- Heatmap(isl1_cor_data,
                           col = colorRamp2(seq(from = -1, to = 1, by = 0.25), rev(brewer.pal(9, "RdBu"))),
                           row_names_gp = gpar(fontsize = 2),
                           column_names_gp = gpar(fontsize = 2),
                           show_column_dend = F,
                           column_title = "Isl1+ gene correlation",
                           column_title_gp = gpar(fontface = "bold", fontsize = 8),
                           name = "nkx_cor",
                           heatmap_legend_param = list("at" = c(-1, 0, 1), "color_bar" = "continuous", "legend_direction" = "horizontal", "title" = "Pearson correlation", "title_position" = "leftcenter", "title_gp" = gpar(fontsize = 4), "labels_gp" = gpar(fontsize = 4), "grid_height" = unit(2, "mm")))

#
# Smoothed expression heatmap
#
isl1_data <- isl1_data[isl1_i,]
rownames(isl1_data) <- fData(c1_subset[isl1_de_genes, ][isl1_i, ])$symbol

isl1_pseudotime_order <- order(isl1_dpt)
isl1_data_smoothed <- t(apply(isl1_data[, isl1_pseudotime_order], 1, function(r) rollapply(r, 11, mean, fill="extend")))

isl1_branch1 <- which(isl1_cluster[isl1_pseudotime_order] %in% c(2,3,4,5))
isl1_branch2 <- which(isl1_cluster[isl1_pseudotime_order] %in% c(1,5))

write.table(file = file.path(parameters$general$path_supdata, "Source_Data_SupplementaryFigure_7d.txt"),
            x = isl1_data_smoothed,
            row.names = T,
            col.names = T,
            sep = "\t",
            quote = F)

#
# Color mapping function to bring both heatmaps to same scale
#
e_min <- min(isl1_data_smoothed)
e_max <- max(isl1_data_smoothed)

breaks <- seq(e_min, e_max, length.out = 99)
colors <- inferno(99)

col_fun <- colorRamp2(breaks, colors)


# Annotation of columns
isl1_cell_annotation_data <- data.frame(Timepoint = pData(c1_subset[, isl1_cells])$Timepoint[isl1_pseudotime_order], Cluster = pData(c1_subset[, isl1_cells])$cluster[isl1_pseudotime_order], Pseudotime = isl1_dpt[isl1_pseudotime_order])
isl1_cell_annotation_data_branch1 <- isl1_cell_annotation_data[isl1_branch1, ]
isl1_cell_annotation_data_branch2 <- isl1_cell_annotation_data[isl1_branch2, ]

isl1_cell_annotation <- HeatmapAnnotation(isl1_cell_annotation_data, col = list(Timepoint = c("e7.5" = "#008B45", "e8.5" = "#00CD66", "e9.5" = "#00FF7F"),
                                                                              Cluster = c("1" = "#1B9E77", "2" = "#D95F02", "3" = "#7570B3", "4" = "#E7298A", "5" = "#66A61E"),
                                                                              Pseudotime = colorRamp2(c(min(isl1_dpt), mean(isl1_dpt), max(isl1_dpt)), colors = c("#008B45", "#00CD66", "#00FF7F"))),
                                         show_legend = F,
                                         show_annotation_name = T,
                                         annotation_name_gp = gpar(fontsize = 4),
                                         annotation_name_offset = unit(0.8, "mm"),
                                         gap = unit(0.1, "mm"),
                                         annotation_height = unit.c(unit(1.5, "mm"), unit(1.5, "mm"), unit(1.5, "mm")))

isl1_cell_annotation_branch1 <- HeatmapAnnotation(isl1_cell_annotation_data_branch1, col = list(Timepoint = c("e7.5" = "#008B45", "e8.5" = "#00CD66", "e9.5" = "#00FF7F"),
                                                                                Cluster = c("1" = "#1B9E77", "2" = "#D95F02", "3" = "#7570B3", "4" = "#E7298A", "5" = "#66A61E"),
                                                                                Pseudotime = colorRamp2(c(min(isl1_dpt), mean(isl1_dpt), max(isl1_dpt)), colors = c("#008B45", "#00CD66", "#00FF7F"))),
                                          show_legend = F,
                                          show_annotation_name = T,
                                          annotation_name_gp = gpar(fontsize = 4),
                                          annotation_name_offset = unit(0.8, "mm"),
                                          gap = unit(0.1, "mm"),
                                          annotation_height = unit.c(unit(1.5, "mm"), unit(1.5, "mm"), unit(1.5, "mm")))

isl1_cell_annotation_branch2 <- HeatmapAnnotation(isl1_cell_annotation_data_branch2, col = list(Timepoint = c("e7.5" = "#008B45", "e8.5" = "#00CD66", "e9.5" = "#00FF7F"),
                                                                                Cluster = c("1" = "#1B9E77", "2" = "#D95F02", "3" = "#7570B3", "4" = "#E7298A", "5" = "#66A61E"),
                                                                                Pseudotime = colorRamp2(c(min(isl1_dpt), mean(isl1_dpt), max(isl1_dpt)), colors = c("#008B45", "#00CD66", "#00FF7F"))),
                                          show_legend = F,
                                          show_annotation_name = F,
                                          annotation_name_gp = gpar(fontsize = 4),
                                          annotation_name_offset = unit(0.8, "mm"),
                                          gap = unit(0.1, "mm"),
                                          annotation_height = unit.c(unit(1.5, "mm"), unit(1.5, "mm"), unit(1.5, "mm")))

isl1_smooth_heatmap <- Heatmap(isl1_data_smoothed,
                              col = col_fun,
                              cluster_columns = F,
                              row_names_gp = gpar(fontsize = 4),
                              column_names_gp = gpar(fontsize = 4),
                              show_column_names = F,
                              name = "isl1_smooth",
                              top_annotation = isl1_cell_annotation,
                              heatmap_legend_param = list("color_bar" = "continuous", "legend_direction" = "horizontal", "title" = "Expression", "title_position" = "leftcenter", "title_gp" = gpar(fontsize = 4), "labels_gp" = gpar(fontsize = 4), "grid_height" = unit(2, "mm")))

isl1_smooth_heatmap_branch1 <- Heatmap(isl1_data_smoothed[, isl1_branch1],
                               col = col_fun,
                               cluster_columns = F,
                               row_names_gp = gpar(fontsize = 1.5),
                               column_names_gp = gpar(fontsize = 1.5),
                               show_column_names = F,
                               name = "isl1_smooth_b1",
                               top_annotation = isl1_cell_annotation_branch1,
                               heatmap_legend_param = list("color_bar" = "continuous", "legend_direction" = "horizontal", "title" = "Expression", "title_position" = "leftcenter", "title_gp" = gpar(fontsize = 4), "labels_gp" = gpar(fontsize = 4), "grid_height" = unit(2, "mm")))

isl1_smooth_heatmap_branch2 <- Heatmap(isl1_data_smoothed[, isl1_branch2],
                               col = col_fun,
                               cluster_columns = F,
                               row_names_gp = gpar(fontsize = 1.5),
                               column_names_gp = gpar(fontsize = 1.5),
                               show_column_names = F,
                               show_row_names = F,
                               name = "isl1_smooth_b2",
                               top_annotation = isl1_cell_annotation_branch2,
                               heatmap_legend_param = list("color_bar" = "continuous", "legend_direction" = "horizontal", "title" = "Expression", "title_position" = "leftcenter", "title_gp" = gpar(fontsize = 4), "labels_gp" = gpar(fontsize = 4), "grid_height" = unit(2, "mm")))

# CLassification
isl1_segments <- factor(apply(isl1_data_smoothed[, isl1_branch1], 1, function(g) {
  smntn <- segmentation(as.numeric(g), method = "PELT", plot = F, verbose = F)
  s <- smntn$segment
  s$level <- ifelse(s$means < 1, "off", ifelse(s$means < 5, "low", ifelse(s$means < 8, "moderate", ifelse(s$means < 10, "high", "very high"))))
  s$level <- factor(s$level, levels = c("off", "low", "moderate", "high", "very high"), ordered = T)
  n <- nrow(s)

  class = "other"
  if(nrow(s) > 1) {
    direction <- paste(sapply(2:nrow(s), function(i) {
      diff <- as.numeric(s$level[i]) - as.numeric(s$level[i-1])
      if(diff >= 1) return("up")
      if(diff <= -1) return("down")
    }), collapse = "-")
    if (startsWith(direction, "up") & endsWith(direction, "up")) class <- "upregulated"
    if (startsWith(direction, "down") & endsWith(direction, "down")) class <- "downregulated"
    #if (startsWith(direction, "up") & endsWith(direction, "down")) class <- "intermediate"
    #if (startsWith(direction, "up") & endsWith(direction, "up") & s$level[1] > "off") class <- "upregulated - primed"
    #if (startsWith(direction, "up") & endsWith(direction, "up") & s$level[1] == "off") class <- "upregulated - de novo"
    #if (startsWith(direction, "down") & endsWith(direction, "down")) class <- "downregulated"
  }
  return(class)
}), levels = c("downregulated", "intermediate", "upregulated", "upregulated - primed", "upregulated - de novo", "other"), ordered = T)

isl1_smooth_heatmap_class_b1 <- Heatmap(isl1_data_smoothed[, isl1_branch1],
                                    col = col_fun,
                                    cluster_columns = F,
                                    split = isl1_segments,
                                    row_names_gp = gpar(fontsize = 3),
                                    column_names_gp = gpar(fontsize = 1.5),
                                    row_title_gp = gpar(fontsize = 4),
                                    gap = unit(0.5, "mm"),
                                    show_column_names = F,
                                    show_row_dend = F,
                                    name = "isl1_smooth_b1",
                                    top_annotation = isl1_cell_annotation_branch1,
                                    heatmap_legend_param = list("color_bar" = "continuous", "legend_direction" = "horizontal", "title" = "Expression", "title_position" = "leftcenter", "title_gp" = gpar(fontsize = 4), "labels_gp" = gpar(fontsize = 4), "grid_height" = unit(2, "mm")))

isl1_smooth_heatmap_class_b2 <- Heatmap(isl1_data_smoothed[, isl1_branch2],
                                        col = col_fun,
                                        cluster_columns = F,
                                        split = isl1_segments,
                                        row_names_gp = gpar(fontsize = 3),
                                        column_names_gp = gpar(fontsize = 1.5),
                                        row_title_gp = gpar(fontsize = 4),
                                        gap = unit(0.5, "mm"),
                                        show_column_names = F,
                                        show_row_names = F,
                                        show_row_dend = F,
                                        name = "isl1_smooth_b2",
                                        top_annotation = isl1_cell_annotation_branch2,
                                        heatmap_legend_param = list("color_bar" = "continuous", "legend_direction" = "horizontal", "title" = "Expression", "title_position" = "leftcenter", "title_gp" = gpar(fontsize = 4), "labels_gp" = gpar(fontsize = 4), "grid_height" = unit(2, "mm")))


#
# Create Figures
#
# Fig2 e,f
fig_2e <- grid.grabExpr(draw(nkx_smooth_heatmap_class, heatmap_legend_side = "bottom"))
fig_2f <- grid.grabExpr(draw(isl1_smooth_heatmap_class_b2 + isl1_smooth_heatmap_class_b1, main_heatmap = "isl1_smooth_b1", heatmap_legend_side = "bottom", gap = unit(1, "mm")))

m <- multi_panel_figure(width = 205, height = 230, columns = 7, rows = 7)

m <- fill_panel(m, fig_2e, label = "E", column = c(1:3), row = c(4,5,6,7))
m <- fill_panel(m, fig_2f, label = "F", column = c(4:7), row = c(4,5,6,7))

ggsave(
  plot = m,
  filename = file.path(parameters$general$path_rfigures, "Figure_2_ef.pdf"),
  width = 205,
  height = 230,
  units = "mm",
  dpi = 600)

# Figure 2
m <- multi_panel_figure(width = 205, height = 230, columns = 7, rows = 14)

m <- fill_panel(m, nkx_ic_plot, label = "C", column = 5:7, row = 1:3)
m <- fill_panel(m, isl1_ic_plot, label = "D", column = 5:7, row = 4:6)
m <- fill_panel(m, fig_2e, label = "E", column = 1:3, row = 7:14)
m <- fill_panel(m, fig_2f, label = "F", column = 4:7, row = 7:14)

ggsave(
  plot = m,
  filename = file.path(parameters$general$path_rfigures, "Figure_2.pdf"),
  width = 205,
  height = 230,
  units = "mm",
  dpi = 600)

# Supplementary figure for classified heatmaps
sup_figX_a <- grid.grabExpr(draw(nkx_smooth_heatmap, heatmap_legend_side = "bottom"))
sup_figX_b <- grid.grabExpr(draw(isl1_smooth_heatmap_branch2 + isl1_smooth_heatmap_branch1, main_heatmap = "isl1_smooth_b1", heatmap_legend_side = "bottom", gap = unit(1, "mm"), padding = unit(c(1, 1, 1, 1), "mm")))

m <- multi_panel_figure(width = 205, height = 230, columns = 7, rows = 3)

m <- fill_panel(m, sup_figX_a, label = "A", column = c(1:3), row = 2)
m <- fill_panel(m, sup_figX_b, label = "B", column = c(4:7), row = 2)

ggsave(
  plot = m,
  filename = file.path(parameters$general$path_rfigures, "Supplementary_Figure_smoothHeatmap.pdf"),
  width = 205,
  height = 230,
  units = "mm",
  dpi = 600)

# Supplementary figure for correlation heatmaps
sup_figX_a <- grid.grabExpr(draw(nkx_cor_heatmap, heatmap_legend_side = "bottom"))
sup_figX_b <- grid.grabExpr(draw(isl1_cor_heatmap, heatmap_legend_side = "bottom"))

m <- multi_panel_figure(width = 205, height = 230, columns = 3, rows = 2)
m <- fill_panel(m, sup_figX_a, row = 1, column = c(1,2))
m <- fill_panel(m, sup_figX_b, row = 2, column = c(1,2))

ggsave(
  plot = m,
  filename = file.path(parameters$general$path_rfigures, "Supplementary_Figure_Corheatmaps.pdf"),
  width = 205,
  height = 230,
  units = "mm",
  dpi = 600)
