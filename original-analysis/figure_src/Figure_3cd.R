library(scater)
library(tidyr)
library(ggplot2)
library(singlecellutils)

#
# Load data
#
source("src/parameters.R")
load(file.path(parameters$general$path_supdata, "c1_subset.Rdata"))
load(file.path(parameters$general$path_supdata, "Isl1-DM.Rdata"))

nkx_cells <- which(pData(c1_subset)$Background == "nkx2-5" & pData(c1_subset)$Platform == "C1")
isl1_cells <- which(pData(c1_subset)$Background == "isl1" & pData(c1_subset)$Platform == "C1")

isl1_dm_coords <- cbind(-pData(c1_subset[, isl1_cells])$dm1, -pData(c1_subset[, isl1_cells])$dm2, -pData(c1_subset[, isl1_cells])$dm3)

#
# Graphical settings
#
global.settings <- list(layout.heights=list(main.key.padding = -2, bottom.padding = 0, axis.xlab.padding = -0.5), layout.widths = list(ylab.axis.padding = -0.5))

#
# Color scalings
#
scale_channel <- function(x, min, max, a, b) {
  x[x > max] <-  max
  return(ceiling(((b - a) * (x - min) / (max - min)) + a))
}

min_isl1 <- min(get_exprs(c1_subset["ENSMUSG00000042258", c(isl1_cells, nkx_cells)], "norm_exprs_sf")[1,])
max_isl1 <- max(get_exprs(c1_subset["ENSMUSG00000042258", c(isl1_cells, nkx_cells)], "norm_exprs_sf")[1,])
min_nkx <- min(get_exprs(c1_subset["ENSMUSG00000015579", c(isl1_cells, nkx_cells)], "norm_exprs_sf")[1,])
max_nkx <- max(get_exprs(c1_subset["ENSMUSG00000015579", c(isl1_cells, nkx_cells)], "norm_exprs_sf")[1,])

isl1_color <- scale_channel(get_exprs(c1_subset["ENSMUSG00000042258", c(isl1_cells, nkx_cells)], "norm_exprs_sf")[1, ], min_isl1, 7, 0, 255)
nkx_color <- scale_channel(get_exprs(c1_subset["ENSMUSG00000015579", c(isl1_cells, nkx_cells)], "norm_exprs_sf")[1, ], min_nkx, 7, 0, 255)

#
# Mix the colors of Isl1 and Nkx2-5
#
mix_color <- sapply(1:length(isl1_color), function(i) {
  red <- nkx_color[i]
  green <- isl1_color[i]
  blue <- 0
  rgb(red, green, blue, maxColorValue = 255)
})

#
# Plot clusterings from Figure 1d with coexpressing cells highlighted
#
isl1_isl1 <- ifelse(get_exprs(c1_subset["ENSMUSG00000042258", isl1_cells], "norm_exprs_sf") > 0.2, "Isl1+", "Isl1-")
isl1_nkx <- ifelse(get_exprs(c1_subset["ENSMUSG00000015579", isl1_cells], "norm_exprs_sf") > 0.2, "Nkx2-5+", "Nkx2-5-")
nkx_isl1 <- ifelse(get_exprs(c1_subset["ENSMUSG00000042258", nkx_cells], "norm_exprs_sf") > 0.2, "Isl1+", "Isl1-")
nkx_nkx <- ifelse(get_exprs(c1_subset["ENSMUSG00000015579", nkx_cells], "norm_exprs_sf") > 0.2, "Nkx2-5+", "Nkx2-5-")

co_expression_isl <- which(isl1_isl1 == "Isl1+" & isl1_nkx == "Nkx2-5+")
co_expression_nkx <- which(nkx_isl1 == "Isl1+" & nkx_nkx == "Nkx2-5+")

mix_color_grey <- rep("darkgrey", length(mix_color))
mix_color_grey[co_expression_isl] <- mix_color[co_expression_isl]
mix_color_grey[length(isl1_cells)+co_expression_nkx] <- mix_color[length(isl1_cells)+co_expression_nkx]

plot_mix_clustering_isl1 <- lattice::xyplot(pData(c1_subset)$tsne2[isl1_cells] ~ pData(c1_subset)$tsne1[isl1_cells],
                                             col = mix_color_grey[1:length(isl1_cells)], cex = .5, pch = 16, ylab = list(label = "t-SNE 2", cex = 0.6), main = list(label = "Isl1 cell clustering", cex = 0.75), scales = list(cex = 0.4, tck = c(0.5, 0)), xlab = list(label = "t-SNE 1", cex = 0.6), par.settings = global.settings)
plot_mix_clustering_nkx <- lattice::xyplot(pData(c1_subset)$tsne2[nkx_cells] ~ pData(c1_subset)$tsne1[nkx_cells],
                                            col = mix_color_grey[(length(isl1_cells)+1):length(nkx_cells)], cex = .5, pch = 16, ylab = list(label = "t-SNE 2", cex = 0.6), main = list(label = "Nkx2-5 cell clustering", cex = 0.75), scales = list(cex = 0.4, tck = c(0.5, 0)), xlab = list(label = "t-SNE 1", cex = 0.6), par.settings = global.settings)

#
# Table Co-expressing
#
nkx_cells_early <- which(pData(c1_subset)$Background == "nkx2-5" & pData(c1_subset)$Platform == "C1" & pData(c1_subset)$Timepoint == "e7.5")
isl1_cells_early <- which(pData(c1_subset)$Background == "isl1" & pData(c1_subset)$Platform == "C1" & pData(c1_subset)$Timepoint == "e7.5")
nkx_cells_middle <- which(pData(c1_subset)$Background == "nkx2-5" & pData(c1_subset)$Platform == "C1" & pData(c1_subset)$Timepoint == "e8.5")
isl1_cells_middle <- which(pData(c1_subset)$Background == "isl1" & pData(c1_subset)$Platform == "C1" & pData(c1_subset)$Timepoint == "e8.5")
nkx_cells_late <- which(pData(c1_subset)$Background == "nkx2-5" & pData(c1_subset)$Platform == "C1" & pData(c1_subset)$Timepoint == "e9.5")
isl1_cells_late <- which(pData(c1_subset)$Background == "isl1" & pData(c1_subset)$Platform == "C1" & pData(c1_subset)$Timepoint == "e9.5")

isl1_isl1_early <- ifelse(get_exprs(c1_subset["ENSMUSG00000042258", isl1_cells_early], "norm_exprs_sf") > 0.1, "Isl1+", "Isl1-")
isl1_nkx_early <- ifelse(get_exprs(c1_subset["ENSMUSG00000015579", isl1_cells_early], "norm_exprs_sf") > 0.1, "Nkx2-5+", "Nkx2-5-")
nkx_isl1_early <- ifelse(get_exprs(c1_subset["ENSMUSG00000042258", nkx_cells_early], "norm_exprs_sf") > 0.1, "Isl1+", "Isl1-")
nkx_nkx_early <- ifelse(get_exprs(c1_subset["ENSMUSG00000015579", nkx_cells_early], "norm_exprs_sf") > 0.1, "Nkx2-5+", "Nkx2-5-")

isl1_isl1_middle <- ifelse(get_exprs(c1_subset["ENSMUSG00000042258", isl1_cells_middle], "norm_exprs_sf") > 0.1, "Isl1+", "Isl1-")
isl1_nkx_middle <- ifelse(get_exprs(c1_subset["ENSMUSG00000015579", isl1_cells_middle], "norm_exprs_sf") > 0.1, "Nkx2-5+", "Nkx2-5-")
nkx_isl1_middle <- ifelse(get_exprs(c1_subset["ENSMUSG00000042258", nkx_cells_middle], "norm_exprs_sf") > 0.1, "Isl1+", "Isl1-")
nkx_nkx_middle <- ifelse(get_exprs(c1_subset["ENSMUSG00000015579", nkx_cells_middle], "norm_exprs_sf") > 0.1, "Nkx2-5+", "Nkx2-5-")

isl1_isl1_late <- ifelse(get_exprs(c1_subset["ENSMUSG00000042258", isl1_cells_late], "norm_exprs_sf") > 0.1, "Isl1+", "Isl1-")
isl1_nkx_late <- ifelse(get_exprs(c1_subset["ENSMUSG00000015579", isl1_cells_late], "norm_exprs_sf") > 0.1, "Nkx2-5+", "Nkx2-5-")
nkx_isl1_late <- ifelse(get_exprs(c1_subset["ENSMUSG00000042258", nkx_cells_late], "norm_exprs_sf") > 0.1, "Isl1+", "Isl1-")
nkx_nkx_late <- ifelse(get_exprs(c1_subset["ENSMUSG00000015579", nkx_cells_late], "norm_exprs_sf") > 0.1, "Nkx2-5+", "Nkx2-5-")


double_negative_isl1 <- c(isl1_cells_early[which(isl1_isl1_early == "Isl1-" & isl1_nkx_early == "Nkx2-5-")],
                          isl1_cells_middle[which(isl1_isl1_middle == "Isl1-" & isl1_nkx_middle == "Nkx2-5-")],
                          isl1_cells_late[which(isl1_isl1_late == "Isl1-" & isl1_nkx_late == "Nkx2-5-")])
double_negative_nkx <- c(nkx_cells_early[which(nkx_isl1_early == "Isl1-" & nkx_nkx_early == "Nkx2-5-")],
                         nkx_cells_middle[which(nkx_isl1_middle == "Isl1-" & nkx_nkx_middle == "Nkx2-5-")],
                         nkx_cells_late[which(nkx_isl1_late == "Isl1-" & nkx_nkx_late == "Nkx2-5-")])
non_double_negatives_isl1 <- c(isl1_cells_early[-which(isl1_isl1_early == "Isl1-" & isl1_nkx_early == "Nkx2-5-")],
                           isl1_cells_middle[-which(isl1_isl1_middle == "Isl1-" & isl1_nkx_middle == "Nkx2-5-")],
                           isl1_cells_late[-which(isl1_isl1_late == "Isl1-" & isl1_nkx_late == "Nkx2-5-")])
non_double_negatives_nkx <- c(nkx_cells_early[-which(nkx_isl1_early == "Isl1-" & nkx_nkx_early == "Nkx2-5-")],
                          nkx_cells_middle[-which(nkx_isl1_middle == "Isl1-" & nkx_nkx_middle == "Nkx2-5-")],
                          nkx_cells_late[-which(nkx_isl1_late == "Isl1-" & nkx_nkx_late == "Nkx2-5-")])


tt1 <- gridExtra::ttheme_minimal(base_size = 8,
                                 padding = unit(c(2, 2), "mm"),
                                 core = list(fg_params = list(hjust = 1, x = 0.9)),
                                 rowhead = list(fg_params = list(fontface = 2L, hjust = 1, x = 0.9)),
                                 colhead = list(fg_params = list(fontface = 2L, hjust = 1, x = 0.9)))

table_isl1_early <- gridExtra::tableGrob(table(isl1_isl1_early, isl1_nkx_early), theme=tt1)
table_nkx_early <- gridExtra::tableGrob(table(nkx_isl1_early, nkx_nkx_early), theme=tt1, rows = NULL)
table_isl1_middle <- gridExtra::tableGrob(table(isl1_isl1_middle, isl1_nkx_middle), theme=tt1)
table_nkx_middle <- gridExtra::tableGrob(table(nkx_isl1_middle, nkx_nkx_middle), theme=tt1, rows = NULL)
table_isl1_late <- gridExtra::tableGrob(table(isl1_isl1_late, isl1_nkx_late), theme=tt1)
table_nkx_late <- gridExtra::tableGrob(table(nkx_isl1_late, nkx_nkx_late), theme=tt1, rows = NULL)

#
# Load SOM input genes and predict Nkx2-5 data onto Isl1-SOM
#
het_early <- readLines(file.path(parameters$general$path_supdata, "heterogeneity-isl1.early.txt"))
het_mid <- readLines(file.path(parameters$general$path_supdata, "heterogeneity-isl1.mid.txt"))
het_late <- readLines(file.path(parameters$general$path_supdata, "heterogeneity-isl1.late.txt"))
het_pool <- readLines(file.path(parameters$general$path_supdata, "heterogeneity-isl1.txt"))

isl1_het <- unique(c(het_early, het_mid, het_late, het_pool))

expression <- get_exprs(c1_subset[, c(isl1_cells, nkx_cells)], "norm_exprs_sf")
expression_scaled <- t(scale(t(expression[isl1_het, ])))

som <- singlecellutils::calcSOM(expression_scaled, train = 1:nrow(expression_scaled), num_epochs = 2000, seed = 1004)
som$codes <- som$codes[[1]]

# Clustering

hdbscan <- function(data, min_samples = 7L, min_cluster_size = 9L, outlier = 0, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  h <- reticulate::import("hdbscan")
  cl <- h$HDBSCAN(min_samples = min_samples, min_cluster_size = min_cluster_size)
  labels <- cl$fit_predict(data) + 1
  labels[labels == 0] <- outlier
  return(factor(labels))
}

seed <- 3465
set.seed(seed)
tsne <- Rtsne::Rtsne(t(som$codes), perplexity = 21, theta = 0.05, max_iter = 2000)

p_color_background <- droplevels(pData(c1_subset[, c(isl1_cells, nkx_cells)])$Background)
p_color_background[match(c(double_negative_isl1, double_negative_nkx), c(isl1_cells, nkx_cells))] <- NA
p_cluster_background <- singlecellutils::colorAMap(tsne$Y, colour_by = p_color_background, palette = c("#31a354", "#3182bd"), cex = rep(0.5, nrow(tsne$Y)), na.cex = 0.5, xlab = list(label = "t-SNE 1", cex = 0.6), ylab = list(label = "t-SNE 2", cex = 0.6), main = list(label = ""), scales = list(cex = 0.4, tck = c(0.5, 0)), par.settings = global.settings)

p_color_time <- droplevels(pData(c1_subset[, c(isl1_cells, nkx_cells)])$Timepoint)
p_color_time[match(c(double_negative_isl1, double_negative_nkx), c(isl1_cells, nkx_cells))] <- NA
p_cluster_time <- singlecellutils::colorAMap(tsne$Y, colour_by = p_color_time, palette = parameters$colors$nkx_timepoint_palette, cex = rep(0.5, nrow(tsne$Y)), na.cex = 0.5, xlab = list(label = "t-SNE 1", cex = 0.6), ylab = list(label = ""), main = list(label = ""), scales = list(cex = 0.4, tck = c(0.5, 0)), par.settings = global.settings)

isl1_clustering <- factor(c(pData(c1_subset[, isl1_cells])$cluster, rep(NA, length(nkx_cells))), levels = as.character(1:5))
p_cluster_isl1 <- singlecellutils::colorAMap(tsne$Y, colour_by = isl1_clustering, palette = parameters$colors$isl_cluster_palette, cex = rep(0.5, nrow(tsne$Y)), na.cex = 0.5, xlab = list(label = "t-SNE 1", cex = 0.6), ylab = list(label = ""), main = list(label = ""), scales = list(cex = 0.4, tck = c(0.5, 0)), par.settings = global.settings)

#
# Predict Nkx2-5 data onto Isl1 trajectory
#

# Temp
#nkx_cells <- non_double_negatives_nkx
#isl1_cells <- non_double_negatives_isl1
# End Temp

m <- match(colnames(isl1_dm@data_env$data), rownames(c1_subset))
nkx_data <- get_exprs(c1_subset[m, nkx_cells], "norm_exprs_sf")

nkx_dcs <- dm_predict(isl1_dm, t(nkx_data))
nkx_coords <- matrix(nkx_dcs@x, ncol = nkx_dcs@Dim[2])
isl1_dm_coords_nkx <- rbind(isl1_dm_coords, cbind(-nkx_coords[, 1], -nkx_coords[, 3], -nkx_coords[, 2]))

#
# Prepare plots
#
# Create col
isl1_col <- rep("darkgrey", length(c(isl1_cells, nkx_cells)))
isl1_col[1:length(isl1_cells)] <- mix_color[1:length(isl1_cells)]
nkx_col <- rep("darkgrey", length(c(isl1_cells, nkx_cells)))
nkx_col[(length(isl1_cells)+1):length(nkx_col)] <- mix_color[(length(isl1_cells)+1):length(nkx_col)]

# Create pch
isl1_pch <- rep(4, length(c(isl1_cells, nkx_cells)))
isl1_pch[1:length(isl1_cells)] <- 16
nkx_pch <- rep(4, length(c(isl1_cells, nkx_cells)))
nkx_pch[(length(isl1_cells)+1):length(nkx_pch)] <- 16

# Create cex
isl1_cex <- rep(0.5, length(c(isl1_cells, nkx_cells)))
isl1_cex[1:length(isl1_cells)] <- .5
nkx_cex <- rep(0.5, length(c(isl1_cells, nkx_cells)))
nkx_cex[(length(isl1_cells)+1):length(nkx_cex)] <- .5

# Create data
plot_isl1_data_x <- rev(isl1_dm_coords_nkx[,1])
plot_isl1_data_y <- rev(isl1_dm_coords_nkx[,2])
plot_nkx_data_x <- isl1_dm_coords_nkx[,1]
plot_nkx_data_y <- isl1_dm_coords_nkx[,2]

#
# Create figure
#
table <- gridExtra::grid.arrange(table_isl1_early, table_nkx_early,
                                 table_isl1_middle, table_nkx_middle,
                                 table_isl1_late, table_nkx_late, nrow = 3)
plot_lineage_overlap_isl1 <- lattice::xyplot(plot_isl1_data_y ~ plot_isl1_data_x, cex = rev(isl1_cex), col = rev(isl1_col), pch = rev(isl1_pch), ylab = list(label = "Dimension 3", cex = 0.6), main = list(label = "Isl1 cells in Isl1 lineage", cex = 0.75), scales = list(cex = 0.4, tck = c(0.5, 0)), xlab = list(label = "Dimension 1", cex = 0.6), par.settings = global.settings)
plot_lineage_overlap_nkx <- lattice::xyplot(isl1_dm_coords_nkx[,2] ~ isl1_dm_coords_nkx[,1], cex = nkx_cex, col = nkx_col, pch = nkx_pch, ylab = list(label = "Dimension 3", cex = 0.6), main = list(label = "Nkx2-5 cells in Isl1 lineage", cex = 0.75), scales = list(cex = 0.4, tck = c(0.5, 0)), xlab = list(label = "Dimension 1", cex = 0.6), par.settings = global.settings)
plot_clustering <- gridExtra::grid.arrange(p_cluster_background, p_cluster_time, p_cluster_isl1, nrow = 1)

figure <- multipanelfigure::multi_panel_figure(width = 205, height = 230, columns = 3, rows = 3)
figure <- multipanelfigure::fill_panel(figure, table, scaling = "fit", row = 1, column = 1)

figure <- multipanelfigure::fill_panel(figure, plot_mix_clustering_isl1, row = 1, column = 3)
figure <- multipanelfigure::fill_panel(figure, plot_mix_clustering_nkx, row = 1, column = 2)

figure <- multipanelfigure::fill_panel(figure, plot_clustering, row = 2, column = 1:3)

figure <- multipanelfigure::fill_panel(figure, plot_lineage_overlap_isl1, row = 3, column = 1)
figure <- multipanelfigure::fill_panel(figure, plot_lineage_overlap_nkx, row = 3, column = 2)

# Save figure
ggplot2::ggsave(figure,
                file = "supplementary_figures/99_review_2c.pdf",
                width = 205,
                height = 230,
                units = "mm",
                useDingbats = FALSE)
