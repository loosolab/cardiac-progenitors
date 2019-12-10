library(scran)
library(scater)
library(singlecellutils)
library(RColorBrewer)
library(viridis)
library(dplyr)
library(tidyr)
library(rioja)
library(multipanelfigure)
library(ComplexHeatmap)
library(circlize)
library(zoo)

#
# Load data
#
source("src/parameters.R")
load(file.path(parameters$general$path_supdata, "c1_subset.Rdata"))

theme2 <- theme(plot.background = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                panel.border = element_blank(), panel.background = element_blank(), axis.line.x = element_line(size=.3), axis.line.y = element_line(size=.3), legend.text = element_text(size=6),
                axis.title.x = element_blank(), axis.title.y = element_text(color="black", size=10), plot.title = element_text(face="bold", color="black", size=6), legend.key.size =  unit(2, "mm"), legend.margin=unit(-25, "mm"),
                strip.background = element_blank(), strip.text = element_blank(), axis.text = element_text(color="black", size=8) )

ggMMplot <- function(var1, var2){
  require(ggplot2)
  levVar1 <- length(levels(var1))
  levVar2 <- length(levels(var2))

  jointTable <- prop.table(table(var1, var2))
  plotData <- as.data.frame(jointTable)
  plotData$marginVar1 <- prop.table(table(var1))
  plotData$var2Height <- plotData$Freq / plotData$marginVar1
  plotData$var1Center <- c(0, cumsum(plotData$marginVar1)[1:levVar1 -1]) +
    plotData$marginVar1 / 2

  ggplot(plotData, aes(var1Center, var2Height)) +
    geom_bar(stat = "identity", aes(fill = var2, width = marginVar1), col = "White") +
    geom_label(aes(label = as.character(var1), x = var1Center, y = 1.05), label.r = unit(0, "lines"))
}

#
# Isl1 KO Heatmap Top 30 genes
#
n_top_genes <- 30
isl1_wt_cells <- which(pData(c1_subset)$Background == "isl1" & (pData(c1_subset)$cluster == "2" | pData(c1_subset)$cluster == "1" | pData(c1_subset)$cluster == "5"))
isl1_ko_cells <- which(pData(c1_subset)$Background == "isl1ko")
#isl1_ki67 <- get_exprs(c1_subset, "norm_exprs_sf")["ENSMUSG00000031004", c(isl1_wt_cells, isl1_ko_cells)]

#
# Cycling score
#
load(file.path(parameters$general$path_rextdata, "c2.cpg.whitfield.v6.0.mmusculus.Rdata"))
names(c2_cpg_whitfield) <- gsub("WHITFIELD_CELL_CYCLE_", "", names(c2_cpg_whitfield))

# Identify genes from the sets that correlate with the score
cycle_score_genes <- lapply(c2_cpg_whitfield, function(gs) {
  m <- match(gs, rownames(c1_subset))
  e <- get_exprs(c1_subset[na.omit(m), c(isl1_wt_cells, isl1_ko_cells)], "exprs")
  score <- colMeans(e)
  cor <- cor(t(e), score, method = "spearman")
  u <- which(cor[,1] > 0.4)
  return(names(u))
})

# Calculate scores from genes
cycle_score_l <- lapply(cycle_score_genes, function(g) {
  m <- match(g, rownames(c1_subset))
  e <- get_exprs(c1_subset[na.omit(m), c(isl1_wt_cells, isl1_ko_cells)], "exprs")
  colMeans(e)
})
cycle_score <- as.data.frame(t(scale(t(do.call("rbind", cycle_score_l)))))
cycle_score$set <- rownames(cycle_score)

cycle_score_df <- spread(gather(as.data.frame(cycle_score), key = "cell", value = "score", -set), key = "set", value = "score")
m <- match(cycle_score_df$cell, colnames(c1_subset))
cycle_score_df$Background <- pData(c1_subset)$Background[m]

cycle_score_df %>%
  group_by(cell) %>%
  mutate(max_score = max(G1_S, G2_M)) %>%
  ungroup() %>%
  mutate(rank = dense_rank(desc(max_score)),
         is_cycling = max_score > -0.1) -> cycle_score_df

isl1_ko_cycleplot <- ggplot(cycle_score_df, aes(x = G1_S, y = G2_M)) +
  geom_point(aes(color = is_cycling, shape = Background)) +
  scale_color_manual(values = c("red", "black"), labels =  c("non cycling", "cycling")) +
  scale_shape_manual(values = 16:17, labels = c("Isl1 WT", "Isl1 KO")) +
  xlab("G1/S score") +
  ylab("G2/M score") +
  guides(color = guide_legend(title = "", label = T), shape = guide_legend(title = "")) +
  theme_classic() +
  theme(axis.title = element_text(size = 7), axis.text = element_text(size = 5))

isl1_ko_cyclerankplot <- ggplot(cycle_score_df, aes(x = rank, y = max_score)) +
  geom_point(aes(color = is_cycling, shape = Background)) +
  scale_color_manual(values = c("red", "black")) +
  xlab("Rank") +
  ylab("Maximal score (G1/S, G2/M)") +
  guides(color = FALSE, shape = FALSE) +
  theme_classic() +
  theme(axis.title = element_text(size = 7), axis.text = element_text(size = 5))

g1_s_g2_m_genes <- unlist(cycle_score_genes[c("G1_S", "G2_M")])
g1_s_g2_m_lengths <- unlist(lapply(cycle_score_genes[c("G1_S", "G2_M")], length))
names(g1_s_g2_m_genes) <- rep(names(g1_s_g2_m_lengths), times = g1_s_g2_m_lengths)

g1_s_g2_m_cor <- cor(t(get_exprs(c1_subset[g1_s_g2_m_genes, cycle_score_df$cell[cycle_score_df$is_cycling]], "exprs")))
colnames(g1_s_g2_m_cor) <- fData(c1_subset[colnames(g1_s_g2_m_cor),])$symbol
rownames(g1_s_g2_m_cor) <- fData(c1_subset[rownames(g1_s_g2_m_cor),])$symbol

write.table(file = file.path(parameters$general$path_supdata, "Source_Data_SupplementaryFigure_9c.txt"),
            x = g1_s_g2_m_cor,
            row.names = T,
            col.names = T,
            quote = F,
            sep = "\t")

col_fun <- colorRamp2(seq(-0.9, 0.9, length.out = 10), rev(brewer.pal(10, "RdBu")))

isl1_ko_cycle_genes <- Heatmap(g1_s_g2_m_cor,
                               col = col_fun,
                               name = "Pearson\ncorrelation",
                               split = names(g1_s_g2_m_genes),
                               cluster_columns = F,
                               cluster_rows = F,
                               row_names_gp = gpar(fontsize = 6),
                               column_names_gp = gpar(fontsize = 6))

t <- table(isko = droplevels(cycle_score_df$Background) == "isl1ko", is_cycling = cycle_score_df$is_cycling)[c(2,1), c(2,1)]
tr <- prop.test(t, correct=F, conf.level = 0.95)

isl1_ko_mmplot <- ggMMplot(droplevels(cycle_score_df$Background), cycle_score_df$is_cycling) +
  guides(fill = F) +
  scale_fill_manual(values = c("#A8DBA8", "#3B8686")) +
  ylab("Proportion of cycling cells") +
  xlab("Cells") +
  theme(panel.background = element_rect(fill = "transparent"), axis.line.y = element_line(colour = "black"), axis.line.x = element_line(colour = "black"))#, axis.text.x = element_blank(), axis.title.x = element_blank(), axis.ticks.x = element_blank())

m <- multi_panel_figure(width = c(85,110), height = c(80, 150), label_type = "lower_alpha")
m <- fill_panel(m, isl1_ko_cyclerankplot)
m <- fill_panel(m, isl1_ko_cycleplot)
m <- fill_panel(m, isl1_ko_cycle_genes, column = 1:2)

ggsave(
  plot = m,
  filename = file.path(parameters$general$path_rfigures, "Supplementary_Figure_9.pdf"),
  width = figure_width(m),
  height = figure_height(m),
  units = "mm",
  dpi = 600,
  useDingbats = FALSE)

#
# Find top 30 genes per KO-cluster comparison
#
comparisons <- c("progenitor", "endothelial", "cardiomyocyte")

isl1_ko_data_l <- lapply(comparisons, function(x) {
  d <- read.table(file.path(parameters$general$path_supdata, paste0("differentialExpression/Isl1.", x, "_Isl1KO.txt")), sep="\t", header = T, stringsAsFactors = F)
  d$contrast <- x
  i <- which(d$fdr < parameters$diffusionmaps$de_fdr & (d$lfc.lo > 1.5 | d$lfc.hi < -1.5) & d$biotype == "protein_coding" & !startsWith(d$symbol, "mt") & !startsWith(d$symbol, "Rpl") & !startsWith(d$symbol, "venus"))
  m <- order(abs(d$lfc.hi[i]), decreasing = T)
  na.omit(d$geneID[i][m[1:min(n_top_genes, length(m))]])
})

isl1_de_selected <- na.omit(unique(unlist(isl1_ko_data_l)))
isl1_subset <- c1_subset[isl1_de_selected, c(isl1_wt_cells, isl1_ko_cells)]

#
# Isl1 KO heatmap
#
isl1_ko_cell_order <- order(pData(isl1_subset)$cluster)

isl1_ko_expression <- t(scale(t(get_exprs(isl1_subset, "norm_exprs_sf")[, isl1_ko_cell_order])))
rownames(isl1_ko_expression) <- fData(isl1_subset)$symbol

isl1_ko_dist <- as.dist(1 - cor(isl1_ko_expression, method="pearson"))

isl1_ko_const_clust <- chclust(isl1_ko_dist)
isl1_ko_const_clust_dendro <- as.dendrogram(isl1_ko_const_clust)

# Color mapping
e_min <- min(isl1_ko_expression)
e_max <- max(isl1_ko_expression)

bounds <- max(abs(e_min), abs(e_max))

breaks <- seq(-bounds, bounds, length.out = 99)
#breaks <- seq(e_min, e_max, length.out = 99)
colors <- colorRampPalette(rev(brewer.pal(10, "RdBu")))(99)
#colors <- inferno(99)

col_fun <- colorRamp2(breaks, colors)

#
# Heatmap Annotation
#
isl1_ko_annotation_data <- data.frame(Genotype = paste(pData(isl1_subset)$Background))

isl1_ko_annotation <- HeatmapAnnotation(isl1_ko_annotation_data,
                                        col = list(Genotype = c("isl1" = "#000000", "isl1ko" = "#FF0000")), #008B45, 00CD66
                                        show_legend = F,
                                        show_annotation_name = T,
                                        annotation_name_gp = gpar(fontsize = 6),
                                        annotation_name_offset = unit(0.8, "mm"),
                                        gap = unit(0.1, "mm"),
                                        annotation_height = unit.c(unit(1.5, "mm")))

isl1_ko_heatmap <- Heatmap(isl1_ko_expression,
                          col = col_fun,
                          top_annotation = isl1_ko_annotation,
                          cluster_columns = isl1_ko_const_clust_dendro,
                          cluster_rows = T,
                          show_column_names = F,
                          show_row_dend = F,
                          show_column_dend = F,
                          name = "islko",
                          #column_title = "Nkx2-5+ cells",
                          column_title_gp = gpar(fontface = "bold", fontsize = 8),
                          #split = rep(paste("Marker genes\ncluster", c(1,2,3)), times = sapply(nkx_de_gene_order, length)),
                          gap = unit(0.3, "mm"),
                          row_title_rot = 90,
                          row_title_gp = gpar(fontsize = 4),
                          row_names_gp = gpar(fontsize = 6),
                          heatmap_legend_param = list("plot" = F, "color_bar" = "continuous", "legend_direction" = "horizontal", "title" = "Scaled expression", "title_position" = "leftcenter", "title_gp" = gpar(fontsize = 4), "labels_gp" = gpar(fontsize = 4), "grid_height" = unit(2, "mm")))

fig_3c <- grid.grabExpr(draw(isl1_ko_heatmap, heatmap_legend_side = "bottom", padding = unit(c(0, 0, 0, 0), "mm")))

# Write out for Source Data file
write.table(file = file.path(parameters$general$path_supdata, "Source_Data_Figure5_d.txt"),
            x = isl1_ko_expression,
            col.names = T,
            row.names = T,
            sep = "\t",
            quote = T)

#
# Predict Isl1-KO onto Isl1 trajectory
#
library(destiny)
global.settings <- list(layout.heights=list(main.key.padding = -2, bottom.padding = 0, axis.xlab.padding = -0.5), layout.widths = list(ylab.axis.padding = -0.5, right.padding = 0))

load(file.path(parameters$general$path_supdata, "Isl1-DM.Rdata"))
isl1_cells <- which(pData(c1_subset)$Background == "isl1" & pData(c1_subset)$Platform == "C1")
isl1_dm_coords <- cbind(-pData(c1_subset[, isl1_cells])$dm1, -pData(c1_subset[, isl1_cells])$dm2)

m <- match(colnames(isl1_dm@data_env$data), rownames(c1_subset))
isl1ko_cells <- which(pData(c1_subset)$Background == "isl1ko")
isl1ko_data <- get_exprs(c1_subset[na.omit(m), isl1ko_cells], "norm_exprs_sf")

isl1ko_dcs <- dm_predict(isl1_dm, t(isl1ko_data))
isl1ko_coords <- matrix(isl1ko_dcs@x, ncol = isl1ko_dcs@Dim[2])
isl1_dm_coords_ko <- rbind(isl1_dm_coords, cbind(-isl1ko_coords[, 1], -isl1ko_coords[, 3]))

cluster <- c(pData(c1_subset[, isl1_cells])$cluster, rep("KO", nrow(isl1ko_coords)))
rm_obs <- is.na(cluster)
cluster <- paste0("cluster", cluster[!rm_obs])

p_isl1_ko <- colorAMap(isl1_dm_coords_ko[!rm_obs,], colour_by = factor(cluster), palette = parameters$colors$islko_cluster_palette, ylab = list(label = "Dimension 3", cex = 0.6), main = list(label = "", cex = 0.75), scales = list(cex = 0.4, tck = c(0.5, 0)), xlab = list(label = "Dimension 1", cex = 0.6), par.settings = global.settings)

#
# Assemble figure
#
m <- multi_panel_figure(width = c(69, 69, 64, 5), height = c(65, 165), column_spacing = c(6,0,0,0), panel_label_type = "lower-alpha")
m <- fill_panel(m, "ext_data/Isl1_KO.png", label = "a", row = 1, column = 1)
m <- fill_panel(m, p_isl1_ko, row = 1, column = 2, label = "b")
m <- fill_panel(m, isl1_ko_mmplot, label = "c", row = 1, column = 3:4)
m <- fill_panel(m, fig_3c, label = "d", column = 1:3, row = 2)

ggsave(
  plot = m,
  filename = file.path(parameters$general$path_rfigures, "Figure_3.pdf"),
  width = figure_width(m),
  height = figure_height(m),
  units = "mm",
  dpi = 600,
  useDingbats = FALSE)

