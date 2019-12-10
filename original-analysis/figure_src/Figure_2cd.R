library(scater)
library(dplyr)
library(tidyr)
library(singlecellutils)
library(multipanelfigure)
library(gridExtra)

#
# Load data
#
source("src/parameters.R")
load(file.path(parameters$general$path_supdata, "c1_subset.Rdata"))

theme2 <- theme(plot.background = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                panel.border = element_blank(), panel.background = element_blank(), axis.line.x = element_line(size=.3), axis.line.y = element_line(size=.3), legend.text = element_text(size=6),
                axis.title.x = element_blank(), axis.title.y = element_text(color="black", size=10), plot.title = element_text(face="bold", color="black", size=6), legend.key.size =  unit(2, "mm"), legend.margin=unit(-25, "mm"),
                strip.background = element_blank(), strip.text = element_blank(), axis.text = element_text(color="black", size=8), axis.text.x = element_blank())

theme_violins <- theme(plot.background = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                     panel.border = element_blank(), panel.background = element_blank(), axis.line.y = element_line(size=.3), axis.line.x = element_line(size=.3), legend.text = element_text(size=6),
                     axis.title.y = element_text(color="black", size=10), axis.title.x = element_blank(), plot.title = element_text(face="bold", color="black", size=6), legend.key.size =  unit(2, "mm"), legend.margin=unit(-25, "mm"),
                     strip.background = element_blank(), strip.text = element_blank(), axis.text = element_text(color="black", size=8), axis.text.x = element_blank())

theme_hists <- theme(plot.background = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                     panel.border = element_blank(), panel.background = element_blank(), axis.line.x = element_line(size=.3), axis.line.y = element_blank(), legend.text = element_text(size=6),
                     axis.title.x = element_text(color="black", size=10), axis.title.y = element_blank(), plot.title = element_text(face="bold", color="black", size=6), legend.key.size =  unit(2, "mm"), legend.margin=unit(-25, "mm"),
                     strip.background = element_blank(), strip.text = element_blank(), axis.text = element_text(color="black", size=8), axis.ticks.y = element_blank(), axis.text.y = element_blank() )

theme_hists_2 <- theme(plot.background = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                     panel.border = element_blank(), panel.background = element_blank(), axis.line.x = element_line(size=.3), axis.line.y = element_blank(), legend.text = element_text(size=6),
                     axis.title.x = element_text(color="black", size=10), axis.title.y = element_blank(), plot.title = element_text(face="bold", color="black", size=6), legend.key.size =  unit(2, "mm"), legend.margin=unit(-25, "mm"),
                     strip.background = element_blank(), strip.text = element_blank(), axis.text = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank() )

#
# Nkx2-5 IC analysis
#
nkx_cells <- which(pData(c1_subset)$Background == "nkx2-5" & pData(c1_subset)$Platform == "C1")
nkx_de_genes <- rownames(fData(c1_subset))[which(fData(c1_subset)$nkx_de & fData(c1_subset)$nkx_marker)]
nkx_cluster <- factor(pData(c1_subset[, nkx_cells])$cluster)

nkx_data <- get_exprs(c1_subset[nkx_de_genes, nkx_cells], "norm_exprs_sf")

nkx_ic_df <- as.data.frame(boot.ic(nkx_data, groups = nkx_cluster, R = 1000, n = 30, p.val = 0.05))
colnames(nkx_ic_df) <- paste("cluster", levels(nkx_cluster))

# Create plots
nkx_ic <- gather(nkx_ic_df, key = "cluster", value = "ic")
nkx_ic$cluster <- factor(nkx_ic$cluster)

nkx_ic_plot <- ggplot(nkx_ic, aes(x = cluster, y = ic, fill = cluster)) +
  geom_boxplot() +
  scale_fill_manual(values = parameters$colors$nkx_cluster_palette) +
  guides(fill = F) +
  ylab("Critical transition index") +
  xlab("") +
  theme2 +
  theme(plot.margin = unit(c(0, -5, 2.5, 0), "mm"))

#
# Nkx pairwise statistical tests
#
nkx_ptest_list <- apply(combn(levels(nkx_cluster), 2), 2, function(x) {
  x <- as.numeric(x)
  ks <- ks.test(nkx_ic_df[, x[1]], nkx_ic_df[, x[2]])
  wr <- wilcox.test(nkx_ic_df[, x[1]], nkx_ic_df[, x[2]])

  data.frame(combination = paste(x[1], "vs.", x[2]),
       ks_p_value = ks$p.value,
       ks_statistic = ks$statistic,
       wilcox_p_value = wr$p.value,
       wilcox_statistic = wr$statistic)
})
nkx_ptest <- data.frame(do.call("rbind", nkx_ptest_list))

#
# Calculate transcriptome variance
#
nkx_colors <- parameters$colors$nkx_cluster_palette
names(nkx_colors) <- levels(nkx_cluster)

nkx_data <- get_exprs(c1_subset[nkx_de_genes, nkx_cells], "norm_exprs_sf")

nkx_cluster_noise_list <- lapply(levels(nkx_cluster), function(c) {
  cells <- which(nkx_cluster == c)
  distances <- apply(combn(cells, 2), 2, function(d) {
    p <- cor(nkx_data[, d[1]], nkx_data[, d[2]], method = "spearman")
    sqrt((1-p)/2)
  })
})

nkx_cluster_noise <- data.frame(cluster = rep(levels(nkx_cluster), times = sapply(nkx_cluster_noise_list, length)), distance = unlist(nkx_cluster_noise_list))
p_nkx_cluster_noise <- ggplot(data = nkx_cluster_noise, aes(x = cluster, y = distance, fill = cluster)) +
  geom_violin(trim = F) +
  geom_boxplot(width = 0.1, fill = "white") +
  scale_fill_manual(values = nkx_colors) +
  scale_y_continuous(expand = c(0.05, 0)) +
  guides(fill = FALSE, colour = FALSE) +
  xlab("") + ylab("Cell to cell distance") + theme_violins +
  theme(plot.margin = unit(c(0, -5, 2.5, 0), "mm"))

#
# Nkx2-5 CxC and GxG correlation analysis
#
nkx_colors <- parameters$colors$nkx_cluster_palette
names(nkx_colors) <- levels(nkx_cluster)

nkx_correlation_plot_list <- lapply(levels(nkx_cluster), function(c) {
  cells <- which(nkx_cluster == c)
  ic <- computeIC(nkx_data[, cells], p.val = 1)

  p1 <- ggplot(data.frame(t=ic$between), aes(x=t)) +
    geom_histogram(binwidth = 0.01, fill = nkx_colors[c]) +
    scale_x_continuous(expand = c(0,0), limits = c(0,1)) +
    scale_y_continuous(expand = c(0,0.3))

  p2 <- ggplot(data.frame(t=ic$within), aes(x=t)) +
    geom_histogram(binwidth = 0.01, fill = nkx_colors[c]) +
    scale_x_continuous(expand = c(0,0), limits = c(0,1)) +
    scale_y_continuous(expand = c(0,0.3))

  if( c == levels(nkx_cluster)[length(levels(nkx_cluster))] ) {
    p1 <- p1 + xlab("Pearson's correlation") + theme_hists
    p2 <- p2 + xlab("Pearson's correlation") + theme_hists
  } else {
    p1 <- p1 + xlab("") + theme_hists_2
    p2 <- p2 + xlab("") + theme_hists_2
  }
  p <- grid.arrange(p2, p1, ncol = 2)
})

#
# Isl1 IC analysis
#
isl1_cells <- which(pData(c1_subset)$Background == "isl1" & pData(c1_subset)$Platform == "C1")
isl1_de_genes <- rownames(fData(c1_subset))[which(fData(c1_subset)$isl1_de & fData(c1_subset)$isl1_marker)]
isl1_cluster <- factor(pData(c1_subset[, isl1_cells])$cluster)

isl1_ko_cells <- which((pData(c1_subset)$Background == "isl1" | pData(c1_subset)$Background == "isl1ko") & pData(c1_subset)$Platform == "C1")
isl1_ko_cluster <- pData(c1_subset[, isl1_ko_cells])$cluster
isl1_ko_cluster[which((pData(c1_subset[, isl1_ko_cells])$Background == "isl1ko") & pData(c1_subset[, isl1_ko_cells])$Platform == "C1")] <- "KO"
isl1_ko_cluster <- factor(isl1_ko_cluster)

isl1_data <- get_exprs(c1_subset[isl1_de_genes, isl1_cells], "norm_exprs_sf")
isl1_ko_data <- get_exprs(c1_subset[isl1_de_genes, isl1_ko_cells], "norm_exprs_sf")

isl1_ic_df <- as.data.frame(boot.ic(isl1_data, groups = isl1_cluster, R = 1000, n = 20, p.val = 0.05))
colnames(isl1_ic_df) <- paste("cluster", levels(isl1_cluster))

isl1_ko_ic_df <- as.data.frame(boot.ic(isl1_ko_data, groups = isl1_ko_cluster, R = 1000, n = 20, p.val = 0.05))
colnames(isl1_ko_ic_df) <- paste("cluster", levels(isl1_ko_cluster))

# Create plots
isl1_ic <- gather(isl1_ic_df, key = "cluster", value = "ic")
isl1_ko_ic <- gather(isl1_ko_ic_df, key = "cluster", value = "ic")

isl1_ic_plot <- ggplot(isl1_ic, aes(x = cluster, y = ic, fill = cluster)) +
  geom_boxplot() +
  scale_fill_manual(values = parameters$colors$isl_cluster_palette) +
  guides(fill = F) +
  ylab("") +
  xlab("") +
  theme2

isl1_ko_ic_plot <- ggplot(isl1_ko_ic, aes(x = cluster, y = ic, fill = cluster)) +
  geom_boxplot() +
  scale_fill_manual(values = parameters$colors$isl_cluster_palette) +
  guides(fill = F) +
  ylab("") +
  xlab("") +
  theme2

#
# Isl1 pairwise statistical tests
#
isl1_ptest_list <- apply(combn(levels(isl1_cluster), 2), 2, function(x) {
  x <- as.numeric(x)
  ks <- ks.test(isl1_ic_df[, x[1]], isl1_ic_df[, x[2]])
  wr <- wilcox.test(isl1_ic_df[, x[1]], isl1_ic_df[, x[2]])

  data.frame(combination = paste(x[1], "vs.", x[2]),
       ks.p_value = ks$p.value,
       ks.statistic = ks$statistic,
       wilcox.p_value = wr$p.value,
       wilcox.statistic = wr$statistic)
})
isl1_ptest <- do.call("rbind", isl1_ptest_list)

#
# Calculate transcriptome variance
#
isl1_colors <- parameters$colors$isl_cluster_palette
names(isl1_colors) <- levels(isl1_cluster)

isl1_cluster_noise_list <- lapply(levels(isl1_cluster), function(c) {
  cells <- which(isl1_cluster == c)
  distances <- apply(combn(cells, 2), 2, function(d) {
    p <- cor(isl1_data[, d[1]], isl1_data[, d[2]], method = "spearman")
    sqrt((1-p)/2)
  })
})

isl1_cluster_noise <- data.frame(cluster = rep(levels(isl1_cluster), times = sapply(isl1_cluster_noise_list, length)), distance = unlist(isl1_cluster_noise_list))
p_isl1_cluster_noise <- ggplot(data = isl1_cluster_noise, aes(x = cluster, y = distance, fill = cluster)) +
  geom_violin(trim = F) +
  geom_boxplot(width = 0.1, fill = "white") +
  scale_fill_manual(values = isl1_colors) +
  scale_y_continuous(expand = c(0.05, 0)) +
  guides(fill = FALSE, colour = FALSE) +
  xlab("") + ylab("") + theme_violins

#
# Isl1 CxC and GxG correlation analysis
#
isl1_ko_cells <- which((pData(c1_subset)$Background == "isl1" | pData(c1_subset)$Background == "isl1ko") & pData(c1_subset)$Platform == "C1")
isl1_ko_cluster <- pData(c1_subset[, isl1_ko_cells])$cluster
isl1_ko_cluster[which((pData(c1_subset)$Background == "isl1ko") & pData(c1_subset)$Platform == "C1")] <- "KO"
isl1_ko_cluster <- factor(isl1_ko_cluster)

isl1_colors <- parameters$colors$isl_cluster_palette
names(isl1_colors) <- levels(isl1_ko_cluster)

isl1_ko_data <- get_exprs(c1_subset[isl1_de_genes, isl1_ko_cells], "norm_exprs_sf")

isl1_correlation_plot_list <- lapply(levels(isl1_ko_cluster), function(c) {
  cells <- which(isl1_ko_cluster == c)
  ic <- computeIC(isl1_ko_data[, cells], p.val = 1)

  p1 <- ggplot(data.frame(t=ic$between), aes(x=t)) +
    geom_histogram(binwidth = 0.01, fill = isl1_colors[c]) +
    scale_x_continuous(expand = c(0,0), limits = c(0,1)) +
    scale_y_continuous(expand = c(0,0.3))
  p2 <- ggplot(data.frame(t=ic$within), aes(x=t)) +
    geom_histogram(binwidth = 0.01, fill = isl1_colors[c]) +
    scale_x_continuous(expand = c(0,0), limits = c(0,1)) +
    scale_y_continuous(expand = c(0,0.3))

  if( c == levels(isl1_ko_cluster)[length(levels(isl1_ko_cluster))] ) {
    p1 <- p1 + xlab("Pearson's correlation") + theme_hists
    p2 <- p2 + xlab("Pearson's correlation") + theme_hists
  } else {
    p1 <- p1 + xlab("") + theme_hists_2
    p2 <- p2 + xlab("") + theme_hists_2
  }
  p <- grid.arrange(p2, p1, ncol = 2)
})

# Fig2 c,d
m <- multi_panel_figure(width = 205, height = 230, columns = 5, rows = 4, column_spacing = c(unit(5, "mm"), unit(5, "mm"), unit(5, "mm"), unit(0, "mm"), unit(0, "mm")))
m <- fill_panel(m, nkx_ic_plot, label = "C", column = 3, row = 1)
m <- fill_panel(m, isl1_ic_plot, label = "", column = 4:5, row = 1)
m <- fill_panel(m, p_nkx_cluster_noise, label = "D", column = 3, row = 2)
m <- fill_panel(m, p_isl1_cluster_noise, label = "", column = 4:5, row = 2)

ggsave(
  plot = m,
  filename = file.path(parameters$general$path_rfigures, "Figure_2_cd.pdf"),
  width = 205,
  height = 230,
  units = "mm",
  dpi = 600)

#
# Save raw data for Source File
#
write.table(file = file.path(parameters$general$path_supdata, "Source_Data_Figure2_c_Nkx.txt"),
            x = nkx_ic,
            row.names = F,
            col.names = T,
            quote = F,
            sep = "\t")

# Supplementary Figure 2c/d
nkx_correlation_plot <- grid.arrange(grobs = nkx_correlation_plot_list, ncol = 1)
isl1_correlation_plot <- grid.arrange(grobs = isl1_correlation_plot_list, ncol = 1)

m <- multi_panel_figure(width = 205, height = 230, columns = 2, rows = 7)
m <- fill_panel(m, nkx_correlation_plot, label = "C", column = 1, row = 5:7)
m <- fill_panel(m, isl1_correlation_plot, label = "D", column = 2, row = 3:7)

ggsave(
  plot = m,
  filename = file.path(parameters$general$path_rfigures, "Supplementary_Figure_2_IC.pdf"),
  width = 205,
  height = 230,
  units = "mm",
  dpi = 600)

# Supplementary Table with p-values
tt1 <- ttheme_minimal(core=list(fg_params=list(hjust=1, x=0.9)),
                      rowhead=list(fg_params=list(hjust=1, x=0.9)),
                      colhead = list(fg_params=list(hjust=1, x=0.9)))

s1_table_nkx <- tableGrob(nkx_ptest, theme=tt1, rows = NULL)
s1_table_isl1 <- tableGrob(isl1_ptest, theme=tt1, rows = NULL)

m <- multi_panel_figure(width = 205, height = 230, columns = 1, rows = 2)
m <- fill_panel(m, s1_table_nkx)
m <- fill_panel(m, s1_table_isl1)

ggsave(
  plot = m,
  filename = file.path(parameters$general$path_rfigures, "Supplementary_Table_IC-pvalues.pdf"),
  width = 205,
  height = 230,
  units = "mm",
  dpi = 600)
