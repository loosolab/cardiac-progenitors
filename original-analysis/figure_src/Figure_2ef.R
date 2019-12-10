library(scater)
library(tidyr)
library(ggplot2)

#
# Load data
#
source("src/parameters.R")
load(file.path(parameters$general$path_supdata, "c1_subset.Rdata"))

#
# Setup
#
# Genes for Isl1 dpt lineplots
isl1_genes <- factor(c("Isl1", "Mef2c", "Nkx2-5", "Tgfb1li1", "Ankrd1", "Myocd", "Smyd1", "Gata2", "Hand1", "Hhex", "Hoxa7", "Hoxa9", "Hoxb6", "Hoxc8", "Msx1", "Pitx1", "Nkx1-2", "Snai1", "Tbx3", "Tbx4", "Ets2", "Sall4"),
                     levels = c("Isl1", "Mef2c", "Nkx2-5", "Tgfb1li1", "Ankrd1", "Myocd", "Smyd1", "Gata2", "Hand1", "Hhex", "Hoxa7", "Hoxa9", "Hoxb6", "Hoxc8", "Msx1", "Pitx1", "Nkx1-2", "Snai1", "Tbx3", "Tbx4", "Ets2", "Sall4"),
                     ordered = T)

# Genes for Nkx2-5 dpt lineplots
nkx_genes <- factor(c("Nkx2-5", "Ankrd1", "Cdkn2d", "Hopx", "Mef2c", "Myocd", "Smyd1", "Tgfb1li1", "Tbx20", "Dnmt3b", "Gata2", "Gata3", "Hand1", "Msx1"),
                    levels = c("Nkx2-5", "Ankrd1", "Cdkn2d", "Hopx", "Mef2c", "Myocd", "Smyd1", "Tgfb1li1", "Tbx20", "Dnmt3b", "Gata2", "Gata3", "Hand1", "Msx1"),
                    ordered = T)

theme <- theme(panel.grid.major = element_blank(),
               panel.background = element_blank(),
               strip.background = element_blank(),
               strip.text = element_text(face = "bold", size = 5, margin = margin(0.5, 1, 0, 1, unit = "mm")),
               panel.border = element_rect(fill = NA, color = "black"),
               panel.spacing = grid::unit(0.5,"mm"),
               axis.text = element_text(size = 4),
               axis.title = element_text(size = 5),
               axis.title.x = element_text(margin = margin(0,0,0,0)),
               axis.title.y = element_text(margin = margin(0,0,0,0))
               )

#
# Generate Isl1 dpt line plots
#
isl1_cells <- which(pData(c1_subset)$Background == "isl1" & pData(c1_subset)$Platform == "C1")

m <- match(isl1_genes, fData(c1_subset)$symbol)
isl1_expression <- as.data.frame(get_exprs(c1_subset[na.omit(m), isl1_cells], "norm_exprs_sf"))
isl1_expression$symbol <- isl1_genes[!is.na(m)]

# Tidyr gather into tidy data
isl1_expression_df <- gather(isl1_expression, key = "cell", value = "expression", -symbol)

# Add cell metadata
m <- match(isl1_expression_df$cell, colnames(c1_subset[, isl1_cells]))
isl1_expression_df$dpt <- pData(c1_subset[, isl1_cells])$dpt[m]
isl1_expression_df$cluster <- factor(pData(c1_subset[, isl1_cells])$cluster[m])

# Exclude endothelial branch
isl1_expression_df <- subset(isl1_expression_df, isl1_expression_df$cluster != 1)

isl1_dpt_lineplot <- ggplot(isl1_expression_df, aes(x = dpt, y = expression)) +
  geom_point(data = subset(isl1_expression_df, !is.na(isl1_expression_df$cluster)), aes(color = cluster), size = 0.3) +
  geom_point(data = subset(isl1_expression_df, is.na(isl1_expression_df$cluster)), color = "grey", size = 0.3) +
  stat_smooth(span = 0.9, method = "loess", n = 30, color = "grey", se = FALSE) +
  scale_color_manual(values = parameters$colors$isl_cluster_palette[2:5]) +
  facet_wrap(~ symbol, ncol = 6) +
  #ylab("Normalized expression") +
  #xlab("Pseudotime") +
  ylab("") +
  xlab("") +
  guides(color = FALSE) +
  theme

#
# Generate Nkx2-5 dpt line plots
#
nkx_cells <- which(pData(c1_subset)$Background == "nkx2-5" & pData(c1_subset)$Platform == "C1")

m <- match(nkx_genes, fData(c1_subset)$symbol)
nkx_expression <- as.data.frame(get_exprs(c1_subset[na.omit(m), nkx_cells], "norm_exprs_sf"))
nkx_expression$symbol <- nkx_genes[!is.na(m)]

# Tidyr gather into tidy data
nkx_expression_df <- gather(nkx_expression, key = "cell", value = "expression", -symbol)

# Add cell metadata
m <- match(nkx_expression_df$cell, colnames(c1_subset[, nkx_cells]))
nkx_expression_df$dpt <- pData(c1_subset[, nkx_cells])$dpt[m]
nkx_expression_df$cluster <- factor(pData(c1_subset[, nkx_cells])$cluster[m])

nkx_dpt_lineplot <- ggplot(nkx_expression_df, aes(x = dpt, y = expression)) +
  geom_point(data = subset(nkx_expression_df, !is.na(nkx_expression_df$cluster)), aes(color = cluster), size = 0.3) +
  geom_point(data = subset(nkx_expression_df, is.na(nkx_expression_df$cluster)), color = "grey", size = 0.3) +
  stat_smooth(span = 0.9, method = "loess", n = 30, color = "grey", se = FALSE) +
  scale_color_manual(values = parameters$colors$nkx_cluster_palette) +
  facet_wrap(~ symbol, nrow = 3) +
  #ylab("Normalized expression") +
  #xlab("Pseudotime") +
  ylab("") +
  xlab("") +
  guides(color = FALSE) +
  theme

#
# Save plots (uncomment to activate)
#
figure <- multipanelfigure::multi_panel_figure(width = 205, height = 230, columns = 2, rows = 2)
figure <- multipanelfigure::fill_panel(figure, nkx_dpt_lineplot, column = 1, row = 2, label = "e")
figure <- multipanelfigure::fill_panel(figure, isl1_dpt_lineplot, column = 2, row = 2, label = "f")

ggsave(plot = figure,
       filename = file.path(parameters$general$path_rfigures, "Figure_2_ef.pdf"),
       width = 205,
       height = 230,
       units = "mm",
       dpi = 600)
# ggsave(plot = isl1_dpt_lineplot,
#        filename = "Figure_2_f.pdf",
#        width = 205,
#        height = 230,
#        units = "mm",
#        dpi = 600)

#
# Save raw tables for Source Data
#
write.table(file = file.path(parameters$general$path_supdata, "Source_Data_Figure_2e.txt"),
            x = isl1_expression_df,
            col.names = T,
            row.names = F,
            sep = "\t",
            quote = F)
write.table(file = file.path(parameters$general$path_supdata, "Source_Data_Figure_2f.txt"),
            x = nkx_expression_df,
            col.names = T,
            row.names = F,
            sep = "\t",
            quote = F)
