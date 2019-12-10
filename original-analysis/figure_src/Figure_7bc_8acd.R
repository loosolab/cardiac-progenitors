load(file = "data/scData_filtered_2.Rda")

#
# Figure 7 b/c
#

plot_lsi_clustering <- scater::plotReducedDim(scData_filtered, use_dimred = "tsne", colour_by = ".cluster_5", shape_by = "batch", add_ticks = F) +
  ggplot2::ggtitle("LSI + HDBSCAN clustering") +
  ggplot2::guides(shape = FALSE) +
  ggplot2::xlab("") + ggplot2::ylab("") +
  ggplot2::theme(plot.title = ggplot2::element_text(face="bold", color="black", size=6), legend.title = ggplot2::element_blank()) +
  ggplot2::scale_color_manual(values = c("1" = "#33a02c", "2" = "#6a3d9a", "3" = "#1f78b4", "4" = "#e31a1c", "5" = "#ff7f00")) +
  ggplot2::scale_shape_manual(values = c(16, 16))

plot_lsi_batch <- scater::plotReducedDim(scData_filtered, use_dimred = "tsne", colour_by = "batch", add_ticks = F) +
  ggplot2::ggtitle("LSI: Cell timepoint") +
  ggplot2::guides(shape = FALSE) +
  ggplot2::xlab("") + ggplot2::ylab("") +
  ggplot2::theme(plot.title = ggplot2::element_text(face="bold", color="black", size=6), legend.title = ggplot2::element_blank())

figure_7bc <- multipanelfigure::multi_panel_figure(rows = 4, columns = 2, width = 205, height = 230)

figure_7bc <- multipanelfigure::fill_panel(figure_7bc, plot_lsi_clustering, column = 1:2, row = 1:2)
figure_7bc <- multipanelfigure::fill_panel(figure_7bc, plot_lsi_batch, column = 1:2, row = 3:4)

#
# Figure 8a
#
plot_chromvar_clustering <- scater::plotReducedDim(scData_filtered, use_dimred = "chromvar_cluster_specific", colour_by = ".cluster_5", shape_by = "batch", add_ticks = F) +
  ggplot2::ggtitle("Clustering on chromVar sample correlation") +
  ggplot2::guides(shape = FALSE, color = FALSE) +
  ggplot2::theme(plot.title = ggplot2::element_text(face="bold", color="black", size=6), legend.title = ggplot2::element_blank()) +
  ggplot2::scale_shape_manual(values = c(16, 16)) +
  ggplot2::scale_color_manual(values = c("3" = "#1f78b4", "2" = "#6a3d9a", "1" = "#33a02c", "4" = "#e31a1c", "5" = "#ff7f00"))

plot_chromvar_batch <- scater::plotReducedDim(scData_filtered, use_dimred = "chromvar_cluster_specific", colour_by = "batch", shape_by = "batch", add_ticks = F) +
  ggplot2::ggtitle("Clustering on chromVar sample correlation") +
  ggplot2::guides(color = FALSE, shape = FALSE) +
  ggplot2::scale_shape_manual(values = c(16, 16)) +
  ggplot2::theme(plot.title = ggplot2::element_text(face="bold", color="black", size=6), legend.title = ggplot2::element_blank())

figure_8a <- multipanelfigure::multi_panel_figure(rows = 4, columns = 2, width = 205, height = 230)

figure_8a <- multipanelfigure::fill_panel(figure_8a, plot_chromvar_clustering, column = 1:2, row = 1:2)
figure_8a <- multipanelfigure::fill_panel(figure_8a, plot_chromvar_batch, column = 1:2, row = 3:4)

#
# Figure 8c/d
#
figure_dpt_tsne <- multipanelfigure::multi_panel_figure(columns = 3, rows = 3, width = 205, height = 230)

dpts <- c("dpt_cardiac", "dpt_endo")
for(d in dpts) {
  f <-scater::plotReducedDim(scData_filtered, use_dimred = "tsne", colour_by = d, shape_by = "batch", add_ticks = F) +
    ggplot2::ggtitle(paste("LSI:", d)) +
    viridis::scale_color_viridis() +
    ggplot2::guides(shape = FALSE) +
    ggplot2::scale_shape_manual(values = c(16, 16)) +
    ggplot2::theme(plot.title = ggplot2::element_text(face="bold", color="black", size=6),
                   legend.title = ggplot2::element_blank(),
                   legend.key.size =  grid::unit(2, "mm"),
                   legend.margin = ggplot2::margin(b = 2),
                   legend.position = c(0.7,1),
                   legend.justification = c(0,0),
                   legend.background = ggplot2::element_rect(fill = "white"),
                   legend.direction = "horizontal",
                   legend.text = ggplot2::element_text(size = 4)) +
    ggplot2::xlab("") + ggplot2::ylab("")
  figure_dpt_tsne <- multipanelfigure::fill_panel(figure_dpt_tsne, f)
}

ggplot2::ggsave(figure_dpt_tsne,
                file = "figures_2/LSI_dpt.pdf",
                height = 230,
                width = 205,
                unit = "mm")

#
# Figure 8f/g
#
