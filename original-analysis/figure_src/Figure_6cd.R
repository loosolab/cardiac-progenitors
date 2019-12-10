library(destiny)
library(scater)
library(singlecellutils)
library(multipanelfigure)

#
# Graphical settings
#
global.settings <- list(layout.heights=list(main.key.padding = -2, bottom.padding = 0, axis.xlab.padding = -0.5), layout.widths = list(ylab.axis.padding = -0.5))

#
# Load neccessary data
#
source("src/parameters.R")
load(file.path(parameters$general$path_supdata, "c1_subset.Rdata"))
load(file.path(parameters$general$path_supdata, "wg_subset.Rdata"))
load(file.path(parameters$general$path_supdata, "Isl1-DM.Rdata"))
load(file.path(parameters$general$path_supdata, "Nkx2-5-DM.Rdata"))

nkx_cells <- which(pData(c1_subset)$Background == "nkx2-5" & pData(c1_subset)$Platform == "C1")
isl1_cells <- which(pData(c1_subset)$Background == "isl1" & pData(c1_subset)$Platform == "C1")
nkx_dm_coords <- cbind(-pData(c1_subset[, nkx_cells])$dm1, -pData(c1_subset[, nkx_cells])$dm2)
isl1_dm_coords <- cbind(-pData(c1_subset[, isl1_cells])$dm1, -pData(c1_subset[, isl1_cells])$dm2, -pData(c1_subset[, isl1_cells])$dm3)

#
# Predict Nkx-OE onto Nkx2-5 trajectory
#
nkxoe_cells <- which(pData(c1_subset)$Background == "nkx2-5oe" & pData(c1_subset)$Platform == "C1")
m <- match(colnames(nkx_dm@data_env$data), rownames(c1_subset))

nkxoe_data <- get_exprs(c1_subset[m, nkxoe_cells], "norm_exprs_sf")

nkxoe_dcs <- dm_predict(nkx_dm, t(nkxoe_data))
nkxoe_coords <- matrix(nkxoe_dcs@x, ncol = nkxoe_dcs@Dim[2])
nkx_dm_coords_oe <- rbind(nkx_dm_coords, cbind(-nkxoe_coords[, 1], -nkxoe_coords[, 2]))

cluster <- c(pData(c1_subset[, nkx_cells])$cluster, rep("OE", nrow(nkxoe_coords)))
rm_obs <- is.na(cluster)
cluster <- paste0("cluster", cluster[!rm_obs])

p_nkx_nkxoe <- colorAMap(nkx_dm_coords_oe[!rm_obs,], colour_by = factor(cluster), palette = parameters$colors$nkxko_cluster_palette, main = list(label = "Nkx2-5 overexpression prediction", cex = 0.75), scales = list(cex = 0.4, tck = c(0.5, 0)), xlab = list(label = "Dimension 1", cex = 0.6), ylab = list(label = "Dimension 2", cex = 0.6), par.settings = global.settings)

#
# Predict Nkx2-5OE onto Isl1 trajectory
#
m <- match(colnames(isl1_dm@data_env$data), rownames(c1_subset))
nkxoe_cells <- which(pData(c1_subset)$Background == "nkx2-5oe" & pData(c1_subset)$Platform == "C1")
nkxoe_data <- get_exprs(c1_subset[m, nkxoe_cells], "norm_exprs_sf")

nkxoe_dcs <- dm_predict(isl1_dm, t(nkxoe_data))
nkxoe_coords <- matrix(nkxoe_dcs@x, ncol = nkxoe_dcs@Dim[2])
isl1_dm_coords_oe <- rbind(isl1_dm_coords[, c(1,2)], cbind(-nkxoe_coords[, 1], -nkxoe_coords[, 3]))

cluster <- c(pData(c1_subset[, isl1_cells])$cluster, rep("OE", nrow(nkxoe_coords)))
rm_obs <- is.na(cluster)
cluster <- paste0("cluster", cluster[!rm_obs])

p_isl1_oe <- colorAMap(isl1_dm_coords_oe[!rm_obs,], colour_by = factor(cluster), palette = parameters$colors$islko_cluster_palette, ylab = list(label = "Dimension 3", cex = 0.6), main = list(label = "Isl1 Nkx2-5 overexpression prediction", cex = 0.75), scales = list(cex = 0.4, tck = c(0.5, 0)), xlab = list(label = "Dimension 1", cex = 0.6), par.settings = global.settings)

#
# Save
#
m <- multi_panel_figure(width = 205, height = 230, columns = 2, rows = 3)
m <- fill_panel(m, p_nkx_nkxoe)
m <- fill_panel(m, p_isl1_oe)

ggsave(
  plot = m,
  filename = file.path(parameters$general$path_rfigures, "Figure6_cd.pdf"),
  width = 205,
  height = 230,
  units = "mm",
  dpi = 600)
