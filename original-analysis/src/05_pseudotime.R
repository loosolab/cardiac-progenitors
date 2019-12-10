library(scater)
library(destiny)
library(singlecellutils)
library(RColorBrewer)
library(multipanelfigure)

#
# Load data
#
source("src/parameters.R")
load(file.path(parameters$general$path_supdata, "c1_subset.Rdata"))

#
# Nkx2-5 diffusion map
#
nkx_cells <- which(pData(c1_subset)$Background == "nkx2-5" & pData(c1_subset)$Platform == "C1")
nkx_de_genes <- which(fData(c1_subset)$nkx_de & fData(c1_subset)$nkx_marker)

data <- get_exprs(c1_subset[nkx_de_genes, nkx_cells], "norm_exprs_sf")
set.seed(1001)
nkx_dm <- DiffusionMap(t(data), distance = 'euclidean')

tip_cell <- find_tips(nkx_dm)[which(pData(c1_subset[, nkx_cells])$cluster[find_tips(nkx_dm)] == 1)]

nkx_dpt <- DPT(nkx_dm, tips = tip_cell)
nkx_pseudotime <- nkx_dpt[tip_cell]

#
# Save dpt and diffusion components to object
#
pData(c1_subset)$dpt <- NA
pData(c1_subset)$dm1 <- NA
pData(c1_subset)$dm2 <- NA

pData(c1_subset)$dpt[nkx_cells] <- nkx_pseudotime
pData(c1_subset)$dm1[nkx_cells] <- nkx_dm@eigenvectors[, 1]
pData(c1_subset)$dm2[nkx_cells] <- nkx_dm@eigenvectors[, 2]

save(nkx_dm, file = file.path(parameters$general$path_supdata, "Nkx2-5-DM.Rdata"))

#
# Isl1 diffusion map
#
isl1_cells <- which(pData(c1_subset)$Background == "isl1" & pData(c1_subset)$Platform == "C1")
isl1_de_genes <- which(fData(c1_subset)$isl1_de & fData(c1_subset)$isl1_marker)

data <- get_exprs(c1_subset[isl1_de_genes, isl1_cells], "norm_exprs_sf")
set.seed(1002)
isl1_dm <- DiffusionMap(t(data), distance = 'e')

tip_cell <- find_tips(isl1_dm)[which(pData(c1_subset[, isl1_cells])$cluster[find_tips(isl1_dm)] == 5)]

isl_dpt <- DPT(isl1_dm, tips = tip_cell)
isl_pseudotime <- isl_dpt[tip_cell,]

#
# Save dpt and diffusion components to object
#
pData(c1_subset)$dpt[isl1_cells] <- isl_pseudotime
pData(c1_subset)$dm1[isl1_cells] <- isl1_dm@eigenvectors[, 1]
pData(c1_subset)$dm2[isl1_cells] <- isl1_dm@eigenvectors[, 3]
pData(c1_subset)$dm3[isl1_cells] <- isl1_dm@eigenvectors[, 2]

save(isl1_dm, file = file.path(parameters$general$path_supdata, "Isl1-DM.Rdata"))
save(c1_subset, file = file.path(parameters$general$path_supdata, "c1_subset.Rdata"))
