library(scran)
library(scater)
library(SCnorm)
library(multipanelfigure)
library(gridExtra)
library(limma)

#
# Load data
#
source("src/parameters.R")
load(file.path(parameters$general$path_supdata, "scd.RData"))

#
# Assign cell cycle stage
#
pairs <- readRDS(system.file("exdata", parameters$normalisation$cycle_markers, package="scran"))
cellcycle <- cyclone(scd, pairs, gene.names = rownames(scd), assay="exprs", verbose=F)

pData(scd)$cellcycle <- factor(cellcycle$phases)
plot(cellcycle$score$G1, cellcycle$score$G2M, pch=16)

cc_tab <- data.frame(count = table(cellcycle$phases))
colnames(cc_tab) <- c("Cell cycle stage", "Cells")

tt1 <- ttheme_minimal(core=list(fg_params=list(hjust=1, x=0.9)),
                      rowhead=list(fg_params=list(hjust=1, x=0.9)),
                      colhead = list(fg_params=list(hjust=1, x=0.9)))

cell_cycle_table <- tableGrob(cc_tab, theme=tt1, rows = NULL)

#
# Perfom lineage-specific TMM normalization with batch correction
#
# Prepare lineages
lineages <- list(isl1 = pData(scd)$Background == "isl1" | pData(scd)$Background == "isl1ko" | pData(scd)$Background == "nkx2-5oe",
                 nkx = pData(scd)$Background == "nkx2-5")

sf.normalized <- lapply(lineages, function(lineage_cells) {
  lineage <- scd[, lineage_cells]

  lineage <- computeSumFactors(lineage, clusters = pData(lineage)$Timepoint)
  lineage <- normalise(lineage)
  get_exprs(lineage, "norm_exprs")
})

tmm.normalized <- lapply(lineages, function(lineage_cells) {
  lineage <- scd[, lineage_cells]

  # Create design matrix
  #design <- matrix(c(rep(1, ncol(lineage)), as.numeric(pData(lineage)$Timepoint)), nrow = ncol(lineage))

  # Remove batch effect
  #corrected <- removeBatchEffect(get_exprs(lineage, "exprs"), design = design, block = pData(lineage)$Platform)

  # Set corrected expression
  #set_exprs(lineage, "exprs") <- corrected

  # Normalize lineages
  lineage <- normalizeExprs(lineage, method="TMM", return_norm_as_exprs = FALSE)

  # Create design matrix
  #design <- matrix(c(rep(1, ncol(lineage)), as.numeric(pData(lineage)$Timepoint)), nrow = ncol(lineage))

  # Remove batch effect
  #tmm <- removeBatchEffect(get_exprs(lineage, "norm_exprs"), design = design, block = pData(lineage)$Platform)

  get_exprs(lineage, "norm_exprs")
})

# Establish subsets
isl1_subset <- scd[, lineages$isl1]
set_exprs(isl1_subset, "norm_exprs_lineage") <- tmm.normalized$isl1
set_exprs(isl1_subset, "norm_exprs_sf") <- sf.normalized$isl1

nkx_subset <- scd[, lineages$nkx]
set_exprs(nkx_subset, "norm_exprs_lineage") <- tmm.normalized$nkx
set_exprs(nkx_subset, "norm_exprs_sf") <- sf.normalized$nkx

#
# Save created objects for further analysis
#
save(isl1_subset, file = file.path(parameters$general$path_supdata, "Isl1_subset.Rdata"))
save(nkx_subset, file = file.path(parameters$general$path_supdata, "Nkx2-5_subset.Rdata"))

#
# Prepare data platform subsets
#
subsets <- list(c1 = pData(scd)$Platform == "C1",
                wg = pData(scd)$Platform == "WG")

#
# Perfom lineage- and platform-specific TMM normalisation
#

sf.normalized <- lapply(subsets, function(platform) {
  isl_lineage <- platform & (pData(scd)$Background == "isl1" | pData(scd)$Background == "isl1ko" | pData(scd)$Background == "nkx2-5oe")
  nkx_lineage <- platform & (pData(scd)$Background == "nkx2-5")

  # Normalize lineages
  isl_subset <- scd[, isl_lineage]
  isl_subset <- computeSumFactors(isl_subset)#, clusters = pData(isl_subset)$Timepoint)
  isl_subset <- normalise(isl_subset)

  nkx_subset <- scd[, nkx_lineage]
  nkx_subset <- computeSumFactors(nkx_subset)#, clusters = pData(nkx_subset)$Timepoint)
  nkx_subset <- normalise(nkx_subset)

  # Bind results and return
  sf <- cbind(get_exprs(isl_subset, "norm_exprs"), get_exprs(nkx_subset, "norm_exprs"))
  m <- match(colnames(scd), colnames(sf))
  sf[, na.omit(m)]
})

c1_subset <- scd[, subsets$c1]
set_exprs(c1_subset, "norm_exprs_sf") <- sf.normalized$c1

wg_subset <- scd[, subsets$wg]
set_exprs(wg_subset, "norm_exprs_sf") <- sf.normalized$wg

tmm.normalized <- lapply(subsets, function(platform) {
  # Establish lineages
  isl_lineage <- platform & (pData(scd)$Background == "isl1" | pData(scd)$Background == "isl1ko" | pData(scd)$Background == "nkx2-5oe")
  nkx_lineage <- platform & (pData(scd)$Background == "nkx2-5")

  # Normalize lineages
  isl_subset <- scd[, isl_lineage]
  isl_subset <- normalizeExprs(isl_subset, method="TMM", return_norm_as_exprs = FALSE)

  nkx_subset <- scd[, nkx_lineage]
  nkx_subset <- normalizeExprs(nkx_subset, method="TMM", return_norm_as_exprs = FALSE)

  # Bind results and return
  tmm <- cbind(get_exprs(isl_subset, "norm_exprs"), get_exprs(nkx_subset, "norm_exprs"))
  m <- match(colnames(scd), colnames(tmm))
  tmm[, na.omit(m)]
})

# Establish subsets
#c1_subset <- scd[, subsets$c1]
set_exprs(c1_subset, "norm_exprs_lineage") <- tmm.normalized$c1

#wg_subset <- scd[, subsets$wg]
set_exprs(wg_subset, "norm_exprs_lineage") <- tmm.normalized$wg

#
# Save created objects for further analysis
#
save(c1_subset, file = file.path(parameters$general$path_supdata, "c1_subset.Rdata"))
save(wg_subset, file = file.path(parameters$general$path_supdata, "wg_subset.Rdata"))

