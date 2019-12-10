library(scater)
library(singlecellutils)
library(Rtsne)
library(multipanelfigure)
library(dplyr)

#
# Load data
#
source("src/parameters.R")
source("src/DE.R")
load(file.path(parameters$general$path_supdata, "c1_subset.Rdata"))

#
# Load heterogeneous genes
#
# Add variable genes from C1 analysis
het_early <- readLines(file.path(parameters$general$path_supdata, "heterogeneity-nkx.early.txt"))
het_mid <- readLines(file.path(parameters$general$path_supdata, "heterogeneity-nkx.mid.txt"))
het_late <- readLines(file.path(parameters$general$path_supdata, "heterogeneity-nkx.late.txt"))
het_pool <- readLines(file.path(parameters$general$path_supdata, "heterogeneity-nkx.txt"))

nkx_het <- unique(c(het_early, het_mid, het_late, het_pool))

het_early <- readLines(file.path(parameters$general$path_supdata, "heterogeneity-isl1.early.txt"))
het_mid <- readLines(file.path(parameters$general$path_supdata, "heterogeneity-isl1.mid.txt"))
het_late <- readLines(file.path(parameters$general$path_supdata, "heterogeneity-isl1.late.txt"))
het_pool <- readLines(file.path(parameters$general$path_supdata, "heterogeneity-isl1.txt"))

isl1_het <- unique(c(het_early, het_mid, het_late, het_pool))

#
# Perform SOM clustering
#
nkx_cells <- which(pData(c1_subset)$Background == "nkx2-5" & pData(c1_subset)$Platform == "C1")
nkx_expression <- get_exprs(c1_subset[, nkx_cells], "norm_exprs_sf")
nkx_expression_scaled <- t(scale(t(nkx_expression[nkx_het, ])))

nkx_som <- calcSOM(nkx_expression_scaled, train = 1:nrow(nkx_expression_scaled), num_epochs = 2000, seed = 1004)
nkx_som$codes <- nkx_som$codes[[1]]

save(nkx_som, file = file.path(parameters$general$path_supdata, "Nkx2-5_SOM.Rdata"))

isl1_cells <- which(pData(c1_subset)$Background == "isl1" & pData(c1_subset)$Platform == "C1")
isl1_cells <- which((pData(c1_subset)$Background == "isl1" | pData(c1_subset)$Background == "isl1ko") & pData(c1_subset)$Platform == "C1")
isl1_expression <- get_exprs(c1_subset[, isl1_cells], "norm_exprs_sf")
isl1_expression_scaled <- t(scale(t(isl1_expression[isl1_het, ])))

isl1_som <- calcSOM(isl1_expression_scaled, train = 1:nrow(isl1_expression_scaled), num_epochs = 2000, seed = 1004)
isl1_som$codes <- isl1_som$codes[[1]]

save(isl1_som, file = file.path(parameters$general$path_supdata, "Isl1_SOM.Rdata"))

hdbscan <- function(data, min_samples = 7L, min_cluster_size = 9L, outlier = 0, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  h <- reticulate::import("hdbscan")
  cl <- h$HDBSCAN(min_samples = min_samples, min_cluster_size = min_cluster_size)
  labels <- cl$fit_predict(data) + 1
  labels[labels == 0] <- outlier
  return(factor(labels))
}

#
# Nkx2-5 dimension reduction and clustering
#
seed <- 2596
set.seed(seed)
nkx_tsne <- Rtsne(t(nkx_som$codes), perplexity = 15, theta = 0.05, max_iter = 2000)

nkx_cl <- hdbscan(nkx_tsne$Y, min.samples = 7, min.cluster.size = 9)
nkx_cl$.cluster[nkx_cl$.cluster == 0] <- NA
nkx_cl$.cluster <- factor(nkx_cl$.cluster)

colorAMap(nkx_tsne$Y, shape_by = factor(pData(c1_subset[, nkx_cells])$Platform), colour_by = nkx_cl$.cluster, pch=16)
colorAMap(nkx_tsne$Y, shape_by = factor(pData(c1_subset[, nkx_cells])$Platform), colour_by = pData(c1_subset[, nkx_cells])$Timepoint, pch=16)
colorAMap(nkx_tsne$Y, colour_by = factor(pData(c1_subset[, nkx_cells])$Platform), shape_by = nkx_cl$.cluster, pch=16)
colorAMap(nkx_tsne$Y, shape_by = factor(pData(c1_subset[, nkx_cells])$Platform), colour_by = nkx_cl$.cluster, pch=16)

pData(c1_subset)$cluster <- NA
pData(c1_subset)$tsne1 <- NA
pData(c1_subset)$tsne2 <- NA
pData(c1_subset)$cluster[nkx_cells] <- nkx_cl$.cluster
pData(c1_subset)$tsne1[nkx_cells] <- nkx_cl$X1
pData(c1_subset)$tsne2[nkx_cells] <- nkx_cl$X2


#
# Isl1 dimension reduction and clustering
#
seed <- 3465
set.seed(seed)
isl1_tsne <- Rtsne(t(isl1_som$codes), perplexity = 15, theta = 0.05, max_iter = 2000)

isl1_cl <- hdbscan(isl1_tsne$Y, min_samples = 7L, min_cluster_size = 9L)
isl1_cl[isl1_cl == 0] <- NA
isl1_cl <- factor(isl1_cl)

# i <- which(pData(c1_subset)$Background == "isl1ko" & pData(c1_subset)$Platform == "C1")
# cluster_ko <- factor(pData(c1_subset[, isl1_cells])$cluster, levels = c("1", "2", "3", "4", "5", "KO"))
# cluster_ko[i] <- "KO"
# colorAMap(isl1_tsne$Y, shape_by = factor(pData(c1_subset[, isl1_cells])$Platform), colour_by = cluster_ko, pch=16, palette = parameters$colors$islko_cluster_palette)
colorAMap(isl1_tsne$Y, shape_by = factor(pData(c1_subset[, isl1_cells])$Platform), colour_by = isl1_cl, pch=16)
colorAMap(isl1_tsne$Y, shape_by = factor(pData(c1_subset[, isl1_cells])$Platform), colour_by = pData(c1_subset[, isl1_cells])$Timepoint, pch=16)
colorAMap(isl1_tsne$Y, shape_by = factor(pData(c1_subset[, isl1_cells])$Platform), colour_by = isl1_cl, pch=16)
colorAMap(isl1_tsne$Y, shape_by = factor(pData(c1_subset[, isl1_cells])$Platform), colour_by = droplevels(pData(c1_subset[, isl1_cells])$Background), pch=16)


pData(c1_subset)$cluster[isl1_cells] <- isl1_cl$.cluster
pData(c1_subset)$tsne1[isl1_cells] <- isl1_cl$X1
pData(c1_subset)$tsne2[isl1_cells] <- isl1_cl$X2

#
# Nkx2-5 differential expression
#
nkx_contrasts <- list(
  cluster1_rest = list(which(nkx_cl$.cluster == 1), which(nkx_cl$.cluster != 1)),
  cluster2_rest = list(which(nkx_cl$.cluster == 2), which(nkx_cl$.cluster != 2)),
  cluster3_rest = list(which(nkx_cl$.cluster == 3), which(nkx_cl$.cluster != 3))
)

nkx_diff_data <- differentialExpression(nkx_expression, contrasts = nkx_contrasts, fData(c1_subset))
names(nkx_diff_data) <- names(nkx_contrasts)

# Nkx2-5 marker
nkx_markers <- get_marker_genes(nkx_expression, nkx_cl$.cluster)

save(nkx_diff_data, file = file.path(parameters$general$path_supdata, "Nkx2-5-diffExprs.Rdata"))
save(nkx_markers, file = file.path(parameters$general$path_supdata, "Nkx2-5-markers.Rdata"))

# Add information to object
nkx_de_genes_l <- lapply(names(nkx_diff_data), function(c) {
  data <- nkx_diff_data[[c]]
  # Add marker gene data
  data %>%
    left_join(nkx_markers, by = "geneID") %>%
    dplyr::rename(marker_pval = pvalue, marker_fdr = fdr.y) -> data_marker
  write.table(data_marker, file = file.path(parameters$general$path_supdata, paste0("differentialExpression/Nkx2-5_",c,".txt")), sep="\t", quote = F, row.names = F, col.names = T)
  i <- which(data$fdr < parameters$diffusionmaps$de_fdr & data$biotype == "protein_coding" & (data$lfc.hi < -parameters$diffusionmaps$de_abs_logfc | data$lfc.lo > parameters$diffusionmaps$de_abs_logfc))
  data$geneID[i]
})
nkx_de_genes <- unique(unlist(nkx_de_genes_l))

m <- match(nkx_de_genes, rownames(fData(c1_subset)))
fData(c1_subset)$nkx_de <- FALSE
fData(c1_subset)$nkx_de[m] <- TRUE

i <- which(nkx_markers$auroc > parameters$diffusionmaps$marker_auroc & nkx_markers$fdr < parameters$diffusionmaps$marker_fdr)
nkx_marker_genes <- nkx_markers$geneID[i]

m <- match(nkx_marker_genes, rownames(fData(c1_subset)))
fData(c1_subset)$nkx_marker <- FALSE
fData(c1_subset)$nkx_marker[m] <- TRUE

#
# Isl1 differential expression
#
isl1_contrasts <- list(
  cluster1_rest = list(which(isl1_cl$.cluster == 1), which(isl1_cl$.cluster != 1)),
  cluster2_rest = list(which(isl1_cl$.cluster == 2), which(isl1_cl$.cluster != 2)),
  cluster3_rest = list(which(isl1_cl$.cluster == 3), which(isl1_cl$.cluster != 3)),
  cluster4_rest = list(which(isl1_cl$.cluster == 4), which(isl1_cl$.cluster != 4)),
  cluster5_rest = list(which(isl1_cl$.cluster == 5), which(isl1_cl$.cluster != 5))
)

isl1_diff_data <- differentialExpression(isl1_expression, contrasts = isl1_contrasts, fData(c1_subset))
names(isl1_diff_data) <- names(isl1_contrasts)

# Isl1 marker
isl1_markers <- get_marker_genes(isl1_expression, isl1_cl$.cluster)

save(isl1_diff_data, file = file.path(parameters$general$path_supdata, "Isl1-diffExprs.Rdata"))
save(isl1_markers, file=file.path(parameters$general$path_supdata, "Isl1-markers.Rdata"))

# Add information to object
isl1_de_genes_l <- lapply(names(isl1_diff_data), function(c) {
  data <- isl1_diff_data[[c]]
  # Add marker gene data
  data %>%
    left_join(nkx_markers, by = "geneID") %>%
    dplyr::rename(marker_pval = pvalue, marker_fdr = fdr.y) -> data_marker
  write.table(data_marker, file = file.path(parameters$general$path_supdata, paste0("differentialExpression/Isl1_",c,".txt")), sep="\t", quote = F, row.names = F, col.names = T)
  i <- which(data$fdr < parameters$diffusionmaps$de_fdr & data$biotype == "protein_coding" & (data$lfc.hi < -parameters$diffusionmaps$de_abs_logfc | data$lfc.lo > parameters$diffusionmaps$de_abs_logfc))
  data$geneID[i]
})
isl1_de_genes <- unique(unlist(isl1_de_genes_l))

m <- match(isl1_de_genes, rownames(fData(c1_subset)))
fData(c1_subset)$isl1_de <- FALSE
fData(c1_subset)$isl1_de[m] <- TRUE

i <- which(isl1_markers$auroc > parameters$diffusionmaps$marker_auroc & isl1_markers$fdr < parameters$diffusionmaps$marker_fdr)
isl1_marker_genes <- isl1_markers$geneID[i]

m <- match(isl1_marker_genes, rownames(fData(c1_subset)))
fData(c1_subset)$isl1_marker <- FALSE
fData(c1_subset)$isl1_marker[m] <- TRUE

#
# Save back subset
#
save(c1_subset, file = file.path(parameters$general$path_supdata, "c1_subset.Rdata"))
