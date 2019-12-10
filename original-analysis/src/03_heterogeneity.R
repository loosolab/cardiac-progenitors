library(biomaRt)
library(singlecellutils)
library(dplyr)
library(ggplot2)
library(cowplot)
library(UpSetR)

#
# Load the C1 data
#
source("src/parameters.R")
load(file = file.path(parameters$general$path_supdata, "c1_subset.Rdata"))
#save(wg_subset, file = file.path(parameters$general$path_supdata, "wg_subset.Rdata"))

het_subset <- c1_subset[fData(c1_subset)$biotype == "protein_coding" & !fData(c1_subset)$is_feature_control_mito, ]

#
# Assemble list of genes to exclude from being an HVG
#
mart <- useMart("ensembl", dataset = parameters$general$biomart_dataset)
go.cellcycle <- getBM(attributes=c('ensembl_gene_id'), filters = 'go', values = 'GO:0007049', mart = mart)
go.translation <- getBM(attributes=c('ensembl_gene_id'), filters = 'go', values = 'GO:0006412', mart = mart)
go.ribosome <- getBM(attributes=c('ensembl_gene_id'), filters = 'go', values = 'GO:0005840', mart = mart)

exclude <- (rownames(het_subset) %in% unique(c(go.cellcycle$ensembl_gene_id, go.translation$ensembl_gene_id, go.ribosome$ensembl_gene_id)))
het_subset <- het_subset[!exclude, ]

#
# Assemble lineages
#
isl1.lineage <- which(pData(c1_subset)$Background == "isl1")
nkx.lineage <- which(pData(c1_subset)$Background == "nkx2-5")

lineages <- list(isl = isl1.lineage, nkx = nkx.lineage)

#
# Lineage and Time-point specific heterogeneity
#
early <- which(pData(het_subset)$Timepoint == "e7.5")
mid <- which(pData(het_subset)$Timepoint == "e8.5")
late <- which(pData(het_subset)$Timepoint == "e9.5")

time.lineages <- list(isl1 = isl1.lineage, isl1.early = intersect(isl1.lineage, early), isl1.mid = intersect(isl1.lineage, mid), isl1.late = intersect(isl1.lineage, late),
                      nkx = nkx.lineage, nkx.early = intersect(nkx.lineage, early), nkx.mid = intersect(nkx.lineage, mid), nkx.late = intersect(nkx.lineage, late))

means <- lapply(time.lineages, function(x) {
  subset <- het_subset[, x]
  data <- 2^get_exprs(subset, "norm_exprs_lineage")-1

  log(rowMeans(data))/log(10)
})

dropouts <- lapply(time.lineages, function(x) {
  subset <- het_subset[, x]
  data <- 2^get_exprs(subset, "norm_exprs_lineage")-1

  apply(data, 1, dropout.fun)
})

cvs <- lapply(time.lineages, function(x) {
  subset <- het_subset[, x]
  data <- 2^get_exprs(subset, "norm_exprs_lineage")-1

  log10(apply(data, 1, cv2.fun))
})

dependencies <- data.frame(dataset = rep(c("Isl1", "Nkx2-5"), each=nrow(het_subset)*4), timepoint = rep(c("all","e7.5","e8.5","e9.5","all","e7.5","e8.5","e9.5"), each=nrow(het_subset)), mean = unlist(means), dropout = unlist(dropouts), cv = unlist(cvs))
rownames(dependencies) <- c(
  paste0("isl1",rep(c("all","e7.5","e8.5","e9.5"), each=nrow(het_subset)), rep(rownames(het_subset), times=4)),
  paste0("nkx",rep(c("all","e7.5","e8.5","e9.5"), each=nrow(het_subset)), rep(rownames(het_subset), times=4)))

het <- lapply(time.lineages, function(l) {
  subset <- het_subset[, l]
  data <- 2^get_exprs(subset, "norm_exprs_lineage")-1

  dro <- apply(data, 1, dropout.fun)
  cvs <- log10(apply(data, 1, cv2.fun))

  dm.do <- heterogeneity(data, statistic = "mean", order_by = dro, normalization = "windows", window = 200)
  dm.cv <- heterogeneity(data, statistic = "mean", order_by = cvs, normalization = "windows", window = 200)

  unique(c(names(dm.do[which(dm.do > parameters$heterogeneity$het_zscore)]), names(dm.cv[which(dm.cv > parameters$heterogeneity$het_zscore)])))
})

lapply(names(het), function(l) {
  write.table(het[[l]], quote = F, row.names = F, col.names = F, file = paste0(parameters$general$path_supdata,"/heterogeneity-",l,".txt"))
})

is.het <- c(paste0("isl1all",het[["isl1"]]), paste0("isl1e7.5",het[["isl1.early"]]), paste0("isl1e8.5",het[["isl1.mid"]]), paste0("isl1e9.5",het[["isl1.late"]]),
            paste0("nkxall",het[["nkx"]]), paste0("nkxe7.5",het[["nkx.early"]]), paste0("nkxe8.5",het[["nkx.mid"]]), paste0("nkxe9.5",het[["nkx.late"]]))
dependencies$het <- ifelse(rownames(dependencies) %in% is.het, T, F)

#
# Supplementary Figure S2
#
theme2 <- theme(plot.background = element_blank(), panel.grid.major = element_line(size=.2, colour = "grey"),panel.grid.minor = element_line(size=.1, colour = "grey"),
                panel.border = element_blank(),panel.background = element_blank(),axis.line.x = element_line(size=.3),axis.line.y = element_line(size=.3), legend.text = element_text(size=6),
                #axis.title.x = element_blank(), axis.title.y = element_blank(),
                plot.title = element_text(face="bold", color="black", size=6), legend.key.size =  unit(2, "mm"), legend.margin=unit(-25, "mm"), axis.text = element_text(size=6))

s2_isl_do <- dependencies %>%
  filter(dataset == "Isl1", timepoint == "all") %>%
  ggplot(aes_string(x="mean", y="dropout", color = "het")) +
  geom_point(size=0.2, aes(alpha=0.3)) +
  scale_colour_manual(values = c("black", "red")) +
  #ggtitle(paste0("Isl1 heterogeneity")) +
  #xlab("Mean gene expression") +
  #ylab("Dropout rate") +
  theme2 + guides(alpha=FALSE, color=FALSE, size=FALSE)
s2_isl_cv <- dependencies %>%
  filter(dataset == "Isl1", timepoint == "all") %>%
  ggplot(aes_string(x="mean", y="cv", color = "het")) +
  geom_point(size=0.2, aes(alpha=0.3)) +
  scale_colour_manual(values = c("black", "red")) +
  #ggtitle(paste0("Isl1 heterogeneity")) +
  #xlab("Mean gene expression") +
  #ylab("Coefficient of variation") +
  theme2 + guides(alpha=FALSE, color=FALSE, size=FALSE)

s2_isl_het <- plot_grid(s2_isl_cv, s2_isl_do)

s2_nkx_do <- dependencies %>%
  filter(dataset == "Isl1", timepoint == "all") %>%
  ggplot(aes_string(x="mean", y="dropout", color="het")) +
  geom_point(size=0.2, aes(alpha=0.3)) +
  scale_colour_manual(values = c("black", "red")) +
  #ggtitle(paste0("Nkx2-5 heterogeneity")) +
  theme2 + guides(alpha=FALSE, color=FALSE, size=FALSE)
s2_nkx_cv <- dependencies %>%
  filter(dataset == "Nkx2-5", timepoint == "all") %>%
  ggplot(aes_string(x="mean", y="cv", color="het")) +
  geom_point(size=0.2, aes(alpha=0.3)) +
  scale_colour_manual(values = c("black", "red")) +
  #ggtitle(paste0("Nkx2-5 heterogeneity")) +
  theme2 + guides(alpha=FALSE, color=FALSE, size=FALSE)

s2_nkx_het <- plot_grid(s2_nkx_cv, s2_nkx_do)
p <- grid.grabExpr(upset(fromList(het[c("isl1", "isl1.early", "isl1.late", "isl1.mid")]), order.by = "freq"), wrap = T)
q <- grid.grabExpr(upset(fromList(het[c("nkx", "nkx.early", "nkx.late", "nkx.mid")]), order.by = "freq"), wrap = T)

m <- multi_panel_figure(width = 205, height = 230, columns = 2, rows = 4)
m <- fill_panel(m, s2_isl_het, row = 1, col = 1)
m <- fill_panel(m, s2_nkx_het, row = 1, col = 2)
m <- fill_panel(m, p, row = c(2,3,4), col = 1)
m <- fill_panel(m, q, row = c(2,3,4), col = 2)

#
# Load the WG data
#
source("src/parameters.R")
load(file = file.path(parameters$general$path_supdata, "wg_subset.Rdata"))
#save(wg_subset, file = file.path(parameters$general$path_supdata, "wg_subset.Rdata"))

het_subset <- wg_subset[fData(wg_subset)$biotype == "protein_coding" & !fData(wg_subset)$is_feature_control_mito, ]

#
# Assemble list of genes to exclude from being an HVG
#
mart <- useMart("ensembl", dataset = parameters$general$biomart_dataset)
go.cellcycle <- getBM(attributes=c('ensembl_gene_id'), filters = 'go', values = 'GO:0007049', mart = mart)
go.translation <- getBM(attributes=c('ensembl_gene_id'), filters = 'go', values = 'GO:0006412', mart = mart)
go.ribosome <- getBM(attributes=c('ensembl_gene_id'), filters = 'go', values = 'GO:0005840', mart = mart)

exclude <- (rownames(het_subset) %in% unique(c(go.cellcycle$ensembl_gene_id, go.translation$ensembl_gene_id, go.ribosome$ensembl_gene_id)))
het_subset <- het_subset[!exclude, ]

#
# Assemble lineages
#
isl1.lineage <- which(pData(wg_subset)$Background == "isl1")
nkx.lineage <- which(pData(wg_subset)$Background == "nkx2-5")

#
# Lineage and Time-point specific heterogeneity
#
early <- which(pData(het_subset)$Timepoint == "e7.5")
mid <- which(pData(het_subset)$Timepoint == "e8.5")
late <- which(pData(het_subset)$Timepoint == "e9.5")

time.lineages <- list(isl1 = isl1.lineage, isl1.early = intersect(isl1.lineage, early), isl1.late = intersect(isl1.lineage, late),
                      nkx = nkx.lineage, nkx.early = intersect(nkx.lineage, early), nkx.mid = intersect(nkx.lineage, mid), nkx.late = intersect(nkx.lineage, late))

means <- lapply(time.lineages, function(x) {
  subset <- het_subset[, x]
  data <- 2^get_exprs(subset, "norm_exprs_lineage")-1

  log(rowMeans(data))/log(10)
})

dropouts <- lapply(time.lineages, function(x) {
  subset <- het_subset[, x]
  data <- 2^get_exprs(subset, "norm_exprs_lineage")-1

  apply(data, 1, dropout.fun)
})

cvs <- lapply(time.lineages, function(x) {
  subset <- het_subset[, x]
  data <- 2^get_exprs(subset, "norm_exprs_lineage")-1

  log10(apply(data, 1, cv2.fun))
})

dependencies <- data.frame(dataset = c(rep("Isl1", times=nrow(het_subset)*3), rep("Nkx2-5", times=nrow(het_subset)*4)),
                           timepoint = rep(c("all","e7.5","e9.5","all","e7.5","e8.5","e9.5"), each=nrow(het_subset)),
                           mean = unlist(means),
                           dropout = unlist(dropouts),
                           cv = unlist(cvs)
                           )
rownames(dependencies) <- c(
  paste0("isl1",rep(c("all","e7.5","e9.5"), each=nrow(het_subset)), rep(rownames(het_subset), times=3)),
  paste0("nkx",rep(c("all","e7.5","e8.5","e9.5"), each=nrow(het_subset)), rep(rownames(het_subset), times=4)))

het <- lapply(time.lineages, function(l) {
  subset <- het_subset[, l]
  data <- 2^get_exprs(subset, "norm_exprs_lineage")-1

  dro <- apply(data, 1, dropout.fun)
  cvs <- log10(apply(data, 1, cv2.fun))

  dm.do <- heterogeneity(data, statistic = "mean", order_by = dro, normalization = "windows", window = 200)
  dm.cv <- heterogeneity(data, statistic = "mean", order_by = cvs, normalization = "windows", window = 200)

  unique(c(names(dm.do[which(dm.do > parameters$heterogeneity$het_zscore)]), names(dm.cv[which(dm.cv > parameters$heterogeneity$het_zscore)])))
})

lapply(names(het), function(l) {
  write.table(het[[l]], quote = F, row.names = F, col.names = F, file = paste0(parameters$general$path_supdata,"/heterogeneity-wg-",l,".txt"))
})

is.het <- c(paste0("isl1all",het[["isl1"]]), paste0("isl1e7.5",het[["isl1.early"]]), paste0("isl1e8.5",het[["isl1.mid"]]), paste0("isl1e9.5",het[["isl1.late"]]),
            paste0("nkxall",het[["nkx"]]), paste0("nkxe7.5",het[["nkx.early"]]), paste0("nkxe8.5",het[["nkx.mid"]]), paste0("nkxe9.5",het[["nkx.late"]]))
dependencies$het <- ifelse(rownames(dependencies) %in% is.het, T, F)
