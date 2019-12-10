library(dplyr)
library(biomaRt)
library(scater)
library(scran)

source("src/parameters.R")
#
# Load Wu et. al., Nkx2-5 WT data
#
e85_attributes <- read.table("ext_data/e8.5_wt_cell_attributes2.csv", sep=";", header=T)
e95_attributes <- read.table("ext_data/e9.5_wt_cell_attributes.csv", sep=",", header=T)
e105_attributes <- read.table("ext_data/e10.5_wt_cell_attributes.csv", sep=",", header=T)

wt_metadata <- bind_rows(e85_attributes, e95_attributes, e105_attributes)
wt_metadata$timepoint <- c(rep("e8.5", nrow(e85_attributes)), rep("e9.5", nrow(e95_attributes)), rep("e10.5", nrow(e105_attributes)))

e85_counts <- read.table("ext_data/e8.5_wt_counts.csv", sep=",", header=T, row.names = 1)
e95_counts <- read.table("ext_data/e9.5_wt_counts.csv", sep=",", header=T, row.names = 1)
e105_counts <- read.table("ext_data/e10.5_wt_counts.csv", sep=",", header=T, row.names = 1)

wt_counts <- counts <- cbind(e85_counts, e95_counts, e105_counts)

#
# Load Wu et. al., Nkx2-5 KO data
#
ko_counts <- read.csv("ext_data/e9.5_Nkx2.5mut_counts.csv", row.names = 1)
ko_metadata <- read.csv("ext_data/e9.5_Nkx2.5mut_cell_attributes.csv")

ko_metadata$timepoint <- "e9.5"
ko_metadata$batch <- "KO"

#
# Merge counts and metadata
#
assertthat::are_equal(rownames(wt_counts), rownames(ko_counts))
wu_counts <- cbind(wt_counts, ko_counts)

wu_metadata <- bind_rows(wt_metadata, ko_metadata)
rownames(wu_metadata) <- make.names(wu_metadata$index)

assertthat::are_equal(rownames(wu_metadata), colnames(wu_counts))

# Filter cells based on quality in metadata
wu_counts <- wu_counts[, wu_metadata$quality]
wu_metadata <- wu_metadata[wu_metadata$quality, ]

# Filter out genes that are not in our data
load(file.path(parameters$general$path_supdata, "c1_subset.Rdata"))
m <- match(rownames(c1_subset), rownames(wu_counts))

wu_counts <- wu_counts[na.omit(m), ]

#
# Fetch featureData from bioMart
#
mart <- useMart("ensembl", dataset = parameters$general$biomart_dataset)
bm <- getBM(c("ensembl_gene_id","external_gene_name","gene_biotype","description"), filters = "ensembl_gene_id", values = rownames(wu_counts), mart = mart)

#
# Assemble SCESet for Wu data
#
m <- match(rownames(wu_counts), bm$ensembl_gene_id)
data <- data.frame(bm[na.omit(m),])
rownames(data) <- data$ensembl_gene_id
colnames(data) <- c("gene_id", "symbol","biotype","description")

# Add neccessary metadata
wu_metadata$Cell <- make.names(wu_metadata$index)

featureData <- new("AnnotatedDataFrame", data = data[,c(2:ncol(data))])
phenoData <- new("AnnotatedDataFrame", data = wu_metadata)

wu <- newSCESet(countData = wu_counts, featureData = featureData, phenoData = phenoData)

#
# Assign cell cycle
#
pairs <- readRDS(system.file("exdata", parameters$normalisation$cycle_markers, package="scran"))
cellcycle <- cyclone(wu, pairs, gene.names = rownames(wu), assay="counts", verbose=F)

pData(wu)$cellcycle <- factor(cellcycle$phases)

#
# sumfactor normalisation
#
clusters <- quickCluster(wu, method = "hclust")

wu <- computeSumFactors(wu, clusters = clusters, sizes = seq(40, 115, 5))
wu <- normalize(wu, return_norm_as_exprs = F)

set_exprs(wu, "norm_exprs_sf") <- get_exprs(wu, "norm_exprs")

#
# TMM normalisation
#
wu <- normalizeExprs(wu, method="TMM", return_norm_as_exprs = FALSE)
set_exprs(wu, "norm_exprs_lineage") <- get_exprs(wu, "norm_exprs")

# Save
save(wu, file = file.path(parameters$general$path_supdata, "wu.Rdata"))
