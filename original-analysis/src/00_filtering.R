library(scater)
library(Rtsne)
library(multipanelfigure)
library(viridis)
library(dplyr)
library(tidyr)
library(gridExtra)
library(singlecellutils)
library(genefilter)
library(RColorBrewer)

#
# Load data
#
source("src/parameters.R")
load(file.path(parameters$general$path_rdata, "sceData.RData"))

#
# Create map of QC data
#
qcdata <- pData(sceData)[,parameters$filtering$qc_columns]
qcdata$platform <- factor(pData(sceData)$Platform)
qcdata$background <- factor(pData(sceData)$Background)

if (!is.null(parameters$filtering$log_qc_columns)) {
  qcdata[, parameters$filtering$log_qc_columns] <- log2(qcdata[, parameters$filtering$log_qc_columns] + 1)
}

is_c1 <- pData(sceData)$Platform == "C1"

set.seed(5411)
t_c1 <- Rtsne(scale(qcdata[is_c1, parameters$filtering$qc_columns]), perplexity = 40, max_iter = 2000)
set.seed(5411)
t_wg <- Rtsne(scale(qcdata[!is_c1, parameters$filtering$qc_columns]), perplexity = 40, max_iter = 2000)

# Add coordinates to qcdata
qcdata$tsne.x <- NA
qcdata$tsne.y <- NA
qcdata$tsne.x[is_c1] <- t_c1$Y[,1]
qcdata$tsne.x[!is_c1] <- t_wg$Y[,1]
qcdata$tsne.y[is_c1] <- t_c1$Y[,2]
qcdata$tsne.y[!is_c1] <- t_wg$Y[,2]

theme1 <- theme(plot.background = element_blank(),panel.grid.major = element_line(size=.2, colour = "grey"),panel.grid.minor = element_line(size=.1, colour = "grey"),
                panel.border = element_blank(),panel.background = element_blank(),axis.line.x = element_line(size=.3),axis.line.y = element_line(size=.3), legend.text = element_text(size=8),
                axis.title.x = element_blank(),axis.title.y = element_blank(), plot.title = element_text(face="bold", color="black", size=8),legend.key.size =  unit(2, "mm"), legend.margin=unit(-25, "mm"),
                strip.background = element_blank(), strip.text = element_text(face="bold", color="black", size=6, margin = margin(0,0,-1,0)), axis.text = element_text(color = "black", size = 8))

m <- multi_panel_figure(width = 205, height = 230, columns = 2, rows = 5)
for(c in parameters$filtering$qc_columns) {
  f <- ggplot(qcdata, aes_string(x="tsne.x", y="tsne.y", color=eval(c))) +
    geom_point(size=0.5) +
    scale_color_viridis(name="") +
    facet_wrap(~platform, ncol=2) +
    xlab("") + ylab("") +
    ggtitle(c) +
    theme1
  m <- fill_panel(m, f)
}

#m

ggsave(
  plot = m,
  filename = file.path(parameters$general$path_rfigures, "Supplementary-Cell-QC-map.png"),
  width = 205,
  height = 230,
  units = "mm",
  dpi = 600)

#
# Filter cells based on defined criteria
#
# We construct a binary QC matrix, with TRUE meaning that the cell failed for this QC measure.
binary.qc.l <- lapply(names(parameters$filtering$cell_outlier_columns), function(c) {
  res <- rep(NA, nrow(qcdata))
  c1o <- isOutlier(qcdata[is_c1, c], nmads = as.numeric(parameters$filtering$cell_outlier_columns[[c]][1]), type=parameters$filtering$cell_outlier_columns[[c]][2], log=as.logical(parameters$filtering$cell_outlier_columns[[c]][3]))
  wgo <- isOutlier(qcdata[!is_c1, c], nmads = as.numeric(parameters$filtering$cell_outlier_columns[[c]][1]), type=parameters$filtering$cell_outlier_columns[[c]][2], log=as.logical(parameters$filtering$cell_outlier_columns[[c]][3]))
  res[!is_c1] <- wgo
  res[is_c1] <- c1o
  return(res)
})

binary.qc <- do.call("cbind", binary.qc.l)
colnames(binary.qc) <- names(parameters$filtering$cell_outlier_columns)
rownames(binary.qc) <- rownames(qcdata)

# Determine failed cells as those with more than 1 failed criterium
failed.cells <- (rowSums(binary.qc) > parameters$filtering$cell_outlier_maxfailed)
qcdata$qc.exclude <- failed.cells

write.table(x=qcdata, file=file.path(parameters$general$path_supdata, "QC-metrics.txt"), sep="\t", quote = F, row.names = T, col.names = T)
pData(sceData)$qc.exclude <- failed.cells

#
# Supplementary Figure 1: Cell QC
#
theme2 <- theme(plot.background = element_blank(),panel.grid.major = element_line(size=.2, colour = "grey"),panel.grid.minor = element_line(size=.1, colour = "grey"),
                panel.border = element_blank(), panel.background = element_blank(), axis.line.x = element_line(size=.3), axis.line.y = element_line(size=.3), legend.text = element_text(size=6),
                axis.title.x = element_blank(), axis.title.y = element_blank(), plot.title = element_text(face="bold", color="black", size=6), legend.key.size =  unit(2, "mm"), legend.margin=unit(-25, "mm"),
                strip.background = element_blank(), strip.text = element_blank(), axis.text = element_text(color="black", size=8) )

xintercept_colors <- c("#8c510a", "#bf812d","#01665e", "#35978f")
xintercept_platforms <- c("C1", "WG")
get_min <- function(c) {
  c(min(qcdata[!binary.qc[,c] & is_c1,c]), min(qcdata[!binary.qc[,c] & !is_c1,c]))
}
get_max <- function(c) {
  c(max(qcdata[!binary.qc[,c] & is_c1,c]), max(qcdata[!binary.qc[,c] & !is_c1,c]))
}

# Mitochondrial content
c <- 'pct_counts_feature_controls_mito'
xintercept_mins <- get_min(c)
xintercept_maxs <- get_max(c)
intercept_data <- data.frame(platform = xintercept_platforms, min = xintercept_mins, max = xintercept_maxs)
s1_mito <- ggplot(qcdata, aes_string(x=eval(c))) +
  geom_histogram(binwidth = 1) +
  geom_vline(data = intercept_data, aes(xintercept=min), color="#8c510a") +
  geom_vline(data = intercept_data, aes(xintercept=max), color="#01665e") +
  facet_wrap(~platform, ncol=2) +
  xlab("") + ylab("") + guides(color=FALSE) +
  xlim(0,15) +
  #ggtitle(c) +
  theme2

# Total detected features
c <- 'total_features'
xintercept_mins <- get_min(c)
xintercept_maxs <- get_max(c)
intercept_data <- data.frame(platform = xintercept_platforms, min = xintercept_mins, max = xintercept_maxs)
s1_features <- ggplot(qcdata, aes_string(x=eval(c))) +
  geom_histogram(binwidth = 0.05) +
  geom_vline(data = intercept_data, aes(xintercept=min), color="#8c510a") +
  geom_vline(data = intercept_data, aes(xintercept=max), color="#01665e") +
  facet_wrap(~platform, ncol=2) +
  xlab("") + ylab("") + guides(color=FALSE) +
  xlim(6,8) +
  #ggtitle(c) +
  theme1

# Dropout
c <- 'pct_dropout'
xintercept_mins <- get_min(c)
xintercept_maxs <- get_max(c)
intercept_data <- data.frame(platform = xintercept_platforms, min = xintercept_mins, max = xintercept_maxs)
s1_dropout <- ggplot(qcdata, aes_string(x=eval(c))) +
  geom_histogram(binwidth = 0.01) +
  geom_vline(data = intercept_data, aes(xintercept=min), color="#8c510a") +
  geom_vline(data = intercept_data, aes(xintercept=max), color="#01665e") +
  facet_wrap(~platform, ncol=2) +
  xlab("") + ylab("") + guides(color=FALSE) +
  #ggtitle(c) +
  theme2

# Percent of genes
c <- 'percent_genes'
xintercept_mins <- get_min(c)
xintercept_maxs <- get_max(c)
intercept_data <- data.frame(platform = xintercept_platforms, min = xintercept_mins, max = xintercept_maxs)
s1_genes <- ggplot(qcdata, aes_string(x=eval(c))) +
  geom_histogram(binwidth = 1) +
  geom_vline(data = intercept_data, aes(xintercept=min), color="#8c510a") +
  geom_vline(data = intercept_data, aes(xintercept=max), color="#01665e") +
  facet_wrap(~platform, ncol=2) +
  xlab("") + ylab("") + guides(color=FALSE) +
  #ggtitle(c) +
  theme2

# Expression of Rplp0
c <- 'exprs_feature_controls_housekeeping'
xintercept_mins <- get_min(c)
xintercept_maxs <- get_max(c)
intercept_data <- data.frame(platform = xintercept_platforms, min = xintercept_mins, max = xintercept_maxs)
s1_rplp <- ggplot(qcdata, aes_string(x=eval(c))) +
  geom_histogram(binwidth = 0.1) +
  geom_vline(data = intercept_data, aes(xintercept=min), color="#8c510a") +
  geom_vline(data = intercept_data, aes(xintercept=max), color="#01665e") +
  facet_wrap(~platform, ncol=2) +
  xlab("") + ylab("") + guides(color=FALSE) +
  xlim(2.5, 4) +
  #ggtitle(c) +
  theme2

# QC map with features
s1_qc_features <- ggplot(qcdata, aes_string(x="tsne.x", y="tsne.y", color=eval(c))) +
  geom_point(size=0.5) +
  scale_color_viridis(name="") +
  facet_wrap(~platform, ncol=2) +
  xlab("") + ylab("") +
  #ggtitle("Total features (log2 scale)") +
  theme1 +
  theme(legend.position="bottom")

s1_qc_fail <- ggplot(qcdata, aes_string(x="tsne.x", y="tsne.y", color="qc.exclude")) +
  geom_point(size=1) +
  scale_colour_manual(values = c("black", "red"), name="") +
  facet_wrap(~platform, ncol=2) +
  xlab("") + ylab("") +
  #ggtitle("QC failed cells") +
  guides(color=FALSE) +
  theme2

write.table(file = file.path(parameters$general$path_supdata, "Source_Data_SupplementaryFigure_1ae.txt"),
            x = qcdata[, c("platform", "total_features", "pct_dropout", "pct_counts_feature_controls_mito", "percent_genes", "exprs_feature_controls_housekeeping", "qc.exclude")],
            quote = F,
            row.names = T,
            col.names = T,
            sep = "\t")

qcdata %>%
  select(Lineage = background, platform, qc.exclude) %>%
  group_by(Lineage, platform, qc.exclude) %>%
  summarise(cells = n()) %>%
  unite("key", platform, qc.exclude, sep=".") %>%
  spread(key, cells) %>%
  mutate(WG.TRUE = ifelse(is.na(WG.TRUE), 0, WG.TRUE), WG.FALSE = ifelse(is.na(WG.FALSE), 0, WG.FALSE)) %>%
  mutate(C1.sum = C1.TRUE + C1.FALSE, WG.sum = WG.TRUE + WG.FALSE, raw.sum = C1.sum + WG.sum, clean.sum = C1.FALSE + WG.FALSE) %>%
  mutate(Raw = paste0(raw.sum, " (", C1.sum, "/", WG.sum, ")")) %>%
  mutate(Clean = paste0(clean.sum, " (", C1.FALSE, "/", WG.FALSE, ")")) %>%
  select(Lineage, Raw, Clean) -> filter.table.raw

tt1 <- ttheme_minimal(core=list(fg_params=list(hjust=1, x=0.9)),
                      rowhead=list(fg_params=list(hjust=1, x=0.9)),
                      colhead = list(fg_params=list(hjust=1, x=0.9)))

s1_table <- tableGrob(filter.table.raw, theme=tt1, rows = NULL)

m <- multi_panel_figure(width = 205, height = 230, columns = 2, rows = 5)
m <- fill_panel(m, s1_features, row = 1, col = 1)
m <- fill_panel(m, s1_dropout, row = 2, col = 1)
m <- fill_panel(m, s1_mito, row = 3, col = 1)
m <- fill_panel(m, s1_genes, row = 4, col = 1)
m <- fill_panel(m, s1_rplp, row = 5, col = 1)
m <- fill_panel(m, s1_qc_features, row = c(1,2), column = 2)
m <- fill_panel(m, s1_qc_fail, row = c(3,4), column = 2)
m <- fill_panel(m, s1_table, row = 5, column = 2)

ggsave(
  plot = m,
  filename = file.path(parameters$general$path_rfigures, "Supplementary_Figure_1.pdf"),
  width = 205,
  height = 230,
  units = "mm",
  dpi = 600)

#
# Gene filtering
#
isl1.lineage <- pData(sceData)$Background == "isl1"
nkx.lineage <- pData(sceData)$Background == "nkx2-5"

# Build gene filters
aggregate.expr <- EpOverA(A = parameters$filtering$EpOverA_A)
k.cell.expr <- kOverA(k=parameters$filtering$kOverA_k, A=parameters$filtering$kOverA_A)

# Make sure both platform counts satisfies the criteria
ffun <- filterfun(aggregate.expr, k.cell.expr)
keep_wg_crit <- genefilter(counts(sceData[,(pData(sceData)$is.sc & !pData(sceData)$qc.exclude & !is_c1)]), ffun)
keep_c1_crit <- genefilter(counts(sceData[,(pData(sceData)$is.sc & !pData(sceData)$qc.exclude & is_c1)]), ffun)

keep_isl1_crit <- genefilter(counts(sceData[,(pData(sceData)$is.sc & !pData(sceData)$qc.exclude & is_c1 & isl1.lineage)]), ffun)
keep_nkx_crit <- genefilter(counts(sceData[,(pData(sceData)$is.sc & !pData(sceData)$qc.exclude & is_c1 & nkx.lineage)]), ffun)

# Make sure that the corresponding count is not 0
one <- kOverA(k=1,A=0)

# Make sure that the wafergene count is also not 0
keep_wg <- genefilter(counts(sceData[,(pData(sceData)$is.sc & !pData(sceData)$qc.exclude & !is_c1)]), one)

# Make sure that the C1 count is also not 0
keep_c1 <- genefilter(counts(sceData[,(pData(sceData)$is.sc & !pData(sceData)$qc.exclude & is_c1)]), one)

keep <- (keep_isl1_crit | keep_nkx_crit) & keep_c1 & keep_wg
sum(keep)
scd <- sceData[keep,(pData(sceData)$is.sc & !pData(sceData)$qc.exclude)]
scd <- calculateQCMetrics(scd)

#
# Save object for further downstream analysis
#
save(scd, file = file.path(parameters$general$path_supdata, "scd.Rdata"))
