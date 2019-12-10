library(scater)
library(ggplot2)
library(multipanelfigure)

#
# Load data
#
source("src/parameters.R")
load(file.path(parameters$general$path_supdata, "wu.Rdata"))

#
# Smooth muscle marker genes
#
timepoint <- which(pData(wu)$timepoint == "e9.5")
sm_marker_genes <- c("ENSMUSG00000015579", "ENSMUSG00000032085", "ENSMUSG00000001349", "ENSMUSG00000035783", "ENSMUSG00000029761", "ENSMUSG00000022836", "ENSMUSG00000048878", "ENSMUSG00000045667")

subset <- wu[sm_marker_genes, timepoint]
rownames(subset) <- fData(subset)$symbol

#plotExpression(subset, x = "type", features = 1:nrow(subset), colour_by = "batch", exprs_values = "norm_exprs_sf", log2_values = F)

#
# SM cell inference
#
marker_expression <- get_exprs(subset, "norm_exprs_sf")

is_smc_cl <- list(
  as.numeric(marker_expression["Nkx2-5",]) < 1,
  as.numeric(marker_expression["Tagln",]) > 2,
  as.numeric(marker_expression["Cnn1",]) > 2,
  as.numeric(marker_expression["Acta2",]) > 2,
  as.numeric(marker_expression["Cald1",]) > 2,
  as.numeric(marker_expression["Mylk",]) > 2,
  #as.numeric(marker_expression["Hexim1",]) > 2,
  as.numeric(marker_expression["Smtnl2",]) > 2
)

is_smc_c <- do.call("rbind", is_smc_cl)
is_smc <- (colSums(is_smc_c[2:nrow(is_smc_c), ]) >= 5) & is_smc_c[1, ]

t <- table(isko = pData(subset)$batch == "KO", is_smc)[c(2,1), c(2,1)]

tr <- prop.test(t, correct=T, conf.level = 0.99)

ggMMplot <- function(var1, var2){
  require(ggplot2)
  levVar1 <- length(levels(var1))
  levVar2 <- length(levels(var2))

  jointTable <- prop.table(table(var1, var2))
  plotData <- as.data.frame(jointTable)
  plotData$marginVar1 <- prop.table(table(var1))
  plotData$var2Height <- plotData$Freq / plotData$marginVar1
  plotData$var1Center <- c(0, cumsum(plotData$marginVar1)[1:levVar1 -1]) +
    plotData$marginVar1 / 2

  ggplot(plotData, aes(var1Center, var2Height)) +
    geom_bar(stat = "identity", aes(fill = var2, width = marginVar1), col = "White") +
    geom_label(aes(label = as.character(var1), x = var1Center, y = 1.05), label.r = unit(0, "lines"))
}

plot_data <- data.frame(Celltype = factor(ifelse(is_smc, "SMC", "non-SMC"), levels = c("SMC", "non-SMC")), Genotype = factor(ifelse(pData(subset)$batch == "KO", "KO", "WT"), levels = c("WT", "KO")))

plot <- ggMMplot(plot_data$Genotype, plot_data$Celltype) +
  guides(fill = F) +
  scale_fill_manual(values = c("#A8DBA8", "#3B8686")) +
  ylab("SMC Proportion") +
  xlab("Cells") +
  theme(panel.background = element_rect(fill = "transparent"), axis.line.y = element_line(colour = "black"), axis.line.x = element_line(colour = "black"))#, axis.text.x = element_blank(), axis.title.x = element_blank(), axis.ticks.x = element_blank())

m <- multi_panel_figure(width = 205, height = 230, columns = 2, rows = 5)
m <- fill_panel(m, plot, row = c(1,2))
ggsave(
  plot = m,
  filename = file.path(parameters$general$path_rfigures, "Supplementary_Figure_SMCprop.pdf"),
  width = 205,
  height = 230,
  units = "mm",
  dpi = 600)
