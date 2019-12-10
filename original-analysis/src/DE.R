library(MAST)
library(ROCR)

differentialExpression <- function(expression, contrasts, fData) {
  diff_data <- lapply(names(contrasts), function(l) {
    cells <- contrasts[[l]]

    data <- cbind(expression[, cells[[1]]], expression[, cells[[2]]])
    cond_A <- unlist(strsplit(l, "_"))[1]
    cond_B <- unlist(strsplit(l, "_"))[2]

    cond <- c(rep(cond_A, length(cells[[1]])), rep(cond_B, length(cells[[2]])))

    cdat <- data.frame(wellKey = colnames(data), condition = factor(cond), stringsAsFactors = F)
    fdat <- data.frame(primerid = rownames(data), stringsAsFactors = F)
    sca <- MAST::FromMatrix(class = "SingleCellAssay",
                            exprsArray=data,
                            cData = cdat,
                            fData = fdat)

    zlm <- MAST::zlm(~ condition, sca, method = "bayesglm", ebayes = TRUE, ebayesControl = list(method = "MLE", model = "H1"))
    s <- MAST::summary(zlm, doLRT = paste0('condition', cond_B))$datatable
    res <- merge(s[contrast==paste0('condition', cond_B) & component=='H',.(primerid, `Pr(>Chisq)`)], #hurdle P values
                 s[contrast==paste0('condition', cond_B) & component=='logFC', .(primerid, coef, ci.hi, ci.lo)], by='primerid') #logFC coefficients
    colnames(res) <- c("geneID", "pval", "lfc", "lfc.hi", "lfc.lo")
    # Calculation of FDR
    res$fdr <- p.adjust(res$pval, method = "fdr")

    # Calculation of basemeanA expressions
    basemeanA <- rowMeans(data[, cond == cond_A])
    m <- match(res$geneID, names(basemeanA))
    res$basemeanA <- basemeanA[m]

    # Calculation of basemeanB expressions
    basemeanB <- rowMeans(data[, cond == cond_B])
    m <- match(res$geneID, names(basemeanB))
    res$basemeanB <- basemeanB[m]

    colnames(res[, c("basemeanA", "basemeanB")]) <- c(cond_A, cond_B)

    # Adding number of cells with detectable expression
    # n_exprsA <- rowSums(counts(subset[, cond == cond_A]) > 10, na.rm = T)
    # m <- match(res$geneID, names(n_exprsA))
    # res$n_exprsA <- n_exprsA[m]
    # n_exprsB <- rowSums(counts(subset[, cond == cond_B]) > 10, na.rm = T)
    # m <- match(res$geneID, names(n_exprsB))
    # res$n_exprsB <- n_exprsB[m]

    m <- match(res$geneID, rownames(fData))
    res <- cbind(res, fData[m, c("symbol", "biotype", "description")])

    res
  })
}

get_auroc <- function(gene, labels) {
  gene <- gene[!is.na(labels)]
  labels <- na.omit(labels)
  score <- rank(gene)
  # Get average score for each cluster
  ms <- aggregate(score ~ labels, FUN = mean)
  # Get cluster with highest average score
  posgroup <- ms[ms$score == max(ms$score), ]$labels
  # Return NAs if there is a tie for cluster with highest average score (by definition this is
  # not cluster specific)
  if (length(posgroup) > 1) {
    return(c(NA, NA, NA))
  }
  # Create 1/0 vector of truths for predictions, cluster with highest average score vs
  # everything else
  truth <- as.numeric(labels == posgroup)
  # Make predictions & get auc using RCOR package.
  pred <- prediction(score, truth)
  val <- unlist(ROCR::performance(pred, "auc")@y.values)
  pval <- suppressWarnings(wilcox.test(score[truth == 1], score[truth == 0])$p.value)
  return(c(val, posgroup, pval))
}

get_marker_genes <- function(dataset, labels) {
  res <- apply(dataset, 1, get_auroc, labels = labels)
  res <- data.frame(matrix(unlist(res), ncol = 3, byrow = T))
  colnames(res) <- c("auroc", "clusts", "pvalue")
  res$fdr <- p.adjust(res$pvalue)
  res$geneID <- rownames(dataset)
  return(res)
}
