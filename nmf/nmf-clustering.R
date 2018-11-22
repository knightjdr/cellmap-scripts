#!/usr/local/bin/Rscript

library(gplots)
library(reshape2)

args <- commandArgs(trailingOnly = TRUE)
basis = read.csv(args[1], check.names=FALSE, header = TRUE, sep = ",", row.names = 1)
basis <- as.matrix(basis)
scores = read.csv(args[2], check.names=FALSE, header = TRUE, sep = ",", row.names = 1)
scores <- as.matrix(scores)
if (length(args) < 3) {
  corrCutoff = 0.7
} else {
  corrCutoff = args[3]
}

# get new heatmap
basis_ordered <- heatmap.2(basis)
basis_reordered = basis[rev(basis_ordered$rowInd), basis_ordered$colInd]
scores_ordered <- heatmap.2(scores)
scores_reordered = scores[rev(scores_ordered$rowInd), scores_ordered$colInd]

# print
write.table(basis_reordered, file = "clustered_basis.csv", sep = ",", quote = FALSE, col.names = NA)
write.table(scores_reordered, file = "clustered_scores.csv", sep = ",", quote = FALSE, col.names = NA)

# correlation
cc = cor(t(basis), method='pearson')
# convert to lower triangular
cc[lower.tri(cc, diag=TRUE)] <- 0
# convert to dataframe
ccDf = setNames(melt(cc), c('source', 'target', 'edge'))
ccDfReduced = ccDf[ccDf$edge >= corrCutoff & ccDf$source != ccDf$target, ]
write.table(ccDfReduced, file = paste('cytoscape_nmf_correlation_', corrCutoff, '.tsv', sep=""), quote = FALSE, row.names = FALSE, sep = "\t")
