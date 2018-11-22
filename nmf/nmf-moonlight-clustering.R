#!/usr/local/bin/Rscript

library(gplots)
library(reshape2)

args <- commandArgs(trailingOnly = TRUE)
basis = read.csv(args[1], check.names=FALSE, header = TRUE, sep = ",", row.names = 1)
basis <- as.matrix(basis)

# get new heatmap
basis_ordered <- heatmap.2(basis)
basis_reordered = basis[rev(basis_ordered$rowInd), basis_ordered$colInd]

# print
write.table(basis_reordered, file = "clustered_moonlighting_basis.csv", sep = ",", quote = FALSE, col.names = NA)
