#!/usr/local/bin/Rscript

# libraries
library(gplot)
library(RColorBrewer)
library(reshape2)

args = commandArgs(trailingOnly = TRUE)
saint = read.delim(args[1], sep="\t", header=T, as.is=T)

# set parameters
if (!is.null(args[2])) {
  fdr = as.numeric(args[2])
} else {
  fdr = 0.01
}

# jaccard distance
jaccard = function(x, y) {
  return(length(intersect(x, y)) / length(union(x, y)))
}

# filter dataset
saint = saint[saint$BFDR <= fdr,]

# unique baits
baits = unique(saint$Bait)

# initialize dataframe for storing distances
distances = data.frame(matrix(ncol = 3, nrow = length(baits) * length(baits)))
colnames(distances) = c("source", "target", "distance")

rowIndex = 1
for(i in 1:length(baits)) {
  print(paste(i, baits[i], sep=" "))

  # get preys for source bait
  sourcePreys = saint$PreyGene[saint$Bait == baits[i]]
  
  for(j in 1:length(baits)) {
    # get rows of SAINT file for target bait
    targetPreys = saint$PreyGene[saint$Bait == baits[j]]
    
    distances$source[rowIndex] = baits[i]
    distances$target[rowIndex] = baits[j]
    distances$distance[rowIndex] = jaccard(sourcePreys, targetPreys)
    rowIndex = rowIndex + 1
  }
}

# output table
write.table(distances, file="bait-overlap.txt", quote = FALSE, row.names = F, sep = "\t")

# distance matrix
matrix = acast(distances, source~target, value.var="distance")
distFunc = function(x) dist(x, method = "euclidean")
hclustFunc <- function(x) hclust(x, method = "complete")

# draw heatmap
fontsize = 0.15
palette = colorRampPalette(c("#ffffff", "#0040ff", "#000000"))(101)
pdf("bait-overlap.pdf", width = 6, height = 6)
heatmap.2(
  matrix,
  cexCol = fontsize,
  cexRow = fontsize,
  col = palette,
  distfun = distFunc,
  hclustfun = hclustFunc,
  tracecol = NA
)
dev.off()
