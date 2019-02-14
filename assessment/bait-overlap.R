#!/usr/local/bin/Rscript

# libraries
library(gplots)
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
if (!is.null(args[3])) {
  topPreyNumber = as.numeric(args[3])
} else {
  topPreyNumber = 0
}

# jaccard distance
jaccard = function(x, y) {
  return(length(intersect(x, y)) / length(union(x, y)))
}

# filter dataset
saint = saint[saint$BFDR <= fdr,]

# unique baits
baits = unique(saint$Bait)

# filter by top preys
if (topPreyNumber > 0) {
  # control subtraction
  temp = lapply(lapply(strsplit(saint$ctrlCounts, "\\|"), as.numeric), mean)
  saint$AvgSpec <- round(saint$AvgSpec - unlist(temp), 2)

  saintTop = saint[FALSE,]
  for(i in 1:length(baits)) {
    # get preys for source bait
    baitDF = saint[saint$Bait == baits[i], ]
    
    # prey length normalization
    normValue = median(baitDF$PreySequenceLength)
    baitDF$AvgSpec <- round(baitDF$AvgSpec * (normValue / baitDF$PreySequenceLength), 2)
    
    # filter
    currPreyNumber = nrow(baitDF)
    if(currPreyNumber > topPreyNumber) {
      currPreyNumber = topPreyNumber
    }
    baitDF = baitDF[order(-baitDF$AvgSpec), ]
    baitDF = baitDF[1:currPreyNumber,]
    
    saintTop = rbind(saintTop, baitDF)
  }
  saint = saintTop
}

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
if (topPreyNumber > 0) {
  outfilePrefix = paste('bait-overlap_top', topPreyNumber, sep='')
} else {
  outfilePrefix = "bait-overlap"
}
outfilename = paste(outfilePrefix, '.txt', sep='')
write.table(distances, file = outfilename, quote = FALSE, row.names = F, sep = "\t")

# distance matrix
matrix = acast(distances, source~target, value.var="distance")
distFunc = function(x) dist(x, method = "euclidean")
hclustFunc <- function(x) hclust(x, method = "complete")

# draw heatmap
fontsize = 0.15
palette = colorRampPalette(c("#ffffff", "#0040ff", "#000000"))(101)
outfilename = paste(outfilePrefix, '.pdf', sep='')
pdf(outfilename, width = 6, height = 6)
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
