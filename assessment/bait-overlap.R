#!/usr/local/bin/Rscript

# libraries
library(cba)
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
if (!is.null(args[4])) {
  comparisonType = args[4]
} else {
  comparisonType = "jaccard"
}
if (!is.null(args[5])) {
  corrMethod = args[5]
} else {
  corrMethod = "pearson"
}

distance = function(type, method) {
  if (type == 'jaccard') {
    jaccard = function(x, y) {
      return(length(intersect(x, y)) / length(union(x, y)))
    }
    return(jaccard)
  } else {
    if (method == 'kendall') {
      statistic = 'tau'
    } else if (method == 'pearson') {
      statistic = 'cor'
    } else {
      statistic = 'rho'
    }
    correlation = function(x, y) {
      df = rbind(sourcePreys, targetPreys)
      if (length(unique(df$Bait)) == 1) {
        return(1)
      }
      df = dcast(df, PreyGene ~ Bait, value.var = 'AvgSpec')
      rownames(df) = df[,1]
      df = df[,-1]
      df[is.na(df)] = 0
      result = cor.test(df[,1], df[,2], method = method, exact = FALSE)
      return(result$estimate[[statistic]])
    }
    return(correlation)
  }
}
extract = function(type) {
  if (type == 'jaccard') {
    jaccard = function(df, bait) {
      return(df$PreyGene[df$Bait == bait])
    }
    return(jaccard)
  } else {
    correlation = function(df, bait) {
      return(df[df$Bait == bait, c('Bait', 'PreyGene', 'AvgSpec')])
    }
    return(correlation)
  }
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

comparisonFunc = distance(comparisonType, corrMethod)
extractFunc = extract(comparisonType)
rowIndex = 1
for(i in 1:length(baits)) {
  print(paste(i, baits[i], sep=" "))

  # get preys for source bait
  sourcePreys = extractFunc(saint, baits[i])
  
  for(j in 1:length(baits)) {
    # get rows of SAINT file for target bait
    targetPreys = extractFunc(saint, baits[j])
    
    distances$source[rowIndex] = baits[i]
    distances$target[rowIndex] = baits[j]
    distances$distance[rowIndex] = comparisonFunc(sourcePreys, targetPreys)
    rowIndex = rowIndex + 1
  }
}

# output table
outfilePrefix = paste('bait-overlap-', comparisonType, sep='')
if (comparisonType == "correlation") {
  outfilePrefix = paste(outfilePrefix, '-', corrMethod, sep='')
} 
if (topPreyNumber > 0) {
  outfilePrefix = paste(outfilePrefix, '-top', topPreyNumber, sep='')
}
outfilename = paste(outfilePrefix, '.txt', sep='')
write.table(distances, file = outfilename, quote = FALSE, row.names = F, sep = "\t")

# calculate distance and cluster with optimal ordering
matrix = acast(distances, source~target, value.var="distance")
if (comparisonType == "jaccard") {
  distMatrix = as.dist((matrix - 1) * -1)
} else {
  distMatrix = dist(matrix, method = "euclidean")
}
clust = hclust(distMatrix, method = "complete")
optimal = order.optimal(distMatrix, clust$merge)
optimizedMatrix = matrix[optimal$order, optimal$order, drop = FALSE]

# draw heatmap
fontsize = 0.15
if (comparisonType == "jaccard") {
  palette = colorRampPalette(c("#ffffff", "#0040ff", "#000000"))(101)
} else {
  palette = colorRampPalette(c("#0000ff", "#ffffff", "#ff0000"))(101)
}
outfilename = paste(outfilePrefix, '.pdf', sep='')
pdf(outfilename, width = 6, height = 6)
heatmap.2(
  optimizedMatrix,
  dendrogram = 'none',
  cexCol = fontsize,
  cexRow = fontsize,
  col = palette,
  Rowv = FALSE,
  Colv = FALSE,
  tracecol = NA
)
dev.off()
