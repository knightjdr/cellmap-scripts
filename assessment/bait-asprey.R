#!/usr/local/bin/Rscript

# libraries
library(cba)
library(gplots)
library(plyr)
library(reshape2)

# variables
fdr = 0.01

# arguments
args <- commandArgs(trailingOnly = TRUE)
baitMap <- read.delim(args[1], sep="\t", header=T, as.is=T)
saint <- read.delim(args[2], sep="\t", header=T, as.is=T)

# control subtract
temp = lapply(lapply(strsplit(saint$ctrlCounts, "\\|"), as.numeric), mean)
saint$AvgSpec <- round(saint$AvgSpec - unlist(temp), 2)
saint$AvgSpec[saint$AvgSpec < 0] <- 0
saint <- saint[saint$AvgSpec > 0, ]

# get baits from SAINT file and create new map from it
baits <- unique(saint$Bait)
newBaitMap <- baitMap[baitMap$bait %in% baits, ]
newBaitMap <- newBaitMap[order(newBaitMap$bait), ]

# count reciprocal interactions
saintReciprocal <- saint[saint$PreyGene %in% newBaitMap$gene & saint$BFDR <= fdr, c('Bait', 'PreyGene', 'BFDR')]
saintReciprocal <- setNames(saintReciprocal, c('Bait', 'Prey', 'BFDR'))
saintReciprocal$Bait <- mapvalues(saintReciprocal$Bait, from = newBaitMap$bait, to = newBaitMap$gene)
saintReciprocal <- unique(saintReciprocal[ , 1:2])
matchedPairs <- c()
baitBaitTotal = 0
recipMatches = 0
recipTotal = 0
recipBaitList <- unique(saintReciprocal$Bait)
for(i in 1:length(recipBaitList)) {
  matchedPreys <- saintReciprocal[saintReciprocal$Bait == recipBaitList[i], c('Prey')]
  # for each prey as a bait, check if its preys match current bait
  for(j in 1:length(matchedPreys)) {
    baitBaitTotal = baitBaitTotal + 1
    if (
      !(paste(recipBaitList[i], matchedPreys[j], sep='-') %in% matchedPairs)
    ) {
      recipTotal = recipTotal + 1
      reciprocalBait <- which(matchedPreys[j] == recipBaitList)[1]
      recipMatchedPreys <- saintReciprocal[saintReciprocal$Bait == recipBaitList[reciprocalBait], c('Prey')]
      if (recipBaitList[i] %in% recipMatchedPreys) {
        recipMatches = recipMatches + 1
        matchedPairs <- c(matchedPairs, paste(recipBaitList[reciprocalBait], recipBaitList[i], sep='-'))
      }
    }
  }
}
print(paste('Reciprocal matches:', recipMatches, 'of', recipTotal, sep=' '))

# grab all rows for preys that are also baits and reduce dataframe, and cluster
saintBaitsAsPrey <- saint[saint$PreyGene %in% newBaitMap$gene, c('Bait', 'PreyGene', 'AvgSpec', 'BFDR')]
saintBaitsAsPrey <- setNames(saintBaitsAsPrey, c('Bait', 'Prey', 'AvgSpec', 'BFDR'))

# add zeros when a bait is not found as a prey
for(i in 1:length(newBaitMap$bait)) {
  for(j in 1:length(newBaitMap$gene)) {
    if (!(newBaitMap$gene[j] %in% saintBaitsAsPrey[saintBaitsAsPrey$Bait == newBaitMap$bait[i], c('Prey')])) {
      saintBaitsAsPrey <- rbind(saintBaitsAsPrey, c(newBaitMap$bait[i], newBaitMap$gene[j], 0, 0))
    }
  }
}
saintBaitsAsPrey$AvgSpec <- as.numeric(saintBaitsAsPrey$AvgSpec)
saintBaitsAsPrey$BFDR <- as.numeric(saintBaitsAsPrey$BFDR)

# order data frame by bait name then by bait gene name
#saintBaitsAsPrey <- saintBaitsAsPrey[order(saintBaitsAsPrey$Bait, saintBaitsAsPrey$Prey), ]

# convert data frame to matrix for heatmap
saintBaitsAsPrey.matrix <- acast(saintBaitsAsPrey, Bait ~ Prey, value.var="AvgSpec", fun.aggregate = mean)

# perform clustering and optimize
dist_bait <- dist(saintBaitsAsPrey.matrix, method = 'euclidean')
dist_prey<- dist(t(saintBaitsAsPrey.matrix), method = 'euclidean')
hc_bait <- hclust(dist_bait, method = 'complete')
hc_prey <- hclust(dist_prey, method = 'complete')
opt_order_bait <- order.optimal(dist_bait, hc_bait$merge)
opt_order_prey <- order.optimal(dist_prey, hc_prey$merge)
saintBaitsAsPrey.matrix <- saintBaitsAsPrey.matrix[opt_order_bait$order, opt_order_prey$order, drop = FALSE]
bait_order <- factor(row.names(saintBaitsAsPrey.matrix), levels = row.names(saintBaitsAsPrey.matrix))
prey_order <- factor(names(saintBaitsAsPrey.matrix[1, ]), levels = names(saintBaitsAsPrey.matrix[1, ]))

# create heatmap and get order from it
#saintBaitsAsPrey_ordered <- heatmap.2(saintBaitsAsPrey.matrix)
#bait_order <- row.names(saintBaitsAsPrey.matrix)[rev(saintBaitsAsPrey_ordered$rowInd)]
#prey_order <- colnames(saintBaitsAsPrey.matrix)[saintBaitsAsPrey_ordered$colInd]

# convert matrix back to dataframe
outputDf <- saintBaitsAsPrey[ , c('Prey', 'Bait', 'AvgSpec', 'BFDR')]
outputDf <- setNames(outputDf, c('row', 'column', 'value', 'score'))
# re-order
outputDf$row <- factor(outputDf$row, levels = prey_order)
outputDf$column <- factor(outputDf$column, levels = bait_order)
outputDf <- outputDf[order(outputDf$row, outputDf$column), ]
# output
outputDf$params <- ""
outputDf$params[1] <- paste('{"type": "dotplot", "kind": "Bait gene vs Bait", "xAxis": "Bait gene", "yAxis": "Bait", "filterType": 1, "primary": 0.95, "secondary": 0.9, "score": "BFDR", "abundance": "AvgSpec"}', sep="")
write.table(outputDf, file = paste('baits-asprey.tsv', sep=""), quote = FALSE, row.names = FALSE, sep = "\t")
