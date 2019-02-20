#!/usr/local/bin/Rscript

library('reshape2')

#read in file
args = commandArgs(trailingOnly = TRUE)
data = read.delim(args[1], sep="\t", header=T, as.is=T)
if (length(args) < 2) {
  fdr = 0.01
} else {
  fdr = args[2]
}

#preys to omit from analysis (highlight abundant frequent flyers)
omit = c()

# filter data by preys
sigPreys = unique(data$PreyGene[data$BFDR <= fdr & !(data$PreyGene %in% omit)])
fData = data[data$PreyGene %in% sigPreys, ]

# control subtraction
temp = lapply(lapply(strsplit(fData$ctrlCounts, "\\|"), as.numeric), mean)
fData$AvgSpec <- round(fData$AvgSpec - unlist(temp), 2)

fData <- fData[order(fData$Bait),]

#convert dataframe to matrix
scData = acast(fData, Bait ~ PreyGene, value.var="AvgSpec", fun.aggregate = mean)
scData[is.na(scData)] = 0
scData[scData < 0] = 0

#normalize each prey column to 1
for(i in 1:ncol(scData)) {
  scData[, i] = scData[, i] / max(scData[, i])
}

#export
write.csv(scData, file="sc-matrix.csv", quote = FALSE)
