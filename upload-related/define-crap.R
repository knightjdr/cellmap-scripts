#!/usr/local/bin/Rscript

args = commandArgs(trailingOnly = TRUE)

# get data
saint = read.delim(args[1], sep="\t", header=T, as.is=T)

# parameters
fdr = 0.01 # the cutoff for significant interactors
baitTheshold = 20 # the number of bait a prey must be seen with to be considered "crap"

# filter SAINT by AvgP cutoff
saintFiltered = saint[saint$BFDR <= fdr, ]

# control subtraction
temp = lapply(lapply(strsplit(saintFiltered$ctrlCounts, "\\|"), as.numeric), mean)
saintFiltered$AvgSpec = round(saintFiltered$AvgSpec - unlist(temp), 2)
saintFiltered$AvgSpec[saintFiltered$AvgSpec < 0] = 0

# get prey list
preys = unique(saintFiltered$PreyGene)

# define crap
crap = c()
crapDF <- data.frame(gene=character(), occurence=integer(), median=double(), stringsAsFactors = FALSE)
for(i in 1:length(preys)) {
  median = median(saintFiltered$AvgSpec[saintFiltered$PreyGene == preys[i]])
  occurence = length(saintFiltered$PreyGene[saintFiltered$PreyGene == preys[i]])
  if (occurence >= baitTheshold) {
    crap = c(crap, preys[i])
  }
  newrow = c(preys[i], occurence)
  crapDF = rbind(crapDF, data.frame("gene" = preys[i], "occurence" = occurence, "median" = median))
}
crap = sort(crap)
crapDF = crapDF[with(crapDF, order(-occurence)), ]

# write to file
write(crap, file = "crap-list.txt")
write.table(crapDF, file = "crap-details.txt", quote = FALSE, row.names = FALSE, sep = "\t")
