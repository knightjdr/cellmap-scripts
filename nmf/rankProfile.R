#!/usr/local/bin/Rscript

args <- commandArgs(trailingOnly = TRUE)
basis = read.table(args[1], sep=",", header = TRUE, row.names = 1)
scores = read.table(args[2], sep=",", header = TRUE, row.names = 1)
goCat = args[3] # one of GO:CC, GO:BP or GO:MF
if (length(args) < 4) {
  resultsFolder = './';
} else {
  resultsFolder = paste(args[4], "/", sep="");
}

#libraries
library(gProfileR)
library(openxlsx)

#variables
ranks = ncol(basis)
genes = nrow(basis)
baits = nrow(scores)
maxGenes = 100
minValue = 0.25 # a prey must have an NMF value at or above this to be used for enrichment
percentageMax = 0.75 # if a prey has an NMF value within this % of max, it can be used in those ranks

#filter basis so that a gene only occurs in its top ranks
fData = basis
for(i in 1:genes) {
  currMax = max(fData[i, ])
  for(j in 1:ranks) {
    if(
      fData[i,j] < percentageMax * currMax |
      fData[i,j] < minValue
    ) {
      fData[i,j] = 0
    }
  }
}

#get genes that define each rank
rankGenes = list()
for(i in 1:ranks) {
  nonZeros = sum(fData[,i] > 0)
  #for top maxgenes
  if(nonZeros > maxGenes) {
    nonZeros = maxGenes
  }
  rankGenes[[i]] = rownames(fData)[order(fData[,i], decreasing=TRUE)][1:nonZeros]
}

# set background
background = rownames(fData)

#gProfiler and output rank data to summary.xlsx
profiles = list()
wb = createWorkbook()
for(i in 1:ranks) {
  print(i)
  profiles[[i]] = gprofiler(
    rankGenes[[i]],
    organism = "hsapiens",
    ordered_query = F,
    significant = T,
    exclude_iea = F,
    underrep = F,
    evcodes = F,
    region_query = F,
    max_p_value = 0.01,
    min_set_size = 0,
    max_set_size = 0,
    min_isect_size = 0,
    correction_method = "gSCS",
    hier_filtering = "none",
    domain_size = "annotated",
    custom_bg = background,
    numeric_ns = "",
    png_fn = NULL,
    include_graph = F,
    src_filter = c(goCat)
  )
  profiles[[i]] = profiles[[i]][order(profiles[[i]]$p.value), ]
  addWorksheet(wb, paste("rank", i))
  if(nrow(profiles[[i]]['term.name']) > 0) {
    writeData(wb, i, profiles[[i]])
  }
}
saveWorkbook(wb, paste(resultsFolder, "summary.xlsx", sep=""), overwrite = TRUE)

#cluster ranks
distRank <- dist(as.matrix(t(basis)), method = 'euclidean')
hcRank <- hclust(distRank, method = 'complete')
write(paste(hcRank$order, collapse = ", "), file = paste(resultsFolder, "rank_order.txt", sep=""))

#create text file listing all go terms for each rank
cat('rank\tname\tmatched\tpvalue\tid\tgenes\n', file = paste(resultsFolder, "terms_perrank.txt", sep=""))
for(i in 1:ranks) {
  for(j in 1:nrow(profiles[[i]])) {
    cat(i, '\t', sep="",  file = paste(resultsFolder, "terms_perrank.txt", sep=""), append = TRUE)
    cat(profiles[[i]]$term.name[j], '\t', sep="", file = paste(resultsFolder, "terms_perrank.txt", sep=""), append = TRUE)
    cat(profiles[[i]]$overlap.size[j], '\t', sep="", file = paste(resultsFolder, "terms_perrank.txt", sep=""), append = TRUE)
    cat(profiles[[i]]$p.value[j], '\t', sep="", file = paste(resultsFolder, "terms_perrank.txt", sep=""), append = TRUE)
    cat(profiles[[i]]$term.id[j], '\t', sep="", file = paste(resultsFolder, "terms_perrank.txt", sep=""), append = TRUE)
    cat(profiles[[i]]$intersection[j], '\n', sep="", file = paste(resultsFolder, "terms_perrank.txt", sep=""), append = TRUE)
  }
}

#output each preys top rank
cat('gene\trank\tscore\n', file = paste(resultsFolder, "gene-localizations.txt", sep=""))
for(i in 1:genes) {
  cat(rownames(basis)[i], which(basis[i, ] == max(basis[i, ])), max(basis[i, ]), '\n', sep="\t", file = paste(resultsFolder, "gene-localizations.txt", sep=""), append = TRUE)
}
