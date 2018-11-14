#!/usr/local/bin/Rscript

# script takes a SAFE node_properties_annotation-highest.txt file and finds enriched GO
# terms per domain

args <- commandArgs(trailingOnly = TRUE)
nodeProperties = read.table(args[1], sep="\t", header = TRUE)
goCat = args[2] # one of GO:CC, GO:BP or GO:MF
if (length(args) < 3) {
  resultsFolder = './';
} else {
  resultsFolder = paste(args[3], "/", sep="");
}

#libraries
library(gProfileR)
library(openxlsx)
library(stringr)

#variables
noDomains = length(unique(nodeProperties$Domain))

# set background
background = as.vector(nodeProperties$Node.label)

#gProfiler and output rank data to summary.xlsx
profiles = list()
wb = createWorkbook()
for(i in 2:noDomains) {
  print(i)
  currGenes = as.vector(nodeProperties$Node.label[nodeProperties$Domain..predominant. == i])
  profiles[[i - 1]] = gprofiler(
    currGenes,
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
  profiles[[i - 1]] = profiles[[i - 1]][order(profiles[[i - 1]]$p.value), ]
  addWorksheet(wb, paste("domain", i))
  if(nrow(profiles[[i - 1]]['term.name']) > 0) {
    writeData(wb, i - 1, profiles[[i - 1]])
  }
}
saveWorkbook(wb, paste(resultsFolder, "domain-summary.xlsx", sep=""), overwrite = TRUE)
