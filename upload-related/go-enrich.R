#!/usr/local/bin/Rscript

# Input a SAINT file and it will generate a list of enriched GO terms for each bait
# based on preys that pass a specifided FDR. Also input a map of bait names to cell-map database
# IDs so that output JSON will be named via this ID.

# libraries
library(gProfileR)
library(jsonlite)

# parse command line
args <- commandArgs(trailingOnly = TRUE)
saint = read.delim(args[1], sep="\t", header=T, as.is=T)
map = read.delim(args[2], sep="\t", header=T, as.is=T)
fdr = 0.01
goCat = "GO:CC"

# filter SAINT so that only significant preys remain
fData = saint[saint$BFDR <= fdr, ]
baits = unique(fData$Bait)

# find custom background
preyBackground = unique(fData$PreyGene)

# go enrichment and JSON output
cat('[\n', file = "go-enrichment.json")
for(i in 1:length(baits)) {
  print(i)
  id = as.numeric(map$id[map$name == baits[i]])
  cat('\t{\n\t\t"_id": ', id, ',\n', sep="", file = "go-enrichment.json", append = TRUE)
  genes = fData$PreyGene[fData$Bait == baits[i]]
  terms = gprofiler(
    genes,
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
    custom_bg = preyBackground,
    numeric_ns = "",
    png_fn = NULL,
    include_graph = F,
    src_filter = c(goCat)
  )
  terms = terms[order(terms$p.value),]
  drop = c("query.number", "significant", "domain", "subgraph.number", "recall", "precision")
  terms = terms[ , !(names(terms) %in% drop)]
  newNames = c("pvalue", 'termsize', 'querysize', 'overlap', "id", "term", "depth", "intersection")
  colnames(terms) = newNames
  cat('\t\t"details": ', toJSON(terms), '\n', sep="", file = "go-enrichment.json", append = TRUE)
  if(i < length(baits)) {
    cat('\t},\n', file = "go-enrichment.json", append = TRUE)
  } else {
    cat('\t}\n', file = "go-enrichment.json", append = TRUE)
  }
}
cat(']', file = "go-enrichment.json", append = TRUE)
