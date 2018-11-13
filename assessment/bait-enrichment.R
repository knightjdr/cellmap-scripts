#!/usr/local/bin/Rscript

# needed for xlsx memory
options(java.parameters = "-Xmx16000m")

#libraries
library(gProfileR)
library(openxlsx)

# parameters
fdr = 0.01
srcFilter = c('GO:BP', 'GO:CC', 'GO:MF', 'CORUM')
topPreyNumber = 0

args = commandArgs(trailingOnly = TRUE)
saint = read.delim(args[1], sep="\t", header=T, as.is=T)

# set parameters
if (!is.null(args[2])) {
  fdr = as.numeric(args[2])
} 
if (!is.null(args[3])) {
  topPreyNumber = as.numeric(args[3])
}

# control subtraction
temp = lapply(lapply(strsplit(saint$ctrlCounts, "\\|"), as.numeric), mean)
saint$AvgSpec <- round(saint$AvgSpec - unlist(temp), 2)

# get prey background
preyBackground = unique(saint$PreyGene[saint$BFDR <= fdr])

# prey length normalization.
normValue = median(saint$PreySequenceLength)
saint$AvgSpec <- round(saint$AvgSpec * (normValue / saint$PreySequenceLength), 2)

# baits to profile
baits <- unique(saint$Bait)

for(i in 1:length(baits)) {
  print(paste(i, baits[i], sep=" "))
  # get rows of SAINT file for current bait
  baitSaint = saint[
    saint$Bait == baits[i] &
      saint$BFDR <= fdr,
    ]
  # check if bait has less than topPreyNumber preys
  currPreyNumber = nrow(baitSaint)
  #for top maxgenes
  if(
    topPreyNumber > 0 &
    currPreyNumber > topPreyNumber
  ) {
    currPreyNumber = topPreyNumber
  }
  currPreys = baitSaint$PreyGene[order(baitSaint$AvgSpec / baitSaint$PreySequenceLength, decreasing=TRUE)][1:currPreyNumber]
  currProfile = gprofiler(
    currPreys,
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
    src_filter = srcFilter
  )
  if (nrow(currProfile) > 0) {
    # drop some uneeded columns
    currProfile <- currProfile[ , !(names(currProfile) %in% c('query.number', 'significant'))]
    # add column with bait name and move to start
    currProfile$bait <- baits[i]
    currProfile <- currProfile[ , c(ncol(currProfile), 1:(ncol(currProfile) - 1))]
    # add data to data frames for output
    if (exists('allProfile') && is.data.frame(get('allProfile'))) {
      allProfile = rbind(allProfile, currProfile)
      bpProfile = rbind(bpProfile, currProfile[currProfile$domain == 'BP', ])
      ccProfile = rbind(ccProfile, currProfile[currProfile$domain == 'CC', ])
      corumProfile = rbind(corumProfile, currProfile[currProfile$domain == 'cor', ])
      mfProfile = rbind(mfProfile, currProfile[currProfile$domain == 'MF', ])
    } else {
      allProfile = currProfile
      bpProfile = currProfile[currProfile$domain == 'BP', ]
      ccProfile = currProfile[currProfile$domain == 'CC', ]
      corumProfile = currProfile[currProfile$domain == 'cor', ]
      mfProfile = currProfile[currProfile$domain == 'MF', ]
    }
  }
}
# sort
allProfile = allProfile[order(allProfile$bait, allProfile$p.value), ]
bpProfile = bpProfile[order(bpProfile$bait, bpProfile$p.value), ]
ccProfile = ccProfile[order(ccProfile$bait, ccProfile$p.value), ]
corumProfile = corumProfile[order(corumProfile$bait, corumProfile$p.value), ]
mfProfile = mfProfile[order(mfProfile$bait, mfProfile$p.value), ]

# output results
if (topPreyNumber > 0) {
  outfilename = paste('bait-enrichment_top', topPreyNumber, '.xlsx', sep='')
} else {
  outfilename = "bait-enrichment.xlsx"
}
wb = createWorkbook()
addWorksheet(wb, 'all')
addWorksheet(wb, 'BP')
addWorksheet(wb, 'CC')
addWorksheet(wb, 'MF')
addWorksheet(wb, 'Corum')
writeData(wb, 1, allProfile, rowNames = FALSE)
writeData(wb, 2, bpProfile, rowNames = FALSE)
writeData(wb, 3, ccProfile, rowNames = FALSE)
writeData(wb, 4, mfProfile, rowNames = FALSE)
writeData(wb, 5, corumProfile, rowNames = FALSE)
saveWorkbook(wb, outfilename, overwrite = TRUE)
