#!/usr/local/bin/Rscript

# libraries
library(stringr)

# read in file
args = commandArgs(trailingOnly = TRUE)
saint = read.delim(args[1], sep="\t", header=T, as.is=T)

# set parameters
if (!is.null(args[2])) {
  fdr = as.numeric(args[2])
} else {
  fdr = 0.01
}

# get high confidence hits
saintFiltered = saint[saint$BFDR <= fdr, c('Bait', 'PreyGene', 'Spec')]

# get baits
baits = sort(unique(saintFiltered$Bait))

# iterate over baits
pdf('rsquared_plots.pdf')
par(mfrow=c(3,3))
rsquared = data.frame(
  bait = character(length(baits)),
  rsquared = numeric(length(baits)),
  adj_rsquared = numeric(length(baits)),
  stringsAsFactors = FALSE
)
for (i in 1:length(baits)) {
  splitReps = str_split_fixed(saintFiltered$Spec[saintFiltered$Bait == baits[i]], '\\|', 2)
  x = as.numeric(splitReps[ , 1])
  y = as.numeric(splitReps[ , 2])
  model = summary(lm(x~y))
  rsquared$bait[i] = baits[i]
  rsquared$rsquared[i] = round(model$r.squared, 3)
  rsquared$adj_rsquared[i] = round(model$adj.r.squared, 3)
  # plot
  plot(x, y, xlab = "rep 1", xlim=c(0, max(x)), ylab = "rep 2", ylim = c(0, max(y)))
  if (length(x) > 1) {
    abline(lm(y ~ x), lty = 2)
  }
  mtext(substitute(bait~R^2*equal*value, list(bait = baits[i], equal = '=', value = rsquared$rsquared[i])), side = 3)
}
dev.off()

write.table(rsquared, file="regression.txt", quote = FALSE, row.names = F)
