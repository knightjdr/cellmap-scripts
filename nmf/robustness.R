library(ggplot2)
library(reshape2)

data = read.table("robustness.txt", sep="\t", header = TRUE)
columns = c(0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1)
columnNames = c("90", "80", "70", "60", "50", "40", "30", "20", "10")
ranks = 20

pdf('boxplot.pdf')
par(mfrow=c(3,3))
for (i in 1:20) {
  rankData = data[data$rank == i, ]
  rankData$percentile = factor(rankData$percentile, columns)
  title = paste("Rank", i, sep=" ")
  boxplot(
    RBD~percentile,
    data = rankData,
    outline = FALSE,
    las = 2,
    xlab = "Percentage of Genes (%)",
    ylab = "RBD",
    main = title,
    xaxt='n',
    yaxt='n'
  )
  axis(1, at = seq(1,9,1), labels = columnNames)
  axis(2, at = seq(0,1,0.5), labels = c("0", "0.5", "1"))
}
dev.off()
