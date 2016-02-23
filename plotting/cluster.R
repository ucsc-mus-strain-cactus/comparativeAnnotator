args <- commandArgs(trailingOnly=T)

data <- read.csv(args[1], header=T, row.names=1)
title <- args[2]
title <- gsub("_", " ", title)
mat <- sapply(as.data.frame(data), as.logical)
library(pvclust)
fit <- pvclust(mat, method.hclust="ward", method.dist="binary")
pdf(paste(args[3], ".pdf", sep=""))
plot(fit, main=title, xlab="Binary clustering (Ward's Method)", ylab="Distance")
dev.off()