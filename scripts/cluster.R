args <- commandArgs(trailingOnly=T)

data <- read.csv(args[1], header=T, row.names=1)
base_title <- args[2]
genome <- args[3]
num_ok <- args[4]
percent_ok <- args[5]
gencode <- args[6]
biotype <- args[7]

base_title <- gsub("_", " ", base_title)
title <- paste(base_title, sprintf("\ngenome: %s    %s (%s%%) not OK transcripts\nGencode set: %s    Biotype: %s", genome, num_ok, percent_ok, gencode, biotype))


mat <- sapply(as.data.frame(data), as.logical)
library(pvclust)
fit <- pvclust(mat, method.hclust="ward", method.dist="binary")
pdf(args[8])
plot(fit, main=title, xlab="Binary clustering (Ward's Method)", ylab="Distance")
dev.off()