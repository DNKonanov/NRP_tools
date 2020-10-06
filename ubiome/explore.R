## Load other method data
setwd("~/git/nrps3/ubiome")
data <- read.table("all.pid.tsv", header=F, sep="\t")
colnames(data) <- c("query", "pid")

require(ggplot2)
ggplot(data, aes(x=pid)) + geom_histogram() + theme_bw()
