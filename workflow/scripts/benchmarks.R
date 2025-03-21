library(ggplot2)
library(ggrepel)

setwd()

benchmarks <- read.delim(file = "20250314-131457_benchmarks.tsv", sep = "\t", header = TRUE)

plot <- ggplot(benchmarks, aes(cpu_time, max_rss, label = rule)) + geom_point() + geom_text_repel() + labs(x="CPU time (s)", y="Memory usage (MB)")