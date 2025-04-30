library(tidyverse)
library(ggrepel)

benchmarks <- read.delim(file = "workflow/stats/benchmarks_aggregated.tsv", sep = "\t", header = TRUE)

benchmarks <- benchmarks |> drop_na() |> filter(cpu_time>1) |> group_by(rule) |> summarise(cpu_time_mean = mean(cpu_time), cpu_time_stdev = sd(cpu_time), max_rss_mean = mean(max_rss), max_rss_stdev = sd(max_rss))

bench_data <- benchmarks |> filter(str_detect(rule, "qc"))

plot_qc <- ggplot(bench_data, aes(cpu_time_mean, max_rss_mean, label = rule)) + geom_point() + geom_text_repel() + labs(x="CPU time (s)", y="Memory usage (MB)")

bench_data <- benchmarks |> filter(str_detect(rule, "trim/"))

plot_trim <- ggplot(bench_data, aes(cpu_time_mean, max_rss_mean, label = rule)) + geom_point() + geom_text_repel() + labs(x="CPU time (s)", y="Memory usage (MB)")

bench_data <- benchmarks |> filter(str_detect(rule, "hisat|star"))

plot_align <- ggplot(bench_data, aes(cpu_time_mean, max_rss_mean, label = rule)) + geom_point() + geom_text_repel() + labs(x="CPU time (s)", y="Memory usage (MB)")

bench_data <- benchmarks |> filter(str_detect(rule, "expression|salmon|stringtie"))

plot_expr <- ggplot(bench_data, aes(cpu_time_mean, max_rss_mean, label = rule)) + geom_point() + geom_text_repel() + labs(x="CPU time (s)", y="Memory usage (MB)")