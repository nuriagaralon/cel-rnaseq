# Generate benchmark plots from aggregated benchmark files
library(tidyverse)
library(ggrepel)

# Get data (Can be changed to benchmarks_aggregated_B)
benchmarks <- read.delim(file = "project_outputs/benchmarks_aggregated_A.tsv",
                         sep = "\t", header = TRUE)

benchmarks <- benchmarks |>
  drop_na() |>
  filter(cpu_time > 1) |>
  group_by(rule) |>
  summarise(cpu_time_mean = mean(cpu_time),
            cpu_time_stdev = sd(cpu_time),
            max_rss_mean = mean(max_rss),
            max_rss_stdev = sd(max_rss),
            .groups = "drop") |>
  mutate(
    category = case_when(
      str_detect(rule, "qc") ~ "QC",
      str_detect(rule, "trim") ~ "Trimming",
      str_detect(rule, "hisat|star") ~ "Alignment",
      str_detect(rule, "expression|salmon|stringtie|rsem") ~ "Expression"
    )
  )

# Sort facets
benchmarks$category <- factor(
  benchmarks$category,
  levels = c("QC", "Trimming", "Alignment", "Expression")
)

# Plot
ggplot(benchmarks, aes(cpu_time_mean, max_rss_mean, label = rule)) +
  geom_point(size = 2) +
  geom_errorbar(aes(xmin = cpu_time_mean - cpu_time_stdev,
                    xmax = cpu_time_mean + cpu_time_stdev), orientation = "y") +
  geom_errorbar(aes(ymin = max_rss_mean - max_rss_stdev,
                    ymax = max_rss_mean + max_rss_stdev)) +
  geom_text_repel(size = 3, max.overlaps = 50, box.padding = 0.5) +
  labs(x = "CPU time (s)", y = "Peak RSS (MB)") +
  facet_wrap(~ category, scales = "free") +
  theme_bw()
