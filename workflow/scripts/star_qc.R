# Aggregate STAR quality logs

#--------------------------
# Author: Núria Garriga Alonso
# GitHub: @nuriagaralon
# Repository: https://github.com/nuriagaralon/cel-rnaseq
#--------------------------

# log file
log_file <- snakemake@log[[1]]

log_con <- file(log_file, open = "wt")
sink(log_con)
sink(log_con, type = "message")

# Libraries
library(dplyr)
library(readr)
library(purrr)
library(stringr)
library(tidyr)

# Get files
align_stats_files <- snakemake@input[["stats"]]

align_stats <- align_stats_files |>
  map(function(file){
    df <- read.delim(file, sep = "|", 
                     header = FALSE, strip.white = TRUE)
    filename <- basename(file) |> str_remove("_Log.final.out")
    names(df) <- c("Stats", filename)
    df
  })

# Pivot table to get stats per column
align_df <- purrr::reduce(align_stats, full_join, by = "Stats") |> 
  pivot_longer(
    cols = -Stats,
    names_to = "Sample",
    values_to = "value") |>
  pivot_wider(
    names_from = Stats,
    values_from = value) |>
  mutate(across(-Sample, ~ as.numeric(gsub("%", "", .x))))

# Calculate % of interesting stats
align_perc <- align_df |>
  transmute(
    Sample,
    Total_reads = `Number of input reads`,
    Unique_pct = `Uniquely mapped reads number` / Total_reads * 100,
    MultiMapping_pct = (`Number of reads mapped to multiple loci` +
                        `Number of reads mapped to too many loci`) / Total_reads * 100,
    Unmapped_pct = (`Number of reads unmapped: too many mismatches` +
                      `Number of reads unmapped: too short` +
                      `Number of reads unmapped: other`) / Total_reads * 100
  ) 

align_summ <- align_perc |> 
  summarise_if(is.numeric, list(mean = mean,sd = sd, min = min, max = max)) |> 
  pivot_longer(everything(), names_to = "Sample", values_to = "value") |> 
  extract( Sample, into = c("Stats", "Sample"), regex = "^(.*)_([^_]*)$") |>
  pivot_wider(names_from = Stats, values_from = value)

align_res <- rbind(align_perc, align_summ) |>
  mutate(across(-Sample, ~ round(.x, 2)))

# Save file
write_tsv(align_res, file = snakemake@output[["stardata"]])

# Close log
sink(type = "message")
sink()
close(log_con)