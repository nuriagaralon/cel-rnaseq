# Aggregate featureCounts quality logs

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
expr_stats_files <- snakemake@input[["stats"]]

expr_stats <- expr_stats_files |>
  map(function(file){
    df <- read_tsv(file)
    names(df)[2] <- names(df)[2] |>
      basename() |>
      str_remove("_Aligned.sortedByCoord.out.bam")
    df
  })

# Pivot table to get stats per column
expr_df <- purrr::reduce(expr_stats, full_join, by = "Status") |> 
  pivot_longer(
    cols = -Status,
    names_to = "Sample",
    values_to = "value") |>
  pivot_wider(
    names_from = Status,
    values_from = value)

# Calculate % of interesting stats
expr_perc <- expr_df |>
  transmute(
    Sample,
    Total_reads = Assigned + Unassigned_Unmapped +
      Unassigned_MultiMapping + Unassigned_NoFeatures + Unassigned_Ambiguity,
    Assigned_pct = Assigned / Total_reads * 100,
    Unmapped_pct = Unassigned_Unmapped / Total_reads * 100,
    MultiMapping_pct = Unassigned_MultiMapping / Total_reads * 100,
    NoFeatures_pct = Unassigned_NoFeatures / Total_reads * 100,
    Ambiguity_pct = Unassigned_Ambiguity / Total_reads * 100
  )

expr_summ <- expr_perc |> 
  summarise_if(is.numeric, list(mean = mean,sd = sd, min = min, max = max)) |> 
  pivot_longer(everything(), names_to = "Sample", values_to = "value") |> 
  extract( Sample, into = c("Stats", "Sample"), regex = "^(.*)_([^_]*)$") |>
  pivot_wider(names_from = Stats, values_from = value)

expr_res <- rbind(expr_perc, expr_summ) |>
  mutate(across(-Sample, ~ round(.x, 2)))

# Save file
write_tsv(expr_res, file = snakemake@output[["fcdata"]])

# Close log
sink(type = "message")
sink()
close(log_con)
