# Aggregate expression data

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
library(tibble)
library(stringr)

# Feature table mapping
features_df <- read_tsv(snakemake@input[['features']])

features_df <- features_df |>
  select(symbol, GeneID, locus_tag) |>
  distinct() |>
  mutate(across(everything(), as.character))


# Build counts dataframe
files <- snakemake@input[["counts"]]

count_tables <- files |>
  map(function(file){
    df <- read_tsv(file, comment = "#") |>
            select(Geneid, starts_with("results"))
    
    names(df)[2] <- names(df)[2] |>
      basename() |>
      str_remove("_Aligned.sortedByCoord.out.bam")
    
    df
  })

count_df <- purrr::reduce(count_tables, full_join, by = "Geneid")
count_df <- column_to_rownames(count_df, "Geneid")

# Build metadata dataframe
meta_df <- read_tsv(snakemake@input[['metadata']])

# Sample list
count_samples <- names(count_df)
metadata_samples <- meta_df$Sample.ID

# SAFETY CHECK: SAMPLES AND NAMES
if(!setequal(count_samples, metadata_samples)){
  stop(
    paste0(
      "Samples do not match between counts matrix and metadata.\n",
      "Missing in metadata: ",
      paste(setdiff(count_samples, metadata_samples), collapse = ", "), "\n",
      "Missing in counts: ",
      paste(setdiff(metadata_samples, count_samples), collapse = ", ")      
    )
  )
}

# Reorder metadata
meta_df <- meta_df[match(count_samples, metadata_samples), ]

# SAFETY CHECK: SAMPLE ORDER
if(!identical(count_samples, meta_df$Sample.ID)){
  stop("Metadata sample ordering failed after matching")
}

# Save objects
save(
  list = c("count_df", "meta_df", "features_df"), 
  file = snakemake@output[["rdata"]]
)

# Close log
sink(type = "message")
sink()
close(log_con)