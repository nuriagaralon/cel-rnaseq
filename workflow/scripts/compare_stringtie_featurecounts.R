library(tidyverse)

# Load files
fc1 <- read_tsv("Exp_2_3_gene_counts.tsv", comment = "#") |>
  select("Geneid", starts_with("results")) |>
  rename("FC_Exp_2_3" = starts_with("results"))

fc2 <- read_tsv("Neg_1_2_gene_counts.tsv", comment = "#") |>
  select("Geneid", starts_with("results")) |>
  rename("FC_Neg_1_2" = starts_with("results"))

st <- read_csv("gene_count_matrix.csv")

# structure featurecounts
fc <- full_join(fc1, fc2, by = "Geneid")

# Fix gene ID stringtie
st <- st |> mutate(Geneid = sub("\\|.*$", "", gene_id)) |>
  rename("ST_Neg_1_2" = Neg_1_2, "ST_Exp_2_3" = Exp_2_3 ) |>
  select(-gene_id)

# intersect 
common <- inner_join(st, fc, by = "Geneid")

only_fc <- anti_join(fc, st, by = "Geneid")
only_st <- anti_join(st, fc, by = "Geneid")


# Pearson correlation for sample Exp_2_3
cor.test(
  log2(common$FC_Exp_2_3 + 1),
  log2(common$ST_Exp_2_3 + 1),
  method = "pearson"
)

# Pearson correlation for sample Neg_1_2
cor.test(
  log2(common$FC_Neg_1_2 + 1),
  log2(common$ST_Neg_1_2 + 1),
  method = "pearson"
)
