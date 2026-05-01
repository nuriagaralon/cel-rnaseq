
library(tidyverse)
library(tximport)

files <- list.files(path = "expression", pattern = "t_data.ctab", recursive = TRUE)

tmp <- read_tsv(files[1])
tx2gene <- tmp[, c("t_name", "gene_id")]
txi <- tximport(files, type = "stringtie", tx2gene = tx2gene)

# from tximport vignette
sampleTable <- data.frame(condition = factor(rep(c("A", "B"), each = 3)))
rownames(sampleTable) <- colnames(txi$counts)
dds <- DESeqDataSetFromTximport(txi, sampleTable, ~condition)
# dds is now ready for DESeq() see DESeq2 vignette

# https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#transcript-abundance-files-and-tximport-tximeta