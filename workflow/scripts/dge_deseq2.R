# Differential Gene Expression with DESeq2
# Enrichment Analysis with clusterProfiler

#--------------------------
# Author: Núria Garriga Alonso
# GitHub: @nuriagaralon
# Repository: https://github.com/nuriagaralon/cel-rnaseq
#--------------------------

# Libraries
library(tidyverse)
library(ggpubr)
library(DESeq2)
library(pheatmap)
library(EnhancedVolcano)
library(clusterProfiler)
library(org.Ce.eg.db)
library(ggVennDiagram)

# Load data file
load("results/dge/gene_counts_data.RData")

### DESEQ2 MODEL AND ANALYSIS ###
# convert Generation and Replicate data to factors (instead of numeric)
meta_df$Generation <- factor(meta_df$Generation)
meta_df$Replicate <- factor(meta_df$Replicate)

# Create DESeq data object: dds
dds <- DESeqDataSetFromMatrix(
  countData = count_df,
  colData = meta_df,
  design = ~ Sample.Group
)

# Prefiltering (eliminate genes with low or zero counts)
bio_replicates <- meta_df |> 
  group_by(Sample.Group) |> 
  summarise(n = n()) |>  
  pull(n) |> min() 

smallest_groupsize <- bio_replicates #Biological replicates
keep <- rowSums(counts(dds) >= 10) >= smallest_groupsize
dds <- dds[keep, ]

# Releveling to set Control sample, good for default results.
dds$Sample.Group <- relevel(dds$Sample.Group, ref = "Negative_G1")

# Conduct DESeq analysis
dds <- DESeq(dds)

# VST transform
ddsvst <- vst(dds)

### EXPLORATORY DATA ANALYSIS ###
# PCA plot
pcaData <- plotPCA(ddsvst, intgroup = "Sample.Group", returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

pcaData34 <- plotPCA(ddsvst, intgroup = "Sample.Group", returnData = TRUE, pcsToUse = 3:4)
percentVar34 <- round(100 * attr(pcaData34, "percentVar"))

PCA_plot <- function(pcaData, percentVar,
                    pc_x = 1, pc_y = 2,
                    color_by = NULL, shape_by = NULL){
  ggplot(pcaData, aes(.data[[paste0("PC", pc_x)]],
                      .data[[paste0("PC", pc_y)]],
                      color = .data[[color_by]], label = .data[["name"]],
                      shape = .data[[shape_by]])) +
    geom_point(size = 3) +
    xlab(paste0("PC", pc_x, ": ", percentVar[1], "% variance")) +
    ylab(paste0("PC", pc_y, ": ", percentVar[2], "% variance"))
}

a <- PCA_plot(pcaData, percentVar, 1, 2, "Generation", "Type")
b <- PCA_plot(pcaData, percentVar, 1, 2, "Incubator", "Replicate")
c <- PCA_plot(pcaData34, percentVar34, 3, 4, "Generation", "Type")
d <- PCA_plot(pcaData34, percentVar34, 3, 4, "Incubator", "Replicate")

# Label Exp_2_3 and Neg_1_2 
pcaData$label <- ifelse(pcaData$name %in% c("Exp_2_3", "Neg_1_2"),
                        pcaData$name, NA)

a <- a + geom_text(aes(label = pcaData$label), color = "black", na.rm = TRUE, vjust = -1)

pca_all <- ggarrange(ggarrange(a, c, nrow = 2, common.legend = TRUE, legend = "right", labels = c("A", "C")),
                     ggarrange(b, d, nrow = 2, common.legend = TRUE, legend = "right", labels = c("B", "D")))

ggsave(filename = "results/dge/pca_label.pdf", plot = pca_all, width = 11, height = 7)

# Heatmap with pheatmap
metadata <- column_to_rownames(meta_df, "Sample.ID")

ddsvst_mat <- assay(ddsvst) #Get transformed matrix
ddsvst_cor <- cor(ddsvst_mat) #Get correlation values

heatmap_all <- pheatmap(ddsvst_cor,
                        annotation = dplyr::select(metadata, Generation, Type)
                        )$gtable

ggsave(filename = "results/dge/heatmap.pdf",
       plot = heatmap_all, width = 6, height = 5.5)

### DESEQ2 RESULTS ###
#Significant filters: padj < 0.05 and log2FC = 1 (Fold Change = 2)
alph <- 0.05
sig_lfc <- 1

# Contrast function
get_results <- function(dds_object, num, den, pval_a = 0.05,
                        lfc_t = 1, shrink = FALSE){
  res <- results(
    dds_object,
    contrast = c("Sample.Group", num, den),
    alpha = pval_a,
    lfcThreshold = lfc_t
  )
  
  if(shrink){
    res <- lfcShrink(dds_object, res = res, type = "ashr")
  }
  
  res_df <- as.data.frame(res) |>
    rownames_to_column("locus_tag") |>
    mutate(
      Contrast = paste(num, den, sep = "/")
    )
  list(
    res = res,
    res_df = res_df
  )
}

# All pairwise contrasts:
#sample_list <- combn(unique(meta_df$Sample.Group), 2) |>
#  t() |> as.data.frame() |> setNames(c("numerator", "denominator"))

# Interesting pairwise contrasts:
sample_list <- data.frame(
  numerator = c("Exp_G1","Exp_G1","Negative_G1",
                "Exp_G2","Exp_G2","Negative_G2",
                "Exp_G3","Exp_G3","Negative_G3",
                "Exp_G2","Exp_G3",
                "Negative_G2","Negative_G3",
                "Control_G2","Control_G3"),
  denominator = c("Negative_G1","Control_G1","Control_G1",
                  "Negative_G2","Control_G2","Control_G2",
                  "Negative_G3","Control_G3","Control_G3",
                  "Exp_G1","Exp_G1",
                  "Negative_G1","Negative_G1",
                  "Control_G1","Control_G1")
)

## Results, no shrinkage
#res_deseq_list <- pmap(
#  sample_list,
#  function(numerator, denominator){
#    get_results(
#      dds,
#      numerator,
#      denominator,
#      pval_a = alph,
#      lfc_t = sig_lfc,
#      shrink = FALSE
#    )
#  }
#)
#
#deseq_list <- map(res_deseq_list, "res")
#res_list <- map(res_deseq_list, "res_df")
#
#names(deseq_list) <- sapply(res_list, function(x){x$Contrast[1]})
#names(res_list) <- sapply(res_list, function(x){x$Contrast[1]})
#
#rm(res_deseq_list)
#
#res_all <- bind_rows(res_list) |> left_join(features_df, by = "locus_tag")
#
#res_all_filtered <- res_all |> 
#  filter(padj < alph) |> 
#  filter(abs(log2FoldChange) >= sig_lfc)
#
#write.csv(res_all_filtered, 
#          file = paste0("results/dge/res_all_filtered_", alph, "_", sig_lfc, ".csv"))

# Results, ashr shrinkage:
res_deseq_shrink_list <- pmap(
  sample_list,
  function(numerator, denominator){
    get_results(
      dds,
      numerator,
      denominator,
      pval_a = alph,
      lfc_t = sig_lfc,
      shrink = TRUE
    )
  }
)

deseq_shrink_list <- map(res_deseq_shrink_list, "res")
res_shrink_list <- map(res_deseq_shrink_list, "res_df")

names(deseq_shrink_list) <- sapply(res_shrink_list, function(x){x$Contrast[1]})
names(res_shrink_list) <- sapply(res_shrink_list, function(x){x$Contrast[1]})

rm(res_deseq_shrink_list)

res_shrink_all <- bind_rows(res_shrink_list) |> left_join(features_df, by = "locus_tag")

res_shrink_all_filtered <- res_shrink_all |> 
  filter(padj < alph) |> 
  filter(abs(log2FoldChange) >= sig_lfc)

write.csv(res_shrink_all_filtered, 
          file = paste0("results/dge/res_shrink_all_filtered_", alph, "_", sig_lfc, ".csv"))

### RESULTS: MA PLOT ###
contrast_list <- sample_list |>
  filter(grepl("Exp_G3", numerator) & grepl("Exp_G1", denominator) | 
         grepl("Negative_G3", numerator) & grepl("Negative_G1", denominator) | 
         grepl("Control_G3", numerator) & grepl("Control_G1", denominator))

contrast_list <- paste(contrast_list$numerator, contrast_list$denominator, sep = "/")

maplot_list <- map(deseq_shrink_list[contrast_list],
                   ~ ggmaplot(.x, fdr = 0.05, fc = 1, size = 0.2, top = 0))

# Save MA plots
maplots <- ggarrange(plotlist = maplot_list, nrow=1, labels="AUTO")
ggsave(filename = "results/dge/MAplotsG3G1.pdf", plot = maplots, height = 3.25, width = 13)

### RESULTS: VOLCANO PLOT ###
# This is an example. Edit contrast to your convenience
volc_contrast_list <- contrast_list

volc_plots <- map(volc_contrast_list, function(volc_contrast){
  volc_data <- res_shrink_all |> filter(Contrast == volc_contrast)
  
  top_padj <- volc_data |> filter(!is.na(padj)) |>
    arrange(padj) |> slice_head(n = 10) |> pull(symbol)
  
  top_lfc <- volc_data |> filter(!is.na(log2FoldChange)) |>
    arrange(desc(abs(log2FoldChange))) |> slice_head(n = 10) |> pull(symbol)
  
  volc_lab <- union(top_padj, top_lfc)
  
  EnhancedVolcano(volc_data,
                  lab = volc_data$symbol,
                  selectLab = volc_lab,
                  x = "log2FoldChange",
                  y = "padj",
                  pCutoff = 0.05,
                  FCcutoff = 1,
                  pointSize = 1,
                  titleLabSize = 0,
                  subtitleLabSize = 0,
                  subtitle = "",
                  title = "",
                  drawConnectors = TRUE,
                  max.overlaps = Inf
  )
})

volc_save <- ggarrange(plotlist = volc_plots,
                      labels = "AUTO", nrow = 1, common.legend = TRUE,
                      font.label = list(size = 20))

ggsave(filename = "results/dge/volcanos13.pdf",
       width = 17, height = 6.5, plot = volc_save)

### RESULTS: PLOT COUNTS ###
pvalcounts <- plotCounts(dds,
                         gene = "CELE_T07G12.5",
                         intgroup = "Generation",
                         returnData = TRUE) |> mutate(Gene = "svct-1")

lfccounts <- plotCounts(dds,
                        gene = "CELE_C55B7.4",
                        intgroup = "Generation",
                        returnData = TRUE) |> mutate(Gene = "acdh-1")

countsplot <- bind_rows(lfccounts, pvalcounts)

cplot <- ggplot(countsplot,
                aes(x = Generation,
                    y = count,
                    color = Gene)) +
  geom_point(size = 3,
             position = position_jitter(width = 0.1)) +
  scale_y_log10() +
  ylab("Counts")

ggsave(filename = "results/dge/plotcountsG1G3.pdf", plot = cplot, width = 7, height = 5)


### ENRICHMENT ANALYSIS ###
universe_list <- data.frame(locus_tag = rownames(dds)) |>
  left_join(features_df, by = "locus_tag") |>
  pull(GeneID)

enrich_contrast_list <- contrast_list

ego_list <- list()

for(enrich_contrast in enrich_contrast_list){
  enrich_data <- res_shrink_all_filtered |> filter(Contrast == enrich_contrast)
  
  # Enrichment analysis: upregulated
  enrich_up <- enrich_data |> filter(log2FoldChange > 0)
  gene_list_u <- enrich_up$GeneID
  
  ego_u <- enrichGO(gene = gene_list_u, 
                    universe = universe_list,
                    keyType = "ENTREZID", 
                    OrgDb = org.Ce.eg.db,
                    ont="BP",
                    pAdjustMethod = "BH",
                    pvalueCutoff = 0.05,
                    qvalueCutoff = 0.05
  )
  
  # Enrichment analysis: downregulated
  enrich_down <- enrich_data |> filter(log2FoldChange < 0)
  gene_list_d <- enrich_down$GeneID
  
  ego_d <- enrichGO(gene = gene_list_d, 
                    universe = universe_list,
                    keyType = "ENTREZID", 
                    OrgDb = org.Ce.eg.db,
                    ont="BP",
                    pAdjustMethod = "BH",
                    pvalueCutoff = 0.05,
                    qvalueCutoff = 0.05
  )
  
  # Store in list
  ego_list[[paste0("ego_u_", enrich_contrast)]] <- ego_u
  ego_list[[paste0("ego_d_", enrich_contrast)]] <- ego_d
}

ego_results <- map(ego_list, ~ as.data.frame(.x))

ego_dotplots <- map(ego_list, ~ dotplot(.x))

dots <- ggarrange(plotlist = ego_dotplots,
                  nrow = 3, ncol = 2, labels = "AUTO"
                  )

ggsave(filename = "results/dge/dotplots.pdf",
       width = 12, height = 12, plot = dots)

### RESULTS: VENN DIAGRAMS ###
# G3-G1 comparison #
venn_contrasts <- contrast_list

venn_type_G3_vs_G1 <- setNames(
  lapply(venn_contrasts, function(x) {
    res_shrink_all_filtered |>
    filter(Contrast == x) |>
    pull(GeneID)
  }),
  venn_contrasts
)

plot_venn_type_G3_vs_G1 <- ggVennDiagram(venn_type_G3_vs_G1, label = "count") +
  scale_fill_distiller(palette = "Blues", direction = 1) +
  scale_x_continuous(expand = expansion(mult = .22))

# G2-G1 comparison #
venn_contrasts_2 <- sample_list |>
  filter(grepl("Exp_G2", numerator) & grepl("Exp_G1", denominator) | 
         grepl("Negative_G2", numerator) & grepl("Negative_G1", denominator) | 
         grepl("Control_G2", numerator) & grepl("Control_G1", denominator))

venn_contrasts_2 <- paste(venn_contrasts_2$numerator,
                          venn_contrasts_2$denominator, sep = "/")

venn_type_G2_vs_G1 <- setNames(
  lapply(venn_contrasts_2, function(x) {
    res_shrink_all_filtered |>
      filter(Contrast == x) |>
      pull(GeneID)
  }),
  venn_contrasts_2
)

plot_venn_type_G2_vs_G1 <- ggVennDiagram(venn_type_G2_vs_G1, label = "count") +
  scale_fill_distiller(palette = "Blues", direction = 1) +
  scale_x_continuous(expand = expansion(mult = .22))

venns <- ggarrange(plot_venn_type_G2_vs_G1,
                   plot_venn_type_G3_vs_G1, labels = "AUTO")

ggsave(filename = "results/dge/venns.pdf", width = 12, height = 5)

# Enrichment analysis #
Exp_G3G1_only <- setdiff(
  venn_type_G3_vs_G1$`Exp_G3/Exp_G1`,
  union(venn_type_G3_vs_G1$`Negative_G3/Negative_G1`,
        venn_type_G3_vs_G1$`Control_G3/Control_G1`)
)

enrich_G3G1 <- enrichGO(gene = Exp_G3G1_only, 
                        universe = universe_list,
                        keyType = "ENTREZID", 
                        OrgDb = org.Ce.eg.db,
                        ont="BP",
                        pAdjustMethod = "BH",
                        pvalueCutoff = 0.05,
                        qvalueCutoff = 0.05
)


Exp_G2G1_only <- setdiff(
  venn_type_G2_vs_G1$`Exp_G2/Exp_G1`,
  union(venn_type_G2_vs_G1$`Negative_G2/Negative_G1`,
        venn_type_G2_vs_G1$`Control_G2/Control_G1`)
)

enrich_G2G1 <- enrichGO(gene = Exp_G2G1_only, 
                        universe = universe_list,
                        keyType = "ENTREZID", 
                        OrgDb = org.Ce.eg.db,
                        ont="BP",
                        pAdjustMethod = "BH",
                        pvalueCutoff = 0.05,
                        qvalueCutoff = 0.05
)
