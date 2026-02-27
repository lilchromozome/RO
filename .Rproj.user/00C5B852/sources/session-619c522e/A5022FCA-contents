library(ComplexHeatmap)
library(DESeq2)
library(circlize)
library(tidyverse)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(patchwork)
library(ggExtra)
# setwd("C:/Users/willllllli/Documents/Dr. Z lab/RNA seq/RO") # Windows
setwd("~/Documents/Zambidis lab/RNAseq/RO") # MAC

RO_full <- read.table("Counts/count.out",
                     header = T,
                     sep = "\t",
                     comment.char = "#",
                     stringsAsFactors = F)

RO_counts <- RO_full[, 7:ncol(RO_full)]
rownames(RO_counts) <- RO_full$Geneid
colnames(RO_counts) = c("H9_r1_E8", "H9_r2_E8",  "H9_r1_L3i", "H9_r2_L3i",
                            "RUES02_E8", "RUES02_L3i")
RO_counts <- RO_counts[, c("H9_r1_E8", "H9_r2_E8", "RUES02_E8", "H9_r1_L3i", "H9_r2_L3i", "RUES02_L3i")]
featureLength <- RO_full$Length

ens <- rownames(RO_counts)
symbols <- mapIds(org.Hs.eg.db,
                  keys = ens,
                  keytype = "ENSEMBL",
                  column = "SYMBOL",
                  multiVals = "first")
RO_counts$gene_symbol <- symbols
rownames(RO_counts) <- make.unique(ifelse(is.na(symbols), ens, symbols))
RO_counts$gene_symbol <- NULL

undiff_genes <- read.csv("undiff_dorgau.csv")
diff_genes <- read.csv("diff_dorgau.csv")

diff_cell_types <- unique(gsub("\\.(Liu|Dorgau)$", "", names(diff_genes)))

diff_genes_list <- lapply(diff_cell_types, function(ct) {
  matching_cols <- grep(paste0("^", ct, "\\.(Liu|Dorgau)$"), names(diff_genes), value = TRUE)
  unique(na.omit(unlist(diff_genes[matching_cols])))
})

names(diff_genes_list) <- diff_cell_types
max_len <- max(sapply(diff_genes_list, length))

diff_genes_df <- as.data.frame(
  lapply(diff_genes_list, function(x) {
    length(x) <- max_len  # pad with NA
    x
  })
)

# Differential expression ----------------------------------
samples <- colnames(RO_counts)
cellline <- sub("_(E8|L3i)", "", samples)
condition <- sub(".*_(E8|L3i).*", "\\1", samples)
coldata <- data.frame(
  cellline = factor(cellline),
  condition = factor(condition, levels = c("E8", "L3i")),
  row.names = samples
)
all(colnames(RO_counts) == rownames(coldata))
dds <-DESeqDataSetFromMatrix(
  countData = RO_counts,
  colData = coldata,
  design = ~ cellline + condition
)
dds <- dds[rowSums(counts(dds)) > 10, ]
dds <- DESeq(dds)
res <- results(dds, contrast = c("condition", "L3i", "E8"))

# rlog2
rlog <- cbind(
  E8 = -res$log2FoldChange,
  L3i = res$log2FoldChange
)
rownames(rlog) <- rownames(res)


## Z-score --------------------------------
mat <- counts(dds, normalized = TRUE)
mat_scaled <- t(scale(t(mat)))
mat_scaled <- mat_scaled[complete.cases(mat_scaled), , drop = FALSE]

# Heatmap plotting --------------------------------------
source("R scripts/visualization_functions.R")
result_undiff <- make_heatmap(mat_scaled, rlog, undiff_genes, "Undifferentiated Cell Types")
undiff_ht <- result_undiff$ht_combined
undiff_ht_rlog <- result_undiff$ht_combined_rlog

result_diff <- make_heatmap(mat_scaled, rlog, diff_genes, "Differentiated Cell Types")
diff_ht <- result_diff$ht_combined
diff_ht_rlog <- result_diff$ht_combined_rlog
