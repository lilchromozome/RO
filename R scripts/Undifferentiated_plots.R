library(ComplexHeatmap)
library(DESeq2)
library(circlize)
library(tidyverse)
library(patchwork)
library(ggExtra)

undiff <- read.delim("Counts/RO_counts_undifferentiated_gene_symbol.csv", sep=",")
rownames(undiff) <- undiff$X
undiff$X <- NULL
undiff <- undiff[, c("E5C3_E8", "H9_E8", "RUES02_E8", "E5C3_L3i","H9_L3i",  "RUES02_L3i")]

RO <- read.delim("Counts/RO_counts_gene_symbol.csv", sep=",")
rownames(RO) <- RO$X
RO$X <- NULL
colnames(RO) <- c("H9_E8_r1_RO", "H9_E8_r2_RO", "H9_L3i_r1_RO",  "H9_L3i_r2_RO", "RUES02_E8_RO", "RUES02_L3i")
RO <- RO[, c("H9_E8_r1_RO", "H9_E8_r2_RO", "RUES02_E8_RO", "H9_L3i_r1_RO",  "H9_L3i_r2_RO",  "RUES02_L3i")]

all <- cbind(undiff,RO)

undiff_genes <- read.csv("retinal_progenitor_genes.csv")
diff_genes <- read.csv("diff_dorgau.csv")

# Undiff Differential expression ----------------------------------
samples <- colnames(undiff)
cellline <- sub("_(E8|L3i)", "", samples)
condition <- sub(".*_(E8|L3i).*", "\\1", samples)
coldata <- data.frame(
  cellline = factor(cellline),
  condition = factor(condition, levels = c("E8", "L3i")),
  row.names = samples
)
all(colnames(undiff) == rownames(coldata))
dds_undiff <-DESeqDataSetFromMatrix(
  countData = undiff,
  colData = coldata,
  design = ~ cellline + condition
)
dds_undiff <- dds_undiff[rowSums(counts(dds_undiff)) > 10, ]
dds_undiff <- DESeq(dds_undiff)
undiff_res <- results(dds_undiff, contrast = c("condition", "L3i", "E8"))

# rlog2
undiff_rlog <- cbind(
  E8 = -undiff_res$log2FoldChange,
  L3i = undiff_res$log2FoldChange
)
rownames(undiff_rlog) <- rownames(undiff_res)
# Z-score 
undiff_mat <- counts(dds_undiff, normalized = TRUE)
undiff_mat_scaled <- t(scale(t(undiff_mat)))
undiff_mat_scaled <- undiff_mat_scaled[complete.cases(undiff_mat_scaled), , drop = FALSE]


## RO differential expression -------------------------
samples <- colnames(RO)
cellline <- sub("_(E8|L3i)", "", samples)
condition <- sub(".*_(E8|L3i).*", "\\1", samples)
coldata <- data.frame(
  cellline = factor(cellline),
  condition = factor(condition, levels = c("E8", "L3i")),
  row.names = samples
)
all(colnames(RO) == rownames(coldata))
dds_RO <-DESeqDataSetFromMatrix(
  countData = RO,
  colData = coldata,
  design = ~ cellline + condition
)
dds_RO <- dds_RO[rowSums(counts(dds_RO)) > 10, ]
dds_RO <- DESeq(dds_RO)
RO_res <- results(dds_RO, contrast = c("condition", "L3i", "E8"))

# rlog2
RO_rlog <- cbind(
  E8 = -RO_res$log2FoldChange,
  L3i = RO_res$log2FoldChange
)
rownames(RO_rlog) <- rownames(RO_res)
# Z-score 
RO_mat <- counts(dds_RO, normalized = TRUE)
RO_mat_scaled <- t(scale(t(RO_mat)))
RO_mat_scaled <- RO_mat_scaled[complete.cases(RO_mat_scaled), , drop = FALSE]

## Visualizations ----------------------------------
source('R scripts/visualization_functions.R')

make_heatmap(undiff_mat_scaled, undiff_rlog, undiff_genes, 
             origin = "Undifferentiated", title="Stem-Progenitor Cell Types")
make_heatmap(undiff_mat_scaled, undiff_rlog, diff_genes, 
             origin = "Undifferentiated", title="Differentiated Cell Types")

make_heatmap(RO_mat_scaled, RO_rlog, undiff_genes, origin = "RO", title="Stem-Progenitor Cell Types")
make_heatmap(RO_mat_scaled, RO_rlog, diff_genes, origin = "RO", title="Differentiated Cell Types")

cpm <- t(t(all / colSums(all))) * 1e6
log_mat <- log10(cpm+1)
mat_scaled <- t(scale(t(cpm)))
mat_scaled[is.nan(mat_scaled)] <- 0
make_heatmap_undiffAndROCounts(mat_scaled, undiff_genes, title = "Stem-Progenitor Cell Types")
make_heatmap_undiffAndROCounts(mat_scaled, diff_genes, title = "Differentiated Cell Types")



make_box_celllinezscore(undiff_mat_scaled, undiff_genes, order = "", origin = "Undifferentiated", title = "Stem-Progenitor Cell Types")
make_box_celllinezscore(RO_mat_scaled, undiff_genes, order = "", origin = "RO", title = "Stem-Progenitor Cell Types")

make_box_celllinezscore(undiff_mat_scaled, diff_genes, order = "", origin = "Undifferentiated", title = "Differentiated Cell Types")
make_box_celllinezscore(RO_mat_scaled, diff_genes, order = "", origin = "RO", title = "Differentiated Cell Types")




make_box_avgzscore_sideByside(undiff_mat_scaled, undiff_genes, order = "", origin = "Undifferentiated", title = "Stem-Progenitor Cell Types")
make_box_avgzscore_sideByside(RO_mat_scaled, undiff_genes, order = "", origin = "RO", title = "Stem-Progenitor Cell Types")

make_box_avgzscore_sideByside(undiff_mat_scaled, diff_genes, order = "", origin = "Undifferentiated", title = "Differentiated Cell Types")
make_box_avgzscore_sideByside(RO_mat_scaled, diff_genes, order = "", origin = "RO", title = "Differentiated Cell Types")









