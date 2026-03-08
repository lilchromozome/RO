library(ComplexHeatmap)
library(DESeq2)
library(circlize)
library(tidyverse)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(patchwork)
library(ggExtra)
setwd("C:/Users/willllllli/Documents/Dr. Z lab/RNA seq/RO") # Windows
# setwd("~/Documents/Zambidis lab/RNAseq/RO") # MAC

RO <- read.delim("Counts/RO_counts_gene_symbol.csv", sep=",")
rownames(RO) <- RO$X
RO$X <- NULL
colnames(RO) <- c("H9_E8_r1_RO", "H9_E8_r2_RO", "H9_L3i_r1_RO",  "H9_L3i_r2_RO", "RUES02_E8_RO", "RUES02_L3i_RO")
RO <- RO[, c("H9_E8_r1_RO", "H9_E8_r2_RO", "RUES02_E8_RO", "H9_L3i_r1_RO",  "H9_L3i_r2_RO",  "RUES02_L3i_RO")]

dorgau_progenitor_genes <- read.csv("dorgau_progenitor.csv")
liu_progenitor_genes <- read.csv("Liu_progenitor.csv")
diff_genes <- read.csv("dorgau_diff.csv")


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
dds_RO_deseq <- DESeq(dds_RO)
RO_res <- results(dds_RO_deseq, contrast = c("condition", "L3i", "E8"))

# rlog2
RO_log2FC <- cbind(
  E8 = -RO_res$log2FoldChange,
  L3i = RO_res$log2FoldChange
)
rownames(RO_log2FC) <- rownames(RO_res)

rlog_counts_RO <- rlog(dds_RO_deseq)
rlog_mat <- assay(rlog_counts_RO)
rlog_mat <- rlog_mat - rowMeans(rlog_mat)
rlog_var <- apply(rlog_mat, 1, var)
rlog_mat <- rlog_mat[names(rlog_var > 0), ]
rlog_mat_RO <- t(scale(t(rlog_mat)))
RO_rlog <- data.frame(
  E8  = rowMeans(rlog_mat_RO[, grepl("E8", colnames(rlog_mat_RO))]),
  L3i = rowMeans(rlog_mat_RO[, grepl("L3i", colnames(rlog_mat_RO))])
)
# Z-score 
RO_mat <- counts(dds_RO_deseq, normalized = TRUE)
RO_mat_scaled <- t(scale(t(RO_mat)))
RO_mat_scaled <- RO_mat_scaled[complete.cases(RO_mat_scaled), , drop = FALSE]

# Plotting --------------------------------------
source("R scripts/visualization_functions.R")
make_heatmap(RO_mat_scaled, RO_rlog, dorgau_progenitor_genes, origin = "RO", title="Stem-Progenitor Cell Types Dorgau")
make_heatmap(RO_mat_scaled, RO_rlog, liu_progenitor_genes, origin = "RO", title="Stem-Progenitor Cell Types Liu")
make_heatmap(RO_mat_scaled, RO_rlog, diff_genes, origin = "RO", title="Differentiated Cell Types")

make_box_avg_sideByside(RO_rlog, dorgau_progenitor_genes, order = "", 
                        stat_type = 'rlog2', origin = "RO", title = "Stem-Progenitor Cell Types Dorgau")
make_box_avg_sideByside(RO_rlog, liu_progenitor_genes, order = "", 
                        stat_type = 'rlog2', origin = "RO", title = "Stem-Progenitor Cell Types Liu")
make_box_avg_sideByside(RO_rlog, diff_genes, order = "", 
                        stat_type = 'rlog2', origin = "RO", title = "Differentiated Cell Types")
