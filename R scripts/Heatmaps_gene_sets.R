library(ComplexHeatmap)
library(countToFPKM)
library(DESeq2)
library(tidyverse)
library(grid)
library(biomaRt)

# RO_full <- read.table("C:/Users/willllllli/Documents/Dr. Z lab/RNA seq/RO/Counts/count.out",
#                      header = T,
#                      sep = "\t",
#                      comment.char = "#",
#                      stringsAsFactors = F)
# 
# RO_counts <- RO_full[, 7:ncol(RO_full)]
# rownames(RO_counts) <- RO_full$Geneid
# colnames(RO_counts) = c("H9_E8_r1", "H9_E8_r2",  "H9_L3i_r1", "H9_L3i_r2", 
#                             "RUES02_E8", "RUES02_L3i")
# RO_counts <- RO_counts[, c("H9_E8_r1", "H9_E8_r2", "RUES02_E8", "H9_L3i_r1", "H9_L3i_r2", "RUES02_L3i")]
# featureLength <- RO_full$Length
# 
# ff <- fpkm(RO_counts, featureLength, meanFragmentLength = 250)


RO_counts <- read.csv("C:/Users/willllllli/Documents/Dr. Z lab/RNA seq/RO/Counts/RO_counts_gene_symbol.csv")
rownames(RO_counts) <- RO_counts$X
RO_counts$X <- NULL
RO_counts <- RO_counts[, c("H9_E8_r1", "H9_E8_r2", "RUES02_E8", "H9_L3i_r1", "H9_L3i_r2", "RUES02_L3i")]

genes <- read.csv("C:/Users/willllllli/Documents/Dr. Z lab/RNA seq/RO/Retinal_genes.csv")

# Differential expression ----------------------------------
samples <- colnames(RO_counts)
cellline <- sub("_(E8|L3i)?", "", samples)
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
dds <- DESeq(dds)
res <- results(dds, contrast = c("condition", "L3i", "E8"))
res <- res[!is.na(res$padj) & res$padj <0.05, ]
mat <- counts(dds, normalized = TRUE)
# mat <- mat[rownames(res), , drop = FALSE]
mat_scaled <- t(scale(t(mat)))
mat_scaled <- mat_scaled[complete.cases(mat_scaled), , drop = FALSE]



# Heatmap plotting --------------------------------------
celltypes <- names(genes)

all_gene_sets <- list()
slice_labels <- c()

for (i in seq_along(celltypes)) {
  celltype <- celltypes[i]
  gene_set <- genes[[i]]
  sub <- mat_scaled[rownames(mat_scaled) %in% gene_set, ]

  E8_cols  <- grepl("_E8", colnames(sub))
  L3i_cols <- grepl("_L3i", colnames(sub))
  
  heatmap_mat <- cbind(
    sub[, E8_cols],
    sub[, L3i_cols],
    E8_mean  = rowMeans(sub[, E8_cols]),
    L3i_mean = rowMeans(sub[, L3i_cols])
  )
  all_gene_sets[[i]] <- heatmap_mat
  slice_labels <- c(slice_labels, rep(celltype, nrow(heatmap_mat)))
}

combined_mat <- do.call(rbind, all_gene_sets)

ht_combined <- Heatmap(
  combined_mat,
  name = "z-score",
  show_row_names = TRUE, 
  show_column_names = TRUE, 
  cluster_columns = FALSE, 
  cluster_rows = TRUE,
  cluster_row_slices = FALSE,
  row_split = factor(slice_labels, levels = celltypes),  # This creates the slices
  row_title_rot = 0,
  row_gap = unit(2, "mm"),
  row_names_gp = gpar(fontsize = 9),
  column_names_gp = gpar(fontsize = 12),
  row_title_gp = gpar(fontsize = 14, fontface = "bold")
)

draw(ht_combined,
     show_heatmap_legend = T)
