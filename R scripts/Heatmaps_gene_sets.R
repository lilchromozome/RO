library(ComplexHeatmap)
library(DESeq2)
library(circlize)
library(tidyverse)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(patchwork)

# RO_full <- read.table("C:/Users/willllllli/Documents/Dr. Z lab/RNA seq/RO/Counts/count.out", # Windows
RO_full <- read.table("~/Documents/Zambidis lab/RNAseq/RO/Counts/count.out",
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

undiff_genes <- read.csv("~/Documents/Zambidis lab/RNAseq/RO/Retinal_progenitor_genes.csv") # MAC
diff_genes <- read.csv("~/Documents/Zambidis lab/RNAseq/RO/Retinal_differentiated_genes.csv") # MAC

undiff_cell_types <- unique(gsub("\\.(Liu|Dorgau)$", "", names(diff_genes)))

diff_genes_list <- lapply(undiff_cell_types, function(ct) {
  matching_cols <- grep(paste0("^", ct, "\\.(Liu|Dorgau)$"), names(diff_genes), value = TRUE)
  unique(na.omit(unlist(diff_genes[matching_cols])))
})

names(diff_genes_list) <- undiff_cell_types
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
dds <- DESeq(dds)
res <- results(dds, contrast = c("condition", "L3i", "E8"))
res <- res[!is.na(res$padj) & res$padj <0.05, ]
mat <- counts(dds, normalized = TRUE)
# mat <- mat[rownames(res), , drop = FALSE]
mat_scaled <- t(scale(t(mat)))
mat_scaled <- mat_scaled[complete.cases(mat_scaled), , drop = FALSE]

# Heatmap plotting --------------------------------------
make_heatmap <- function(mat_scaled, genes, title){
  all_gene_sets <- list()
  slice_labels <- c()
  celltypes <- names(genes)
  
  for (i in seq_along(celltypes)) {
    celltype <- celltypes[i]
    gene_set <- genes[[i]]
    sub <- mat_scaled[rownames(mat_scaled) %in% gene_set, ]
    for (g in gene_set){
      if (!g %in% gene_set) {print(g)}
    }
    
    # sub <- tpm[rownames(tpm) %in% gene_set, ]
    
    E8_cols  <- grepl("_E8", colnames(sub))
    L3i_cols <- grepl("_L3i", colnames(sub))
    
    heatmap_mat <- cbind(
      # sub[, E8_cols],
      # sub[, L3i_cols],
      E8_mean  = rowMeans(sub[, E8_cols]),
      L3i_mean = rowMeans(sub[, L3i_cols])
    )
    gene_ord <- order(heatmap_mat[, 'L3i_mean'], decreasing = TRUE)
    heatmap_mat <- heatmap_mat[gene_ord, ]
    all_gene_sets[[i]] <- heatmap_mat
    slice_labels <- c(slice_labels, rep(celltype, nrow(heatmap_mat)))
  }
  
  combined_mat <- do.call(rbind, all_gene_sets)
  
  ht_combined <- Heatmap(
    combined_mat,
    name = 'z-score',
    show_row_names = TRUE, 
    show_column_names = TRUE, 
    cluster_columns = FALSE, 
    cluster_rows = FALSE,
    cluster_row_slices = FALSE,
    row_split = factor(slice_labels, levels = celltypes),
    row_title_rot = 0,
    row_gap = unit(2, "mm"),
    row_names_gp = gpar(fontsize = 8),
    column_names_gp = gpar(fontsize = 12),
    row_title_gp = gpar(fontsize = 14, fontface = "bold")
  )
  
  draw(ht_combined,
       show_heatmap_legend = T)
  return(ht_combined)
}

undiff_ht <- make_heatmap(mat_scaled, undiff_genes)
diff_ht <- make_heatmap(mat_scaled, diff_genes_df)

