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

undiff_genes <- read.csv("Retinal_progenitor_genes.csv")
diff_genes <- read.csv("Retinal_differentiated_genes.csv")

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
mat <- mat[rownames(res), , drop = FALSE]
mat_scaled <- t(scale(t(mat)))
mat_scaled <- mat_scaled[complete.cases(mat_scaled), , drop = FALSE]

# boxplot plotting ----------------------------

make_box_celllinezscore <- function(mat_scaled, genes, order = "fwd", title = "") {
  cellline <- sub("_(E8|L3i)", "", colnames(mat_scaled))
  condition <- sub(".*_(E8|L3i).*", "\\1", colnames(mat_scaled))
  
  celltypes <- names(genes)
  boxplot_df_list <- list()
  
  for (i in seq_along(celltypes)) {
    celltype <- celltypes[i]
    gene_set <- genes[[i]]
    sub <- mat_scaled[rownames(mat_scaled) %in% gene_set, , drop = FALSE]
    
    for (g in rownames(sub)) {
      for (j in seq_len(ncol(sub))) {
        boxplot_df_list[[length(boxplot_df_list) + 1]] <- data.frame(
          gene = g,
          zscore = sub[g, j],
          condition = condition[j],
          cellline = cellline[j],
          celltype = celltype
        )
      }
    }
  }
  boxplot_df <- do.call(rbind, boxplot_df_list)
  
  # Calculate mean difference per celltype
  celltype_order <- boxplot_df %>%
    group_by(celltype) %>%
    summarize(diff = mean(zscore[condition == "L3i"]) - mean(zscore[condition == "E8"])) %>%
    arrange(desc(diff)) %>%
    pull(celltype)
  if(order == "rev"){
    celltype_order <- rev(celltype_order)
  }
  
  # Reorder 
  boxplot_df$celltype <- factor(boxplot_df$celltype, levels = celltype_order)
  boxplot_df$condition <- factor(boxplot_df$condition, levels = c("L3i", "E8"))
  
  ggplot(boxplot_df, aes(x = condition, y = zscore, fill = condition)) +
    geom_boxplot(alpha = 0.5) +
    geom_jitter(aes(shape = cellline), width = 0.15, size = 2) +
    facet_wrap(~ celltype, nrow = 1, strip.position = "bottom") +
    scale_fill_manual(values = c("E8" = "red", "L3i" = "blue")) +
    scale_shape_manual(values = c(16, 17, 15)) +
    labs(x = NULL, y = "z-score", shape = "Cell Line",
         fill = "Cell Type", title = title) +
    theme_bw() +
    removeGrid() +
    theme(
      strip.text = element_text(size = 10, face = "bold", angle = 90, hjust = 1),
      strip.placement = "outside",
      strip.background = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank()
    )
}

make_box_avgzscore <- function(mat_scaled, genes, order = "fwd", title = "") {
  cellline <- sub("_(E8|L3i)", "", colnames(mat_scaled))
  condition <- sub(".*_(E8|L3i).*", "\\1", colnames(mat_scaled))
  
  celltypes <- names(genes)
  boxplot_df_list <- list()
  
  for (i in seq_along(celltypes)) {
    celltype <- celltypes[i]
    gene_set <- genes[[i]]
    sub <- mat_scaled[rownames(mat_scaled) %in% gene_set, , drop = FALSE]
    
    L3i_cols <- which(condition == "L3i")
    
    for (g in rownames(sub)) {
      boxplot_df_list[[length(boxplot_df_list) + 1]] <- data.frame(
        gene = g,
        zscore = mean(sub[g, L3i_cols]),
        celltype = celltype
      )
    }
  }
  boxplot_df <- do.call(rbind, boxplot_df_list)
  
  celltype_order <- boxplot_df %>%
    group_by(celltype) %>%
    summarize(mean_z = mean(zscore)) %>%
    arrange(desc(mean_z)) %>%
    pull(celltype)
  if (order == "rev") {
    celltype_order <- rev(celltype_order)
  }
  
  boxplot_df$celltype <- factor(boxplot_df$celltype, levels = celltype_order)
  
  ggplot(boxplot_df, aes(x = celltype, y = zscore, fill = celltype)) +
    geom_boxplot(alpha = 0.5, outlier.shape = NA) +
    geom_jitter(width = 0.15, size = 2) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    labs(x = NULL, y = "Mean Z-Score (L3i)", title = title) +
    theme_bw() +
    removeGrid() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
      legend.position = "none"
    )
}


make_box_rlog <- function(res, genes, order = "fwd", title = "") {
  celltypes <- names(genes)
  boxplot_df_list <- list()
  
  for (i in seq_along(celltypes)) {
    gene_set <- genes[[i]]
    gene_set <- gene_set[gene_set %in% rownames(res)]
    
    boxplot_df_list[[i]] <- data.frame(
      gene = gene_set,
      log2FC = res[gene_set, "log2FoldChange"],
      celltype = celltypes[i]
    )
  }
  
  boxplot_df <- do.call(rbind, boxplot_df_list)
  
  # Order celltypes by mean log2FC
  celltype_order <- boxplot_df %>%
    group_by(celltype) %>%
    summarize(diff = mean(log2FC)) %>%
    arrange(desc(diff)) %>%
    pull(celltype)
  if(order == "rev"){
    celltype_order <- rev(celltype_order)
  }
  
  boxplot_df$celltype <- factor(boxplot_df$celltype, levels = celltype_order)
  
  p <- ggplot(boxplot_df, aes(x = celltype, y = log2FC, fill = celltype)) +
    geom_boxplot(alpha = 0.5) +
    geom_jitter(width = 0.15, size = 2) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    labs( y = "rlog2 (L3i/E8)", x = NULL, title = title) +
    theme_bw() +
    removeGrid() +
    theme(
      strip.text = element_text(size = 12, face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )
  
  return(list(df = boxplot_df, plot = p))
}

make_box_rlog(res, undiff_genes, title = "Undifferentiated Cell Types")
make_box_rlog(res, diff_genes, title = "Differentiated Cell Types")

make_box_avgzscore(mat_scaled, undiff_genes, title = "Undifferentiated Cell Types")
make_box_avgzscore(mat_scaled, diff_genes, title = "Differentiated Cell Types")


make_box_celllinezscore(mat_scaled, undiff_genes, title = "Undifferentiated Cell Types")
make_box_celllinezscore(mat_scaled, diff_genes, title = "Differentiated Cell Types")

log_mat <- log10(mat)
make_box_celllinezscore(log_mat, undiff_genes, title = "Undifferentiated Cell Types")
make_box_celllinezscore(log_mat, diff_genes, title = "Differentiated Cell Types")

