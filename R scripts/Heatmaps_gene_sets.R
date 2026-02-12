library(ComplexHeatmap)
library(DESeq2)
library(circlize)
library(tidyverse)
library(org.Hs.eg.db)
library(AnnotationDbi)

# RO_full <- read.table("C:/Users/willllllli/Documents/Dr. Z lab/RNA seq/RO/Counts/count.out", # Windows
RO_full <- read.table("~/Documents/Zambidis lab/RNAseq/RO/Counts/count.out",
                     header = T,
                     sep = "\t",
                     comment.char = "#",
                     stringsAsFactors = F)

RO_counts <- RO_full[, 7:ncol(RO_full)]
rownames(RO_counts) <- RO_full$Geneid
colnames(RO_counts) = c("H9_E8_r1", "H9_E8_r2",  "H9_L3i_r1", "H9_L3i_r2",
                            "RUES02_E8", "RUES02_L3i")
RO_counts <- RO_counts[, c("H9_E8_r1", "H9_E8_r2", "RUES02_E8", "H9_L3i_r1", "H9_L3i_r2", "RUES02_L3i")]
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

genes <- read.csv("~/Documents/Zambidis lab/RNAseq/RO/Retinal_genes.csv") # MAC
genes <- genes[, c("Retinal.Stem.Cell", "Ganglion", "HA.cell", "Photoreceptor", "Bipolar","Muller") ]

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

# TPM calculation -------------------

rpk <- RO_counts / (featureLength / 1000)  # Reads per kilobase
tpm <- t(t(rpk) / colSums(rpk) * 1e6)  # TPM normalization

# RO_norm <- RO_counts/colSums(RO_counts) * 1e6
# rlog_counts <- log2(RO_norm/(rowMeans(RO_norm) + 1 ))

rld <- rlog(dds, blind = FALSE)
rlog_counts <- assay(rld)


# Heatmap plotting --------------------------------------
celltypes <- names(genes)

all_gene_sets <- list()
slice_labels <- c()

for (i in seq_along(celltypes)) {
  celltype <- celltypes[i]
  gene_set <- genes[[i]]
  sub <- mat_scaled[rownames(mat_scaled) %in% gene_set, ]
  # sub <- tpm[rownames(tpm) %in% gene_set, ]

  E8_cols  <- grepl("_E8", colnames(sub))
  L3i_cols <- grepl("_L3i", colnames(sub))
  
  heatmap_mat <- cbind(
    sub[, E8_cols],
    sub[, L3i_cols],
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
  name = "z-score",
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
