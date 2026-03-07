library(ComplexHeatmap)
library(DESeq2)
library(circlize)
library(tidyverse)
library(patchwork)
library(ggExtra)
library(limma)

undiff <- read.delim("Counts/RO_counts_undifferentiated_gene_symbol.csv", sep=",")
rownames(undiff) <- undiff$X
undiff$X <- NULL
undiff <- undiff[, c("E5C3_E8", "H9_E8", "RUES02_E8", "E5C3_L3i","H9_L3i",  "RUES02_L3i")]

RO <- read.delim("Counts/RO_counts_gene_symbol.csv", sep=",")
rownames(RO) <- RO$X
RO$X <- NULL
colnames(RO) <- c("H9_E8_r1_RO", "H9_E8_r2_RO", "H9_L3i_r1_RO",  "H9_L3i_r2_RO", "RUES02_E8_RO", "RUES02_L3i_RO")
RO <- RO[, c("H9_E8_r1_RO", "H9_E8_r2_RO", "RUES02_E8_RO", "H9_L3i_r1_RO",  "H9_L3i_r2_RO",  "RUES02_L3i_RO")]

all <- cbind(undiff,RO)

dorgau_progenitor_genes <- read.csv("retinal_progenitor_genes.csv")
Liu_progenitor_genes <- read.csv("retinal_progenitor_genes.csv")
diff_genes <- read.csv("diff_dorgau.csv")
malika_genes <- data.frame(Malika_genes = c('NANOG', 'POU5F1', 'SOX2', 'MECOM', 'RELN', 'COL9A1', 'PCDH7',
                 'ATOH7', 'PAX6', 'VSX2', 'SOX1', 'RAX', 'LHX2', 'SIX3', 'SIX6',
                 'ASCL1', 'HES6', 'ATOH7', 'DLL3', 'NRL', 'CRX', 'OTX2', 'OLIG2',
                 'MITF', 'GJA1'))

#===========================================
## Undiff Differential expression ----------------------------------
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
dds_undiff_deseq <- DESeq(dds_undiff)
undiff_res <- results(dds_undiff_deseq, contrast = c("condition", "L3i", "E8"))

# rlog2
undiff_log2FC <- cbind(
  E8 = -undiff_res$log2FoldChange,
  L3i = undiff_res$log2FoldChange
)
rownames(undiff_log2FC) <- rownames(undiff_res)

rlog_counts_undiff <- rlog(dds_undiff)
rlog_mat <- assay(rlog_counts_undiff)
rlog_mat <- rlog_mat - rowMeans(rlog_mat)
rlog_var <- apply(rlog_mat, 1, var)
rlog_mat <- rlog_mat[names(rlog_var > 0), ]
rlog_mat_undiff <- t(scale(t(rlog_mat)))
head(rlog_mat)
undiff_rlog <- data.frame(
  E8  = rowMeans(rlog_mat_undiff[, grepl("E8", colnames(rlog_mat_undiff))]),
  L3i = rowMeans(rlog_mat_undiff[, grepl("L3i", colnames(rlog_mat_undiff))])
)

# Z-score 
undiff_mat <- counts(dds_undiff_deseq, normalized = TRUE)
undiff_mat_scaled <- t(scale(t(undiff_mat)))
undiff_mat_scaled <- undiff_mat_scaled[complete.cases(undiff_mat_scaled), , drop = FALSE]


#===========================================
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

#===========================================
## Combined differential expression-----------------------------------
samples <- colnames(all)
is_RO <- grepl("RO", samples)
media <- sub(".*_(E8|L3i).*", "\\1", samples)
group <- paste0(media, ifelse(is_RO, "_RO", ""))

coldata <- data.frame(
  group = factor(group, levels = c("E8", "L3i", "E8_RO", "L3i_RO")),
  row.names = samples
)

all(colnames(all) == rownames(coldata))

dds_all <- DESeqDataSetFromMatrix(
  countData = all,
  colData = coldata,
  design = ~ group
)
dds_all <- dds_all[rowSums(counts(dds_all)) > 10, ]
dds_all_deseq <- DESeq(dds_all)

# rlog2
rlog_counts_all <- rlog(dds_all_deseq)
rlog_mat_all <- assay(rlog_counts_all)
rlog_mat_all <- removeBatchEffect(rlog_mat_all, batch = coldata$cellline)
rlog_mat_all <- rlog_mat_all - rowMeans(rlog_mat_all)
rlog_var_all <- apply(rlog_mat_all, 1, var)
rlog_mat_all <- rlog_mat_all[names(rlog_var_all > 0), ]
rlog_mat_all <- t(scale(t(rlog_mat_all)))

# Z-score 
all_mat <- counts(dds_all_deseq, normalized = TRUE)
all_mat_scaled <- t(scale(t(all_mat)))
all_mat_scaled <- all_mat_scaled[complete.cases(all_mat_scaled), , drop = FALSE]

## Visualizations ----------------------------------
source('R scripts/visualization_functions.R')
R2 <- undiff[,grepl("RUES02", colnames(undiff))]
cpm <- R2 / colSums(R2) * 1e6
make_box_normCounts(cpm, undiff_genes, origin = "undifferentiated")



make_heatmap(undiff_mat_scaled, undiff_rlog, undiff_genes, 
             origin = "Undifferentiated", title="Stem-Progenitor Cell Types")
make_heatmap(undiff_mat_scaled, undiff_rlog, diff_genes, 
             origin = "Undifferentiated", title="Differentiated Cell Types")
make_heatmap(RO_mat_scaled, RO_rlog, undiff_genes, origin = "RO", title="Stem-Progenitor Cell Types")
make_heatmap(RO_mat_scaled, RO_rlog, diff_genes, origin = "RO", title="Differentiated Cell Types")

make_heatmap_undiffAndROCounts(all_mat_scaled, undiff_genes, stat_type="z-score", title = "Stem-Progenitor Cell Types")
make_heatmap_undiffAndROCounts(all_mat_scaled, diff_genes, stat_type="z-score", title = "Differentiated Cell Types")
make_heatmap_undiffAndROCounts(rlog_mat_all, undiff_genes, stat_type="rlog2", title = "Stem-Progenitor Cell Types")
make_heatmap_undiffAndROCounts(rlog_mat_all, diff_genes, stat_type="rlog2", title = "Differentiated Cell Types")


make_box_cellline(rlog_mat_undiff, undiff_genes, order = "", stat_type = 'rlog2',
                  origin = "Undifferentiated", title = "Stem-Progenitor Cell Types")
make_box_cellline(rlog_mat_RO, undiff_genes, order = "", stat_type = 'rlog2',
                  origin = "RO", title = "Stem-Progenitor Cell Types")

make_box_cellline(rlog_mat_undiff, diff_genes, order = "rev", stat_type = 'rlog2',
                  origin = "Undifferentiated", title = "Differentiated Cell Types")
make_box_cellline(rlog_mat_RO, diff_genes, order = "rev", stat_type = 'rlog2',
                  origin = "RO", title = "Differentiated Cell Types")


make_box_cellline(rlog_mat_undiff, undiff_genes, order = "", origin = "Undifferentiated", title = "Stem-Progenitor Cell Types")
make_box_cellline(rlog_mat_RO, undiff_genes, order = "", origin = "RO", title = "Stem-Progenitor Cell Types")

make_box_cellline(rlog_mat_undiff, diff_genes, order = "rev", origin = "Undifferentiated", title = "Differentiated Cell Types")
make_box_cellline(rlog_mat_RO, diff_genes, order = "rev", origin = "RO", title = "Differentiated Cell Types")




make_box_avg_sideByside(rlog_mat_undiff, undiff_genes, order = "", origin = "Undifferentiated", title = "Stem-Progenitor Cell Types")
make_box_avg_sideByside(RO_mat_scaled, undiff_genes, order = "", origin = "RO", title = "Stem-Progenitor Cell Types")

make_box_avg_sideByside(rlog_mat_undiff, diff_genes, order = "rev", origin = "Undifferentiated", title = "Differentiated Cell Types")
make_box_avg_sideByside(RO_mat_scaled, diff_genes, order = "rev", origin = "RO", title = "Differentiated Cell Types")


make_boxplot_undiffAndROCounts(rlog_mat_all, undiff_genes, stat_type = "rlog2", title = "Stem-Progenitor Cell Types")
make_boxplot_undiffAndROCounts(rlog_mat_all, diff_genes, stat_type = "rlog2", title = "Differentiated Cell Types")


# Malika genes -----------------
make_heatmap(RO_mat_scaled, RO_rlog, malika_genes, origin = "RO", title="Malika genes")
make_heatmap(undiff_mat_scaled, undiff_rlog, malika_genes, origin = "undifferentiated", title="Malika genes")

R2 <- undiff
# R2 <- undiff[,grepl("RUES02", colnames(undiff))]
R2 <- R2 / colSums(R2) * 1e6
# R2 <- R2/unlist(R2["ACTB", ])

sub <- R2[rownames(R2) %in% malika_genes$Malika_genes, , drop = FALSE]
sub <- data.frame(
  E8 = rowMeans(sub[, grepl("E8", colnames(sub))]),
  L3i = rowMeans(sub[, grepl("L3i", colnames(sub))])
)

df <- data.frame(
  gene = rownames(sub),
  log10FC = -log10((sub[, 2]+1)/(sub[,1]+1))  # adjust column index if needed
)
df$gene <- factor(df$gene, levels = df$gene[order(df$log10FC)])

df <- data.frame(
  row.names = rownames(sub),
  E8 = -log2(sub[, 2]/(sub[,1]+1)),
  L3i = log2(sub[, 2]/(sub[,1]+1))
)

ggplot(df, aes(x = gene, y = log10FC, fill = log10FC > 0)) +
  geom_bar(stat = "identity") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  scale_fill_manual(values = c("TRUE" = "blue", "FALSE" = "red"), guide = "none") +
  labs(x = NULL, y = "-log10 Fold Change (L3i/E8)", title = "Mean undifferentiated normalized to total counts") +
  theme_bw() + 
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    legend.position = "none"
  )

make_heatmap(R2, df, malika_genes, origin = "undifferentiated", title="Malika genes")



