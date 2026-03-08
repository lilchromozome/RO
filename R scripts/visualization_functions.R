library(ComplexHeatmap)
library(ggplot2)
library(ggsignif)
library(patchwork)
library(ggExtra)

#' Generate Paired Z-Score and Rlog Heatmaps for Gene Sets
#'
#' @param mat_scaled Numeric matrix of z-scored expression values of L3i vs E8.
#' @param rlog Numeric matrix of rlog-transformed expression values
#' @param genes Named list of character vectors
#' @param title Character string for the heatmap column title (default: \code{""}).
#'
#' @return A list with two elements:
#'   \describe{
#'     \item{$ht_combined}mean E8 and mean L3i z-scores per gene.
#'     \item{$ht_combined_rlog} per-sample rlog values E8 and L3i.
#'   }
#'
make_heatmap <- function(mat_scaled, rlog, genes, origin = "", title=""){
  all_gene_sets <- list()
  all_gene_sets_rlog <- list()
  slice_labels <- c()
  celltypes <- names(genes)
  
  for (i in seq_along(celltypes)) {
    celltype <- celltypes[i]
    gene_set <- genes[[i]]
    sub <- mat_scaled[rownames(mat_scaled) %in% gene_set, , drop = FALSE]
    sub_rlog <- rlog[rownames(rlog) %in% gene_set, , drop = FALSE]
    for (g in gene_set){
      # if (!g %in% rownames(mat_scaled)) {print(g)}
    }
    
    E8_cols  <- grepl("_E8", colnames(sub))
    L3i_cols <- grepl("_L3i", colnames(sub))
    
    heatmap_mat <- cbind(
      # sub[, E8_cols],
      # sub[, L3i_cols],
      E8_mean  = rowMeans(sub[, E8_cols]),
      L3i_mean = rowMeans(sub[, L3i_cols])
    )
    
    heatmap_rlog <- sub_rlog
    
    gene_ord <- order(heatmap_mat[, 'L3i_mean'], decreasing = TRUE)
    heatmap_mat <- heatmap_mat[gene_ord, ]
    gene_ord_rlog <- order(heatmap_rlog[, 'L3i'], decreasing = TRUE)
    heatmap_rlog <- heatmap_rlog[gene_ord_rlog, ]
    all_gene_sets[[i]] <- heatmap_mat
    all_gene_sets_rlog[[i]] <- heatmap_rlog
    slice_labels <- c(slice_labels, rep(celltype, nrow(heatmap_mat)))
  }
  
  combined_mat <- do.call(rbind, all_gene_sets)
  combined_rlog <- do.call(rbind, all_gene_sets_rlog)
  
  ht_combined <- Heatmap(
    combined_mat,
    name = 'z-score',
    column_labels = c(paste("E8", origin), paste("L3i", origin)),
    show_row_names = TRUE, 
    show_column_names = TRUE, 
    cluster_columns = FALSE, 
    cluster_rows = FALSE,
    cluster_row_slices = FALSE,
    row_split = factor(slice_labels, levels = celltypes),
    # row_title = NULL,
    row_title_rot = 0,
    row_gap = unit(2, "mm"),
    row_names_gp = gpar(fontsize = 8),
    column_names_gp = gpar(fontsize = 12),
    row_title_gp = gpar(fontsize = 14, fontface = "bold"),
    column_title = title,
    column_title_gp = gpar(fontsize = 16, fontface = "bold"),
  )
  ht_combined_rlog <- Heatmap(
    combined_rlog,
    name = 'rlog2',
    # col = colorRamp2(c(-4, 0, 4), c("blue", "white", "red")),
    column_labels = c(paste("E8", origin), paste("L3i", origin)),
    show_row_names = TRUE, 
    show_column_names = TRUE, 
    cluster_columns = FALSE, 
    cluster_rows = FALSE,
    cluster_row_slices = FALSE,
    row_split = factor(slice_labels, levels = celltypes),
    # row_title = NULL,
    row_title_rot = 0,
    row_gap = unit(2, "mm"),
    row_names_gp = gpar(fontsize = 8),
    column_names_gp = gpar(fontsize = 12),
    row_title_gp = gpar(fontsize = 14, fontface = "bold"),
    column_title = title,
    column_title_gp = gpar(fontsize = 16, fontface = "bold"),
  ) 
  
  draw(ht_combined,
       show_heatmap_legend = T)
  draw(ht_combined_rlog,
       show_heatmap_legend = T)
  return(list(ht_combined = ht_combined, ht_combined_rlog = ht_combined_rlog))
}

#' Generate Paired Heatmaps for Gene Sets
#'
#' @param mat_scaled Numeric matrix of z-scored expression values of L3i vs E8.
#' @param genes Named list of character vectors
#' @param stat_type name of the statistic in matrix
#' @param title Character string for the heatmap column title (default: \code{""}).
#'
#' @return A list with two elements:
#'   \describe{
#'     \item{$ht_combined}mean E8 and mean L3i z-scores per gene.
#'     \item{$ht_combined_rlog} per-sample rlog values E8 and L3i.
#'   }
#'
make_heatmap_undiffAndROCounts <- function(mat_scaled, genes, stat_type = "", title = "") {
  all_gene_sets <- list()
  slice_labels <- c()
  celltypes <- names(genes)
  
  for (i in seq_along(celltypes)) {
    celltype <- celltypes[i]
    gene_set <- genes[[i]]
    sub <- mat_scaled[rownames(mat_scaled) %in% gene_set, , drop = FALSE]
    
    E8_cols     <- grepl("_E8", colnames(sub)) & !grepl("RO", colnames(sub))
    L3i_cols    <- grepl("_L3i", colnames(sub)) & !grepl("RO", colnames(sub))
    E8_RO_cols  <- grepl("_E8", colnames(sub)) & grepl("RO", colnames(sub))
    L3i_RO_cols <- grepl("_L3i", colnames(sub)) & grepl("RO", colnames(sub))
    
    heatmap_mat <- cbind(
      E8_mean     = rowMeans(sub[, E8_cols, drop = FALSE]),
      L3i_mean    = rowMeans(sub[, L3i_cols, drop = FALSE]),
      E8_RO_mean  = rowMeans(sub[, E8_RO_cols, drop = FALSE]),
      L3i_RO_mean = rowMeans(sub[, L3i_RO_cols, drop = FALSE])
    )
    
    gene_ord <- order(heatmap_mat[, "L3i_RO_mean"], decreasing = TRUE)
    heatmap_mat <- heatmap_mat[gene_ord, ]
    all_gene_sets[[i]] <- heatmap_mat
    slice_labels <- c(slice_labels, rep(celltype, nrow(heatmap_mat)))
  }
  
  combined_mat <- do.call(rbind, all_gene_sets)
  
  ht_combined <- Heatmap(
    combined_mat,
    name = stat_type,
    # col = colorRamp2(c(min(combined_mat), max(combined_mat)), c("white", "red")),
    column_labels = c("E8", "L3i", "E8_RO", "L3i_RO"),
    show_row_names = TRUE,
    show_column_names = TRUE,
    cluster_columns = FALSE,
    cluster_rows = TRUE,
    cluster_row_slices = FALSE,
    row_split = factor(slice_labels, levels = celltypes),
    row_title_rot = 0,
    row_gap = unit(2, "mm"),
    row_names_gp = gpar(fontsize = 8),
    column_names_gp = gpar(fontsize = 12),
    row_title_gp = gpar(fontsize = 14, fontface = "bold"),
    column_title = title,
    column_title_gp = gpar(fontsize = 16, fontface = "bold")
  )
  
  return(ht_combined)
}

#' Generate Paired boxplots for Gene Sets
#'
#' @param mat_scaled Numeric matrix of z-scored expression values of L3i vs E8.
#' @param genes Named list of character vectors
#' @param title Character string for the heatmap column title (default: \code{""}).
#'
#' @return A list with two elements:
#'   \describe{
#'     \item{$ht_combined}mean E8 and mean L3i z-scores per gene.
#'     \item{$ht_combined_rlog} per-sample rlog values E8 and L3i.
#'   }
#'
make_boxplot_undiffAndROCounts <- function(mat_scaled, genes, stat_type = "", title = "") {
  celltypes <- names(genes)
  all_gene_sets <- list()
  slice_labels <- c()
  
  for (i in seq_along(celltypes)) {
    celltype <- celltypes[i]
    gene_set <- genes[[i]]
    sub <- mat_scaled[rownames(mat_scaled) %in% gene_set, , drop = FALSE]
    if (nrow(sub) == 0) next
    
    E8_cols     <- grepl("_E8", colnames(sub)) & !grepl("RO", colnames(sub))
    L3i_cols    <- grepl("_L3i", colnames(sub)) & !grepl("RO", colnames(sub))
    E8_RO_cols  <- grepl("_E8", colnames(sub)) & grepl("RO", colnames(sub))
    L3i_RO_cols <- grepl("_L3i", colnames(sub)) & grepl("RO", colnames(sub))
    
    heatmap_mat <- cbind(
      E8       = rowMeans(sub[, E8_cols, drop = FALSE]),
      L3i      = rowMeans(sub[, L3i_cols, drop = FALSE]),
      E8_RO    = rowMeans(sub[, E8_RO_cols, drop = FALSE]),
      L3i_RO   = rowMeans(sub[, L3i_RO_cols, drop = FALSE])
    )
    
    all_gene_sets[[i]] <- heatmap_mat
    slice_labels <- c(slice_labels, rep(celltype, nrow(heatmap_mat)))
  }
  
  combined_mat <- do.call(rbind, all_gene_sets)
  
  # Reshape to long format
  boxplot_df <- data.frame(
    gene = rep(rownames(combined_mat), 4),
    value = c(combined_mat[, "E8"], combined_mat[, "L3i"], combined_mat[, "E8_RO"], combined_mat[, "L3i_RO"]),
    condition = rep(c("E8", "L3i", "E8_RO", "L3i_RO"), each = nrow(combined_mat)),
    celltype = rep(slice_labels, 4)
  )
  
  boxplot_df$celltype <- factor(boxplot_df$celltype, levels = celltypes)
  boxplot_df$condition <- factor(boxplot_df$condition, levels = c("E8", "L3i", "E8_RO", "L3i_RO"))
  
  ggplot(boxplot_df, aes(x = condition, y = value, fill = condition)) +
    geom_boxplot(alpha = 0.5, outlier.shape = NA) +
    geom_jitter(width = 0.15, size = 1.5) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    facet_wrap(~ celltype, nrow = 1, strip.position = "bottom") +
    scale_fill_manual(values = c("E8" = "blue", "L3i" = "red", "E8_RO" = "lightblue", "L3i_RO" = "salmon")) +
    labs(x = NULL, y = stat_type, fill = "Condition", title = title) +
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

#' Boxpslot of Z-Scores by Cell Type and Condition, Grouped by Cell Line
#'
#' @param mat_scaled Numeric matrix of z-scored expression values. Columns must contain `_E8` or `_L3i`.
#' @param genes Named list of character vectors. Names = cell type labels, values = gene symbols.
#' @param order `"fwd"` (default) sorts cell types by descending L3i−E8 difference; `"rev"` reverses.
#' @param title Character string for the plot title (default: `""`).
#'
#' @return A ggplot object. Faceted boxplots (one per cell type) with jittered points shaped by cell line.
make_box_cellline <- function(mat_scaled, genes, stat_type="", order = "", origin = "", title = "") {
  cellline <- sub("_(E8|L3i)", "", colnames(mat_scaled))
  condition <- sub(".*_(E8|L3i).*", "\\1", colnames(mat_scaled))
  
  celltypes <- names(genes)
  boxplot_df_list <- list()
  
  for (i in seq_along(celltypes)) {
    celltype <- celltypes[i]
    gene_set <- genes[[i]]
    if (length(gene_set) == 0) next
    sub <- mat_scaled[rownames(mat_scaled) %in% gene_set, , drop = FALSE]
    missing <- gene_set[!gene_set %in% rownames(mat_scaled)]
    if (length(missing) > 0) {
      message("Missing genes in ", celltype, ": ", paste(missing, collapse = ", "))
    }
    
    for (g in rownames(sub)) {
      for (j in seq_len(ncol(sub))) {
        boxplot_df_list[[length(boxplot_df_list) + 1]] <- data.frame(
          gene = g,
          stat = sub[g, j],
          condition = condition[j],
          cellline = cellline[j],
          celltype = celltype
        )
      }
    }
  }
  boxplot_df <- do.call(rbind, boxplot_df_list)
  
  # # Calculate mean difference per celltype
  celltype_order <- boxplot_df %>%
    group_by(celltype) %>%
    summarize(diff = mean(stat[condition == "L3i"]) - mean(stat[condition == "E8"])) %>%
    arrange(desc(diff)) %>%
    pull(celltype)

  if (order == 'fwd'){
    celltype_order <- celltype_order
  }
  else if (order == "rev") {
    celltype_order <- rev(celltype_order)
  }
  else
    celltype_order <- names(genes)
  
  
  # Reorder 
  boxplot_df$celltype <- factor(boxplot_df$celltype, levels = celltype_order)
  boxplot_df$condition <- factor(boxplot_df$condition, levels = c("E8", "L3i"))
  
  ggplot(boxplot_df, aes(x = condition, y = stat, fill = condition)) +
    geom_boxplot(alpha = 0.5) +
    geom_jitter(aes(shape = cellline), width = 0.15, size = 2) +
    geom_signif(
      comparisons = list(c("E8", "L3i")),
      test = "wilcox.test",
      test.args = list(alternative = "two.sided"),  # E8 < L3i
      map_signif_level = TRUE,                 # shows stars
      textsize = 3
    )+
    facet_wrap(~ celltype, nrow = 1, strip.position = "bottom") +
    scale_fill_manual(values = c("E8" = "blue", "L3i" = "red"),
                      labels = c("E8" = paste("E8", origin), "L3i" = paste("L3i", origin))) +
    scale_shape_manual(values = c(16, 17, 15)) +
    labs(x = NULL, y = stat_type, shape = "Cell Line",
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



#' Boxplot of Mean L3i stats per Gene Across Cell Types
#'
#' @param mat_scaled Numeric matrix of expression values. Columns must contain `_E8` or `_L3i`.
#' @param genes Named list of character vectors. Names = cell type labels, values = gene symbols.
#' @param order `"fwd"` (default) sorts cell types by descending mean; `"rev"` reverses.
#' @param origin the tissue of origin, 'Undifferentiated' or 'RO'
#' @param title Character string for the plot title (default: `""`).
#'
#' @return A ggplot object. One boxplot per cell type showing the mean L3i for each gene.
make_box_avg <- function(mat_scaled, genes, stat_type = "", order = "", origin = "", title = "") {
  cellline <- sub("_(E8|L3i)", "", colnames(mat_scaled))
  condition <- sub(".*_(E8|L3i).*", "\\1", colnames(mat_scaled))
  celltypes <- names(genes)
  boxplot_df_list <- list()
  for (i in seq_along(celltypes)) {
    celltype <- celltypes[i]
    gene_set <- genes[[i]]
    if (length(gene_set) == 0) next
    sub <- mat_scaled[rownames(mat_scaled) %in% gene_set, , drop = FALSE]
    L3i_cols <- which(condition == "L3i")
    for (g in rownames(sub)) {
      boxplot_df_list[[length(boxplot_df_list) + 1]] <- data.frame(
        gene = g,
        stat = mean(sub[g, L3i_cols]),
        celltype = celltype
      )
    }
  }
  boxplot_df <- do.call(rbind, boxplot_df_list)
  celltype_order <- boxplot_df %>%
    group_by(celltype) %>%
    summarize(mean_z = mean(stat)) %>%
    arrange(desc(mean_z)) %>%
    pull(celltype)
  if (order == 'fwd') {
    celltype_order <- celltype_order
  } else if (order == "rev") {
    celltype_order <- rev(celltype_order)
  } else {
    celltype_order <- names(genes)
  }
  boxplot_df$celltype <- factor(boxplot_df$celltype, levels = celltype_order)
  ggplot(boxplot_df, aes(x = celltype, y = stat, fill = celltype)) +
    geom_boxplot(alpha = 0.5, outlier.shape = NA) +
    geom_jitter(width = 0.15, size = 2) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    labs(x = NULL, y = paste0("Mean ", stat_type, " (L3i", origin, ")"), title = title) +
    theme_bw() +
    removeGrid() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
      legend.position = "none"
    )
}

#' Boxplot of Mean L3i AND E8 per Gene Across Cell Types
#'
#' @param mat_scaled Numeric matrix of z-scored expression values. Columns must contain `_E8` or `_L3i`.
#' @param genes Named list of character vectors. Names = cell type labels, values = gene symbols.
#' @param order `"fwd"` (default) sorts cell types by descending mean z-score; `"rev"` reverses.
#' @param origin the tissue of origin, 'Undifferentiated' or 'RO'
#' @param title Character string for the plot title (default: `""`).
#'
#' @return A ggplot object. 2 boxplots per cell type showing the mean L3i AND E8 z-scores for each gene.
make_box_avg_sideByside <- function(mat_scaled, genes, stat_type ="", order = "", origin = "", title = "") {
  cellline <- sub("_(E8|L3i)", "", colnames(mat_scaled))
  condition <- sub(".*_(E8|L3i).*", "\\1", colnames(mat_scaled))
  
  celltypes <- names(genes)
  boxplot_df_list <- list()
  
  for (i in seq_along(celltypes)) {
    celltype <- celltypes[i]
    gene_set <- genes[[i]]
    if (length(gene_set) == 0) next
    sub <- mat_scaled[rownames(mat_scaled) %in% gene_set, , drop = FALSE]
    
    L3i_cols <- which(condition == "L3i")
    E8_cols  <- which(condition == "E8")
    
    for (g in rownames(sub)) {
      boxplot_df_list[[length(boxplot_df_list) + 1]] <- data.frame(
        gene = g,
        stat = mean(sub[g, L3i_cols]),
        condition = "L3i",
        celltype = celltype
      )
      boxplot_df_list[[length(boxplot_df_list) + 1]] <- data.frame(
        gene = g,
        stat = mean(sub[g, E8_cols]),
        condition = "E8",
        celltype = celltype
      )
    }
  }
  boxplot_df <- do.call(rbind, boxplot_df_list)
  
  
  celltype_order <- boxplot_df %>%
    filter(condition == "L3i") %>%
    group_by(celltype) %>%
    summarize(mean_z = mean(stat)) %>%
    arrange(desc(mean_z)) %>%
    pull(celltype)
  if (order == 'fwd'){
    celltype_order <- celltype_order
  }
  else if (order == "rev") {
    celltype_order <- rev(celltype_order)
  }
  else
    celltype_order <- names(genes)
  
  boxplot_df$celltype <- factor(boxplot_df$celltype, levels = celltype_order)
  boxplot_df$condition <- factor(boxplot_df$condition, levels = c("E8", "L3i"))
  
  ggplot(boxplot_df, aes(x = condition, y = stat, fill = condition)) +
    geom_boxplot(alpha = 0.5, outlier.shape = NA) +
    geom_jitter(width = 0.15, size = 2) +
    geom_signif(
      comparisons = list(c("E8", "L3i")),
      test = "wilcox.test",
      test.args = list(alternative = "two.sided"),
      map_signif_level = TRUE,               
      textsize = 3
    )+
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    facet_wrap(~ celltype, nrow = 1, strip.position = "bottom") +
    scale_fill_manual(values = c("E8" = "blue", "L3i" = "red"),
                      labels = c("E8" = paste("E8", origin), "L3i" = paste("L3i", origin))) +
    labs(x = NULL, y = stat_type, fill = "Condition", title = title) +
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


#' Boxplot using normalized count matrix NOT FROM DESEQ2
#' 
#'@param mat normalized matrix
#'@param genes
#'@param order
#'@param origin
#'@param title
#'
#'@return boxplot

make_box_normCounts <- function(mat, genes, order = "", stat_type = "", origin = "", title = "") {
  celltypes <- names(genes)
  boxplot_df_list <- list()
  
  for (i in seq_along(celltypes)) {
    gene_set <- genes[[i]]
    sub <- mat[rownames(mat) %in% gene_set, , drop = FALSE]
    if (nrow(sub) == 0) next
    
    E8_mean  <- rowMeans(sub[, grepl("E8", colnames(sub)), drop = FALSE])
    L3i_mean <- rowMeans(sub[, grepl("L3i", colnames(sub)), drop = FALSE])
    
    log2FC <- log10((L3i_mean + 1) / (E8_mean + 1))
    
    boxplot_df_list[[length(boxplot_df_list) + 1]] <- data.frame(
      gene = names(log2FC),
      log2FC = log2FC,
      condition = "L3i",
      celltype = celltypes[i]
    )
    boxplot_df_list[[length(boxplot_df_list) + 1]] <- data.frame(
      gene = names(log2FC),
      log2FC = -log2FC,
      condition = "E8",
      celltype = celltypes[i]
    )
  }
  boxplot_df <- do.call(rbind, boxplot_df_list)
  
  celltype_order <- boxplot_df %>%
    filter(condition == "L3i") %>%
    group_by(celltype) %>%
    summarize(mean_fc = mean(log2FC)) %>%
    arrange(desc(mean_fc)) %>%
    pull(celltype)
  
  if (order == "fwd") {
    celltype_order <- celltype_order
  } else if (order == "rev") {
    celltype_order <- rev(celltype_order)
  } else {
    celltype_order <- names(genes)
  }
  
  boxplot_df$celltype <- factor(boxplot_df$celltype, levels = celltype_order)
  boxplot_df$condition <- factor(boxplot_df$condition, levels = c("E8", "L3i"))
  
  ggplot(boxplot_df, aes(x = condition, y = log2FC, fill = condition)) +
    geom_boxplot(alpha = 0.5, outlier.shape = NA) +
    geom_jitter(width = 0.15, size = 2) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    facet_wrap(~ celltype, nrow = 1, strip.position = "bottom") +
    scale_fill_manual(values = c("E8" = "blue", "L3i" = "red"),
                      labels = c("E8" = paste("E8", origin), "L3i" = paste("L3i", origin))) +
    labs(x = NULL, y = "log2FC", fill = "Condition", title = title) +
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
