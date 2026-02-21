library(DESeq2)
library(org.Hs.eg.db)
library(AnnotationDbi)

# featureCounts -a Homo_sapiens.GRCh38.109.gtf -o count.out -T 8 *.bam 

counts <- read.table("/Users/willli/Documents/Zambidis lab/RNAseq/RO/Counts/count_undifferentiated.out",
                     header = T,
                     sep = "\t",
                     comment.char = "#",
                     stringsAsFactors = F)

counts_matrix <- counts[, 7:ncol(counts)]
rownames(counts_matrix) <- counts$Geneid
colnames(counts_matrix) = c("E5C3_E8", "E5C3_L3i",  "H9_E8", "H9_L3i", 
                            "RUES02_E8", "RUES02_L3i")

write.csv(counts_matrix, "/Users/willli/Documents/Zambidis lab/RNAseq/RO/Counts/RO_undifferentiated_counts.csv")

ens <- rownames(counts_matrix)
symbols <- mapIds(org.Hs.eg.db,
                  keys = ens,
                  keytype = "ENSEMBL",
                  column = "SYMBOL",
                  multiVals = "first")
counts_matrix$gene_symbol <- symbols
rownames(counts_matrix) <- make.unique(ifelse(is.na(symbols), ens, symbols))
counts_matrix$gene_symbol <- NULL

write.csv(counts_matrix, "/Users/willli/Documents/Zambidis lab/RNAseq/RO/Counts/RO_counts_undifferentiated_gene_symbol.csv")
