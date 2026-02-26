library(DESeq2)
library(org.Hs.eg.db)
library(AnnotationDbi)

# featureCounts -a Homo_sapiens.GRCh38.109.gtf -o count.out -T 8 *.bam 

counts <- read.table("/Users/willli/Documents/Zambidis lab/RNAseq/RO/Counts/count.out",
                     header = T,
                     sep = "\t",
                     comment.char = "#",
                     stringsAsFactors = F)

counts_matrix <- counts[, 7:ncol(counts)]
rownames(counts_matrix) <- counts$Geneid
colnames(counts_matrix) = c("H9_E8_r1", "H9_E8_r2",  "H9_L3i_r1", "H9_L3i_r2", 
                            "RUES02_E8", "RUES02_L3i")

write.csv(counts_matrix, "/Users/willli/Documents/Zambidis lab/RNAseq/RO/Counts/RO_counts.csv")

ens <- rownames(counts_matrix)
symbols <- mapIds(org.Hs.eg.db,
                  keys = ens,
                  keytype = "ENSEMBL",
                  column = "SYMBOL",
                  multiVals = "first")
counts_matrix$gene_symbol <- symbols
rownames(counts_matrix) <- make.unique(ifelse(is.na(symbols), ens, symbols))
counts_matrix$gene_symbol <- NULL

write.csv(counts_matrix, "/Users/willli/Documents/Zambidis lab/RNAseq/RO/Counts/RO_counts_gene_symbol.csv")
