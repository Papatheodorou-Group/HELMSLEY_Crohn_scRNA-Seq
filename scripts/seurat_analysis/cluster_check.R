# extract top markers and plot average expression

TOP <- 30

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5) {
    cat("\nUsage: <obj> <cl_mode> <dea> <tsv> <plot>\n\n")
    q()
}

OBJECT <- args[1]
CL_MODE <- args[2]
DEA <- args[3]
TSV <- args[4]
PLOT <- args[5]

library("Seurat")
library("ComplexHeatmap")

# extract top ranked genes
df <- read.table(DEA, sep = "\t", header = TRUE)
idx_fc <- order(df$avg_log2FC, decreasing = TRUE)
idx_pval <- 1:nrow(df)
M <- as.matrix(data.frame(idx_fc, idx_pval))
avg_idx <- apply(M, 1, mean)
top_genes <- df$geneID[order(avg_idx)]
top_genes <- top_genes[1:TOP]

# extract the average expression
object <- readRDS(OBJECT)
Idents(object) <- object@meta.data[[CL_MODE]]
df <- AverageExpression(object, assays = "RNA", features = top_genes)
tsv <- df$RNA

# output
write.table(tsv, file = TSV, sep = "\t", quote = FALSE)
M <- t(scale(t(tsv)))
pdf(PLOT)
Heatmap(M, cluster_rows = FALSE, cluster_columns = FALSE)
dev.off()

q()

