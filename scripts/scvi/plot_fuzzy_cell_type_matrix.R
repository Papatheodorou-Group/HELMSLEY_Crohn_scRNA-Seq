args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
    cat("Usage: <out_matrix>\n")
    q()
}
OUT_PREFIX <- args[1]

# load
out <- paste0(OUT_PREFIX, ".tsv")
out_norm <- paste0(OUT_PREFIX, "_norm.tsv")
out_label <- paste0(OUT_PREFIX, "_label.txt")

M <- read.table(out, sep = "\t", header = TRUE)
M_norm <- read.table(out_norm, sep = "\t", header = TRUE)
label <- read.table(out_label, sep = "\t")[,1]

colnames(M) <- colnames(M_norm) <- rownames(M)
diag(M) <- diag(M_norm) <- rep(NA, nrow(M)) 

library("ComplexHeatmap")
library("circlize")

col <- c("lemonchiffon","blue")
col_fun <- colorRamp2(c(0,0.1), col)

H <- Heatmap(M, col = col, row_split = label, column_split = label)
H_norm <- Heatmap(M_norm, col = col_fun, row_split = label, column_split = label)

plot <- paste0(OUT_PREFIX, ".pdf")
plot_norm <- paste0(OUT_PREFIX, "_norm.pdf")

pdf(plot, width = 10, height = 10)
draw(H)
dev.off()

pdf(plot_norm, width = 10, height = 10)
draw(H_norm)
dev.off()

q()

