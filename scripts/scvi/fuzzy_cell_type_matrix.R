MIN_CELLS <- 10

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
    cat("Usage: <soft_q_lab> <GCA_stats> <out_matrix>\n")
    q()
}
SOFT_Q_LAB <- args[1]
GCA_STATS <- args[2]
OUT_PREFIX <- args[3]

# load
soft_q_lab <- read.table(SOFT_Q_LAB, sep = "\t")
rownames(soft_q_lab) <- soft_q_lab[,1]
colnames(soft_q_lab) <- soft_q_lab[1,]
soft_q_lab <- soft_q_lab[2:nrow(soft_q_lab),2:ncol(soft_q_lab)]
soft_q_lab <- as.matrix(soft_q_lab)
soft_q_lab <- apply(soft_q_lab, 2, as.numeric)

# filter out rare labels
soft_q_lab <- soft_q_lab[,colSums(soft_q_lab) >= MIN_CELLS]

# compute fuzzy cell-type assignment
M <- matrix(NA, ncol = ncol(soft_q_lab), nrow = ncol(soft_q_lab))
rownames(M) <- colnames(M) <- colnames(soft_q_lab)
for (i in 1:nrow(M)) {
    for (j in 1:ncol(M)) {
        M[i,j] <- sum(apply(soft_q_lab[,c(i,j)], 1, min))
    }
}

# normalise
v <- colSums(soft_q_lab)
M_norm <- M
for (i in 1:nrow(M)) {
    for (j in 1:ncol(M)) {
        M_norm[i,j] <- sqrt((M[i,j]*M[i,j])/(v[i]*v[j]))
    }
}

# diag(M) <- diag(M_norm) <- rep(NA, nrow(M))

# assign category
df_cat_celltype <- read.table(GCA_STATS, sep = "\t", header = TRUE)
orig_celltype <- df_cat_celltype$Integrated_05
idx <- match(colnames(M), df_cat_celltype$Integrated_05)
category <- df_cat_celltype$category[idx]

# print
out <- paste0(OUT_PREFIX, ".tsv")
out_norm <- paste0(OUT_PREFIX, "_norm.tsv")
out_label <- paste0(OUT_PREFIX, "_label.txt")
write.table(M, file = out, quote = FALSE, sep = "\t")
write.table(M_norm, file = out_norm, quote = FALSE, sep = "\t")
write(category, file = out_label)

q()

