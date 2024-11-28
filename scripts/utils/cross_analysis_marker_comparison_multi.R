
N_TOP_MARKERS <- 20

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
    cat("error\n")
    q()
}

FILE1 <- args[1]
FILE2 <- args[2]
OUT <- args[3]

df1 <- read.table(FILE1, sep = "\t", header = TRUE)
df2 <- read.table(FILE2, sep = "\t", header = TRUE)

clusters1 <- unique(df1$cluster)
clusters2 <- unique(df2$cluster)
m <- length(clusters1)
n <- length(clusters2)
M <- matrix(0, nrow = m, ncol = n)
names1 <- rownames(M) <- clusters1
names2 <- colnames(M) <- clusters2

l1 <- list()
i <- 1
for (cl1 in clusters1) {
    df_sub <- df1[df1$cluster == cl1,]
    markers1 <- df_sub$gene[order(df_sub$avg_log2FC, decreasing = TRUE)]
    l1[[i]] <- as.character(markers1[1:N_TOP_MARKERS])
    i <- i + 1
}
names(l1) <- names1

l2 <- list()
j <- 1
for (cl2 in clusters2) {
    df_sub <- df2[df2$cluster == cl2,]
    markers2 <- df_sub$gene[order(df_sub$avg_log2FC, decreasing = TRUE)]
    l2[[j]] <- as.character(markers2[1:N_TOP_MARKERS])
    j <- j + 1
}
names(l2) <- names2

for (i in 1:length(l1)) {
    idx1 <- match(names(l1)[i], rownames(M))
    for (j in 1:length(l2)) {
        idx2 <- match(names(l2)[j], colnames(M))
        M[idx1,idx2] <- sum(l1[[i]] %in% l2[[j]])
    }
}

write.table(M, file = OUT, sep = "\t", quote = FALSE)

q()
