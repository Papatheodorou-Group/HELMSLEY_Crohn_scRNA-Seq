
N_TOP_MARKERS <- 20

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
    cat("error\n")
    q()
}

DIR1 <- args[1]
DIR2 <- args[2]
OUT <- args[3]

files1 <- list.files(DIR1)
files2 <- list.files(DIR2)
files1 <- files1[grep("DEG.+-all.tsv", files1)]
files2 <- files2[grep("DEG.+-all.tsv", files2)]

m <- length(files1)
n <- length(files2)
M <- matrix(0, nrow = m, ncol = n)
names1 <- rownames(M) <- gsub("DEG_MAST_cl", "", gsub("-all.tsv", "", files1))
names2 <- colnames(M) <- gsub("DEG_MAST_cl", "", gsub("-all.tsv", "", files2))
M <- M[order(as.numeric(rownames(M))), order(as.numeric(colnames(M)))]

l1 <- list()
i <- 1
for (f1 in files1) {
    df1 <- read.table(file.path(DIR1,f1), sep = "\t", header = TRUE)
    markers1 <- df1$geneID[order(df1$avg_log2FC, decreasing = TRUE)]
    l1[[i]] <- as.character(markers1[1:N_TOP_MARKERS])
    i <- i + 1
}
names(l1) <- names1

l2 <- list()
j <- 1
for (f2 in files2) {
    df2 <- read.table(file.path(DIR2,f2), sep = "\t", header = TRUE)
    markers2 <- df2$geneID[order(df2$avg_log2FC, decreasing = TRUE)]
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
