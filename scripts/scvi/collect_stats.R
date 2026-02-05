
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
    cat("\nUsage: <r_obj_prefix> <out> <meta>\n")
    q()
}

R_OBJ_PREFIX <- args[1]
OUT <- args[2]
META <- args[3:length(args)]

DF1 <- paste0(R_OBJ_PREFIX, "_", META[1], ".tsv")
DF2 <- paste0(R_OBJ_PREFIX, "_", META[2], ".tsv")

df_cat <- read.table(DF1, sep = "\t", header = TRUE)
df_ct <- read.table(DF2, sep = "\t", header = TRUE)

df <- data.frame(df_cat$labels_scanvi, df_ct$labels_scanvi)
df_out <- table(df)

tab <- as.data.frame(df_out)
colnames(tab) <- c(META, "count")
tab <- tab[tab$count > 0,]
tab <- tab[order(tab[,1]),]

write.table(tab, file = OUT, sep = "\t", row.names = FALSE, quote = FALSE)

q()


