args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
    cat("Usage: <tsv> <out_prefix>\n")
    q()
}
TSV <- args[1]
OUT_PREFIX <- args[2]

library("ggplot2")

df <- read.table(TSV, sep = "\t", header = TRUE)
df <- cbind(df, count = rep(1,nrow(df)))

g <- ggplot(data = df, aes(x = Sample.name, y = count, fill = category)) + theme_minimal() + geom_bar(stat="identity", position = "fill")
pdf(paste0(OUT_PREFIX, ".pdf"))
print(g)
dev.off()

categories <- unique(df$category)
for (c in categories) {
    g <- ggplot(data = df[df$category == c,], aes(x = Sample.name, y = count, fill = Integrated_05)) + theme_minimal() + geom_bar(stat="identity", position = "fill")
    pdf(paste0(OUT_PREFIX, "_", c, ".pdf"))
    print(g)
    dev.off()
}

q()
