args <- commandArgs(trailingOnly = TRUE)
GENES <- args[1]
FEATURE_MAPPING <- args[2]

# sample GENES lines 
# ENSG00000243485	MIR1302-2HG	Gene Expression
# ENSG00000237613	FAM138A	Gene Expression
# ENSG00000186092	OR4F5	Gene Expression
genes <- read.table(gzfile(GENES), sep = "\t")

# output fields: symbol, ensembl_ID
symbols_unique <- unique(genes[,2])
for (s in symbols_unique) {
    idx <- which(genes[,2] == s)
    if (length(idx) > 1) {
        cat(paste(s, "\n"))
        v <- paste0(s, c("", paste0(".", 1:(length(idx)-1))))
        genes[idx,2] <- v
    }
}
feature_mapping <- genes
colnames(feature_mapping) <- c("ensembl_ID", "symbol", "feature_type")
write.table(feature_mapping, file = FEATURE_MAPPING, quote = FALSE, sep = "\t", row.names = FALSE)

q()

