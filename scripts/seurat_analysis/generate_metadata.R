args <- commandArgs(trailingOnly = TRUE)
OBJECT <- args[1]
OUT <- args[2]

library("Seurat")
library("Matrix")

object <- readRDS(OBJECT)
write.table(object@meta.data, file = OUT, quote = FALSE, sep = "\t")

q()

