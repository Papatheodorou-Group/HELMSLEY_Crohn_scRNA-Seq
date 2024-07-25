args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
    cat("Usage: <object> <matrix_dir>\n")
    q()
}
OBJECT <- args[1]
MATRIX_DIR <- args[2]

library("Seurat")
library("Matrix")

object <- readRDS(OBJECT)

writeMM(t(object@assays$RNA@counts), file = file.path(MATRIX_DIR, "matrix.mtx")) # transpose for scanpy!!!
write(rownames(object@assays$RNA@counts), file = file.path(MATRIX_DIR, "genes.tsv"))
write(colnames(object@assays$RNA@counts), file = file.path(MATRIX_DIR, "barcodes.tsv"))

q()

