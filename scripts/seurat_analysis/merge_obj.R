# merge the doublet score by doubletFinder and the known doublet cluster to define the final set of doublets
# the final number of doublet will match 10x guidelines, given the number of loaded cells

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
    cat("\nUsage: <dir> <obj1> ... <objn>\n")
    cat("\n<dir>        output objects\n")
    cat("<obj1>       input dirs\n\n")
    q()
}

DIR <- args[1]
N <- length(args)-1

library("Seurat")

obj <- readRDS(args[2])
for (i in 3:length(args)) {
    obj2 <- readRDS(args[i])
    obj <- merge(obj, obj2, merge.data = TRUE)
}

obj <- FindVariableFeatures(obj)
obj <- ScaleData(obj)
obj <- RunPCA(obj)
obj <- RunUMAP(obj, dims=1:20)

Idents(obj) <- obj@meta.data$sample.name

pdf(file.path(DIR, "UMAP_samples.pdf"))
DimPlot(obj, reduction="umap")
dev.off()

saveRDS(obj, file = file.path(DIR, "object.Rds"))


q()

