# select a cluster, subset the object, and do a subclustering

DIMS <- 10

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5) {
    cat("\nUsage: <obj> <cl_mode> <cl_id> <res> <cat> <plot_prefix> <dea>\n\n")
    q()
}

OBJECT <- args[1]
CL_MODE <- args[2]
CL_ID <- args[3]
RES <- as.numeric(args[4])
CAT <- args[5]
PLOT_PREFIX <- args[6]
DEA <- args[7]

CL_LIST <- unlist(strsplit(CL_ID, split = ","))

library("Seurat")

object <- readRDS(OBJECT)
Idents(object) <- object@meta.data[[CL_MODE]] # this is just to ensure that subset() works (random bug in Seurat...)
object_cl <- subset(object, cells = colnames(object)[object@meta.data[[CL_MODE]] %in% CL_LIST])
object_cl <- FindVariableFeatures(object_cl)
object_cl <- ScaleData(object_cl)
object_cl <- RunPCA(object_cl)
object_cl <- FindNeighbors(object_cl, reduction="pca", dims = 1:10)
object_cl <- FindClusters(object_cl, res = RES)
object_cl <- RunUMAP(object_cl, reduction="pca", dims = 1:DIMS)

write.table(Idents(object_cl), paste0(PLOT_PREFIX, "_clusters.tsv"), quote = FALSE, sep = "\t", col.names = FALSE)

pdf(paste0(PLOT_PREFIX, "_clusters.pdf"))
DimPlot(object_cl, reduction = "umap")
dev.off()

pdf(paste0(PLOT_PREFIX, "_samples.pdf"))
DimPlot(object_cl, reduction = "umap", group.by = "sample.name")
dev.off()

pdf(paste0(PLOT_PREFIX, "_cellType-filtByCat.pdf"))
DimPlot(object_cl, reduction = "umap", group.by = "Integrated_05", cells = colnames(object_cl)[object_cl@meta.data$category_inferred == CAT])
dev.off()

dea <- FindAllMarkers(object_cl, test.use = "MAST")
write.table(dea, file = DEA, sep = "\t", quote = FALSE)

q()



