# re-label the clusters

CL_MODE_CT <- "cell_type_final"
CL_MODE_CAT <- "category_final"

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5) {
    cat("\nUsage: <obj> <out_cl_mode> <annotation> <plot_prefix> <reduction>\n\n")
    q()
}

OBJECT <- args[1]
OUT_CL_MODE <- args[2]
ANNOTATION <- args[3]
PLOT_PREFIX <- args[4]
REDUCTION <- args[5]

library("Seurat")

object <- readRDS(OBJECT)
annotation <- read.table(ANNOTATION, sep = "\t")

ncells <- ncol(object)
meta <- as.character(object@meta.data[[OUT_CL_MODE]])

out_meta_cat <- out_meta_ct <- rep(NA,ncells)
for (i in 1:nrow(annotation)) {
    idx <- which(meta == annotation[i,1])
    out_meta_ct[idx] <- rep(annotation[i,2],length(idx))
    out_meta_cat[idx] <- rep(annotation[i,3],length(idx))
}
names(out_meta_ct) <- names(out_meta_cat) <- colnames(object)
object <- AddMetaData(object, metadata = out_meta_ct, col.name = CL_MODE_CT)
object <- AddMetaData(object, metadata = out_meta_cat, col.name = CL_MODE_CAT)

# plot cell types

Idents(object) <- object@meta.data[[CL_MODE_CT]]
cells <- colnames(object)[!is.na(object@meta.data[[CL_MODE_CT]])]

pdf(paste0(PLOT_PREFIX, "_", CL_MODE_CT, ".pdf"), width = 12)
DimPlot(object, reduction = REDUCTION, cells = cells, label = TRUE, repel = TRUE)
dev.off()

# plot categories

object@meta.data[[CL_MODE_CAT]] <- factor(object@meta.data[[CL_MODE_CAT]], levels = sort(unique(object@meta.data[[CL_MODE_CAT]])))
Idents(object) <- object@meta.data[[CL_MODE_CAT]]
cells <- colnames(object)[!is.na(object@meta.data[[CL_MODE_CAT]])]

pdf(paste0(PLOT_PREFIX, "_", CL_MODE_CAT, ".pdf"), width = 12)
DimPlot(object, reduction = REDUCTION, cells = cells, label = TRUE, repel = TRUE)
dev.off()

saveRDS(object, file = OBJECT)

q()



