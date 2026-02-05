# refine the clustering with subcluster information

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
    cat("\nUsage: <obj> <cl_mode> <file_list> <out_cl_mode>\n\n")
    q()
}

OBJECT <- args[1]
CL_MODE <- args[2]
FILE_LIST <- args[3]
OUT_CL_MODE <- args[4]

library("Seurat")

object <- readRDS(OBJECT)
file_list <- read.table(FILE_LIST, sep = "\t")

meta <- as.character(object@meta.data[[CL_MODE]])

for (i in 1:nrow(file_list)) {

    id <- file_list[i,1]
    file <- file_list[i,2]
    data <- read.table(file, sep = "\t")
    
    idx <- which(meta == id)
    idx2 <- match(colnames(object)[idx], data[,1])
    meta[idx] <- paste(id, data[idx2,2], sep = "_")
#    meta[idx] <- factor(meta[idx], levels = unique(meta[idx]))
}

names(meta) <- colnames(object)

object <- AddMetaData(object, metadata = meta, col.name = OUT_CL_MODE)

saveRDS(object, file = OBJECT)

q()



