
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
    cat("\nUsage: <in_obj> <out_obj> <sample_info>\n")
    q()
}

IN_OBJ <- args[1]
OUT_OBJ <- args[2]
SAMPLE_INFO <- args[3]

library("Seurat")
library("dplyr")

obj <- readRDS(IN_OBJ)
sample_info <- read.table(SAMPLE_INFO, sep = "\t")

idx <- which(sample_info[,1] %in% obj@meta.data$sample.name)
sample_info <- sample_info[idx,]

sample_name <- sample_info[,1]
diagnosis <- sample_info[,2]
region_code <- sample_info[,3] 

obj@meta.data <- dplyr::rename(obj@meta.data, "Sample name" = sample.name)
meta_diag <- plyr::mapvalues(obj@meta.data[["Sample name"]], from = sample_name, to = diagnosis)
meta_reg <- plyr::mapvalues(obj@meta.data[["Sample name"]], from = sample_name, to = region_code)

df <- data.frame(meta_diag, meta_reg)
rownames(df) <- colnames(obj)
obj <- AddMetaData(object = obj, metadata = df)

obj@meta.data <- dplyr::rename(obj@meta.data, "Diagnosis" = "meta_diag")
obj@meta.data <- dplyr::rename(obj@meta.data, "Region code" = "meta_reg")

saveRDS(obj, file = OUT_OBJ)


q()

