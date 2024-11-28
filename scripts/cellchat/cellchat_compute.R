# from tutorial https://www.nature.com/articles/s41596-024-01045-4#Sec22
# N.B.: skip the optional steps

# devtools::install_github("jinworks/CellChat")
# devtools::install_github('immunogenomics/presto')
# reticulate::py_install(packages = 'umap-learn')
#

CELL_GROUPS <- c("T cells", "B cells", "Endothelial", "Macrophages", "Mesenchymal")

GenerateCellChatObject <- function(object, db) {
 
    data.input <- object[["RNA"]]@data 
    Seurat::Idents(object) <- object@meta.data$category_final
    cells <- colnames(object)[!(is.na(Seurat::Idents(object)) | Seurat::Idents(object) == "Epithelial inflamed")]
    object <- subset(object, cells = cells) # exclude the NA category and the "Epithelial inflamed", since they are only found in normal
    labels <- Seurat::Idents(object)
    meta <- data.frame(labels = labels, row.names = names(labels))

    object <- Seurat::AddMetaData(object, metadata = object@meta.data$sample.name, col.name="samples")
    cellchat <- createCellChat(object = object, group.by = "ident", assay = "RNA")
    cellchat@DB <- db
    cellchat <- subsetData(cellchat)

    cellchat <- identifyOverExpressedGenes(cellchat)
    cellchat <- identifyOverExpressedInteractions(cellchat) 
    cellchat <- computeCommunProb(cellchat, type = "triMean", trim = NULL, raw.use = TRUE)
    cellchat <- filterCommunication(cellchat, min.cells = 10)

    cellchat <- computeCommunProbPathway(cellchat) 
    sources.use <- CELL_GROUPS
    targets.use <- CELL_GROUPS
    cellchat <- aggregateNet(cellchat, sources.use = sources.use, targets.use = targets.use) 

    return(cellchat)
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
    cat("\nUsage: <obj1> <obj2> <out_dir>\n")
    cat("\n<obj1>     Seurat object for control condition\n")
    cat("\n<obj2>     Seurat object for disease condition\n")
    cat("\n<out_dir>  where to put CellChat results\n\n")
    q()
}

OBJ1 <- args[1]
OBJ2 <- args[2]
OUT_DIR <- args[3]

library("NMF")
library("circlize")
library("ComplexHeatmap")
library("CellChat")
library("patchwork")
future::plan("multisession", workers = 1)

dir.create(OUT_DIR, showWarnings = FALSE)

object_normal <- readRDS(OBJ1)
object_crohn <- readRDS(OBJ2)

CellChatDB <- CellChatDB.human
# showDatabaseCategory(CellChatDB)
# dplyr::glimpse(CellChatDB$interaction) 
CellChatDB.use <- subsetDB(CellChatDB) # exclude non-protein interactions

# run CellChat for each condition
cellchat_normal <- GenerateCellChatObject(object = object_normal, db = CellChatDB.use)
cellchat_crohn <- GenerateCellChatObject(object = object_crohn, db = CellChatDB.use)

saveRDS(cellchat_normal, file = file.path(OUT_DIR, "cellchat_normal.rds"))
saveRDS(cellchat_crohn, file = file.path(OUT_DIR, "cellchat_crohn.rds")) 

# run differential CCI analysis
object.list <- list(N = cellchat_normal, C = cellchat_crohn)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))

save(object.list, file = file.path(OUT_DIR, "cellchat_object.RData"))
save(cellchat, file = file.path(OUT_DIR, "cellchat_merged.RData"))


q()


