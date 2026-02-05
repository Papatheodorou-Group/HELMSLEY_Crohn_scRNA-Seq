# from tutorial https://www.nature.com/articles/s41596-024-01045-4#Sec22
# N.B.: skip the optional steps

# devtools::install_github("jinworks/CellChat")
# devtools::install_github('immunogenomics/presto')
# reticulate::py_install(packages = 'umap-learn')
#
# module purge
# module load r-4.2.2-gcc-11.2.0-oa3uudy gmp-6.2.1-gcc-11.2.0-mneucsf 
# R_LIBS_USER=/hps/software/users/marioni/francesca/R_libs
# export R_LIBS_USER

# conda activate cellchat-env -> NOT NEEDED!

library("NMF")
library("circlize")
library("ComplexHeatmap")
library("CellChat")
library("patchwork")
future::plan("multisession", workers = 1)

dir.create("cellchat", showWarnings = FALSE)

CELL_GROUPS <- c("T cells", "B cells", "Endothelial", "Macrophages", "Mesenchymal")

object_normal <- readRDS("merged/merged_N_filt_2ndRound_TIL/hvg_pca_clust/object.Rds")
object_crohn <- readRDS("merged/merged_C_filt_no010/hvg_pca_clust/object.Rds")

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

CellChatDB <- CellChatDB.human
# showDatabaseCategory(CellChatDB)
# dplyr::glimpse(CellChatDB$interaction) 
CellChatDB.use <- subsetDB(CellChatDB) # exclude non-protein interactions

# run CellChat for each condition
cellchat_normal <- GenerateCellChatObject(object = object_normal, db = CellChatDB.use)
cellchat_crohn <- GenerateCellChatObject(object = object_crohn, db = CellChatDB.use)

saveRDS(cellchat_normal, file = "cellchat/cellchat_normal.rds") 
saveRDS(cellchat_crohn, file = "cellchat/cellchat_crohn.rds") 

# run differential CCI analysis
object.list <- list(N = cellchat_normal, C = cellchat_crohn)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))

save(object.list, file = "cellchat/cellchat_object.RData")
save(cellchat, file = "cellchat/cellchat_merged.RData") 



####################################


# module purge
# module load r-4.2.2-gcc-11.2.0-oa3uudy gmp-6.2.1-gcc-11.2.0-mneucsf 
# R_LIBS_USER=/hps/software/users/marioni/francesca/R_libs
# export R_LIBS_USER
#
# conda activate cellchat-env

library("NMF")
library("circlize")
library("ComplexHeatmap")
library("CellChat")
library("patchwork")
future::plan("multisession", workers = 1)

CELL_GROUPS <- c("T cells", "B cells", "Endothelial", "Macrophages", "Mesenchymal")

load("cellchat/cellchat_object.RData")
load("cellchat/cellchat_merged.RData")

# Comparing the total number of interactions
gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2)) 
pdf("cellchat/gg1.pdf")
gg1
dev.off()
# Comparing the total interaction strength
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight")
pdf("cellchat/gg2.pdf")
gg2
dev.off()

pdf("cellchat/netVisual_diffInteraction.pdf")
netVisual_diffInteraction(cellchat, weight.scale = T) 
dev.off()

pdf("cellchat/netVisual_diffInteraction_weight.pdf")
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight")
dev.off()

pdf("cellchat/netVisual_heatmap.pdf")
netVisual_heatmap(cellchat) 
dev.off()

# ERROR: not loading umap properly - skip this
# cellchat <- computeNetSimilarityPairwise(cellchat, type = "functional")
# cellchat <- netEmbedding(cellchat, type = "functional")
# cellchat <- netClustering(cellchat, type = "functional")

pdf("cellchat/rankNet.pdf")
rankNet(cellchat, slot.name = "netP", mode = "comparison", measure = "weight", sources.use = NULL, targets.use = NULL, stacked = T, do.stat = FALSE)
dev.off()

pdf("cellchat/rankNet_stat.pdf")
rankNet(cellchat, mode = "comparison", measure = "weight", stacked = T, do.stat = TRUE)
dev.off()

for (x in CELL_GROUPS) {
    for (y in CELL_GROUPS) {
        pdf(paste0("cellchat/rankNet_stat_", x, "_", y, ".pdf"), height = 5)
        print(rankNet(cellchat, mode = "comparison", measure = "weight", stacked = T, do.stat = TRUE, sources.use = x, targets.use = y))
        dev.off()
    }
}

pdf("cellchat/rankNet_net_stat.pdf", height = 40)
rankNet(cellchat, slot.name = "net", mode = "comparison", measure = "weight", stacked = T, do.stat = TRUE)
dev.off()

for (x in CELL_GROUPS) {
    for (y in CELL_GROUPS) {
        pdf(paste0("cellchat/rankNet_net_stat_", x, "_", y, ".pdf"), height = 10)
        print(rankNet(cellchat, slot.name = "net", mode = "comparison", measure = "weight", stacked = T, do.stat = TRUE, sources.use = x, targets.use = y))
        dev.off()
    }
}

# to extract the LR relationship with pathways:
# cellchat@LR$N$LRsig$interaction_name is the complex
# cellchat@LR$N$LRsig$pathway_name is the name of the pathway
# cellchat@LR$N$LRsig$ligand is the name of the ligand in the interaction 


q()


