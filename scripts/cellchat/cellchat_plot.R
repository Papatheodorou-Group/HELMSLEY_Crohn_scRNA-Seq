# from tutorial https://www.nature.com/articles/s41596-024-01045-4#Sec22
# N.B.: skip the optional steps

# devtools::install_github("jinworks/CellChat")
# devtools::install_github('immunogenomics/presto')
# reticulate::py_install(packages = 'umap-learn')
#

CELL_GROUPS <- c("T cells", "B cells", "Endothelial", "Macrophages", "Mesenchymal")

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
    cat("\nUsage: <out_dir>\n")
    cat("\n<out_dir>  where to put CellChat results\n\n")
    q()
}

OUT_DIR <- args[1]

library("NMF")
library("circlize")
library("ComplexHeatmap")
library("CellChat")
library("patchwork")
future::plan("multisession", workers = 1)

load(file.path(OUT_DIR, "cellchat_object.RData"))
load(OUT_DIR, "cellchat_merged.RData"))

# Comparing the total number of interactions
pdf(file.path(OUT_DIR, "gg1.pdf"))
compareInteractions(cellchat, show.legend = F, group = c(1,2)) 
dev.off()
# Comparing the total interaction strength
pdf(file.path(OUT_DIR, "gg2.pdf"))
compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight")
dev.off()

pdf(file.path(OUT_DIR, "netVisual_diffInteraction.pdf"))
netVisual_diffInteraction(cellchat, weight.scale = T) 
dev.off()

pdf(file.path(OUT_DIR, "netVisual_diffInteraction_weight.pdf"))
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight")
dev.off()

pdf(file.path(OUT_DIR, "netVisual_heatmap.pdf"))
netVisual_heatmap(cellchat) 
dev.off()

# ERROR: not loading umap properly - skip this
# cellchat <- computeNetSimilarityPairwise(cellchat, type = "functional")
# cellchat <- netEmbedding(cellchat, type = "functional")
# cellchat <- netClustering(cellchat, type = "functional")

pdf(file.path(OUT_DIR, "rankNet.pdf"))
rankNet(cellchat, slot.name = "netP", mode = "comparison", measure = "weight", sources.use = NULL, targets.use = NULL, stacked = T, do.stat = FALSE)
dev.off()

pdf(file.path(OUT_DIR, "rankNet_stat.pdf"))
rankNet(cellchat, mode = "comparison", measure = "weight", stacked = T, do.stat = TRUE)
dev.off()

for (x in CELL_GROUPS) {
    for (y in CELL_GROUPS) {
        pdf(file.path(OUT_DIR, paste0("rankNet_stat_", x, "_", y, ".pdf")), height = 5)
        print(rankNet(cellchat, mode = "comparison", measure = "weight", stacked = T, do.stat = TRUE, sources.use = x, targets.use = y))
        dev.off()
    }
}

pdf(file.path(OUT_DIR,"rankNet_net_stat.pdf"), height = 40)
rankNet(cellchat, slot.name = "net", mode = "comparison", measure = "weight", stacked = T, do.stat = TRUE)
dev.off()

for (x in CELL_GROUPS) {
    for (y in CELL_GROUPS) {
        pdf(file.path(OUT_DIR, paste0("rankNet_net_stat_", x, "_", y, ".pdf")), height = 10)
        print(rankNet(cellchat, slot.name = "net", mode = "comparison", measure = "weight", stacked = T, do.stat = TRUE, sources.use = x, targets.use = y))
        dev.off()
    }
}

# to extract the LR relationship with pathways:
# cellchat@LR$N$LRsig$interaction_name is the complex
# cellchat@LR$N$LRsig$pathway_name is the name of the pathway
# cellchat@LR$N$LRsig$ligand is the name of the ligand in the interaction 


q()


