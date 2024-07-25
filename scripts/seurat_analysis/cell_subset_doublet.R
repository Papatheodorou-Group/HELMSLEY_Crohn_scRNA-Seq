# merge the doublet score by doubletFinder and the known doublet cluster to define the final set of doublets
# the final number of doublet will match 10x guidelines, given the number of loaded cells

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
    cat("\nUsage: <object2> <dr> <cells> <object1> <mode> <cl>\n")
    cat("\n<object1>        input object, containing <mode> metadata\n")
    cat("<mode>          cluster mode as stored in metadata\n")
    cat("<cl>            label of the cluster(s) to be filtered out\n")
    cat("<object2>       containing the doublet score and label\n")
    cat("<dr>            dimansional reduction for object2 (where doublets are computed)\n")
    cat("<cells>         output list of cells retained\n\n")
    q()
}

OBJ2 <- args[1]
DR <- args[2]
CELLS <- args[3]
if (length(args) > 5) {
    OBJ1 <- args[4]
    MODE <- args[5]
    CL <- args[6]
}

library("Seurat")

if (length(args) > 5) { # optionally apply by-cluster filtering
    cls <- unlist(strsplit(CL, split=","))
    object <- readRDS(OBJ1)
    idx <- which(!(object@meta.data[[MODE]] %in% cls))
    cells <- colnames(object)[idx]
} 

object <- readRDS(OBJ2)
doublet_score <- paste0("doublet_score_", DR)
doublet_label <- paste0("doublet_label_", DR)

if (length(args) < 6) {
    cells <- colnames(object)
}

doublet_scores <- object@meta.data[[doublet_score]]
total_doublets <- sum(object@meta.data[[doublet_label]] == "Doublet")
doublets_filt_by_cluster <- sum(!(colnames(object) %in% cells))
doublets_to_filter <- total_doublets - doublets_filt_by_cluster

idx_filtered <- which(colnames(object) %in% cells)
threshold <- sort(doublet_scores[idx_filtered], decreasing = TRUE)[doublets_to_filter]
idx_filtered_2 <- which(doublet_scores[idx_filtered] < threshold)

cells <- colnames(object)[idx_filtered[idx_filtered_2]]

write.table(cells, file = CELLS, col.names = FALSE, row.names = FALSE, quote = FALSE)


q()

