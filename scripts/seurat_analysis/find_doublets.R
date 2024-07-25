EPSILON <- 0.0001
pN <- 0.25


ExtractPCs <- function(object, dr = "pca", pcs = 50, pval = 1e-05, non.rand.sd.frac = 0.5, file) {

    # compute the number of top significant PCs 
    # that explain 50% more variance compared to random (i.e. last 10 PCs)
    # this check is meant to exclude artifacts due to the svd approximation
    # that is used to compute the PCs
    min_stdev <- (1+non.rand.sd.frac)*mean(object@reductions[[dr]]@stdev[(pcs-10):pcs])
    for (nPCs in 0:(length(object@reductions[[dr]]@stdev)-1)) {
        if (object@reductions[[dr]]@stdev[nPCs+1] < min_stdev) 
            break
    }
    # take the highest PC that explains > 0.05 variance with respect to the next
    write(nPCs, file = file)

}

GetNPCs <- function(in_dir) {

    pc_file <- file.path(in_dir,"num_PCs.txt")
    nPCs <- drop(as.matrix(read.table(pc_file)))

    nPCs
}

DoubletFind <- function(object, dr = "pca", pcs = 50, expDoublPerc = 2, out_dir) {

    names(object@reductions)[names(object@reductions) == dr] <- "pca"
    
    pdf(file.path(out_dir, "doubletFinder.pdf"))
    
    sweep.res.list <- paramSweep_v3(object, PCs = 1:pcs, sct = FALSE)
    sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
    bcmvn <- find.pK(sweep.stats)
    optimal_pK <- bcmvn$pK[which(bcmvn$BCmetric > max(bcmvn$BCmetric) - EPSILON)]
    optimal_pK <- as.numeric(as.character(optimal_pK))
    
    expDoubl <- round(ncol(object)*expDoublPerc/100)
    object <- doubletFinder_v3(object, PCs = 1:pcs, pN = pN, pK = optimal_pK, nExp = expDoubl, reuse.pANN = FALSE, sct = FALSE)
    
    dev.off()
    
    id <- paste(pN, optimal_pK, expDoubl, sep = "_")
    score_slot <- paste0("pANN_", id)
    label_slot <- paste0("DF.classifications_", id)
    
    colnames(object@meta.data)[colnames(object@meta.data) == score_slot] <- paste0("doublet_score_", dr)
    colnames(object@meta.data)[colnames(object@meta.data) == label_slot] <- paste0("doublet_label_", dr)
    names(object@reductions)[names(object@reductions) == "pca"] <- dr
    
    return(object)
}

DrPlot <- function(object, dr = "pca", dr_plot = "pca", pcs = 50, out_dir, dr_type = "PCA") {

    doublet_score_slot <- paste0("doublet_score_", dr)
    doublet_label_slot <- paste0("doublet_label_", dr)

    doublets <- colnames(object)[object@meta.data[[doublet_label_slot]] == "Doublet"] 
    pdf(file.path(out_dir, paste0(dr_type, "_plot_doublet_label.pdf")), width = 3.7*2, height = 3*2)
    print(DimPlot(object, cells.highlight = doublets, reduction = dr_plot))
    dev.off()
    
    pdf(file.path(out_dir, paste0(dr_type, "_plot_doublet_score.pdf")), width = 3.7*2, height = 3*2)
    print(FeaturePlot(object, features = doublet_score_slot, reduction = dr_plot))
    dev.off()    

}


args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5) {
    cat("\nUsage: <obj.Robj> <in_dir> <out_dir> <params> <doublPerc>\n")
    cat("\n<obj.Robj>   input Seurat object created at step 0\n")
    cat("<in_dir>     input directory containing the list of genes (from previous step) and the results of PCs selection\n")
    cat("<out_dir>    output directory containing the plots\n")
    cat("<params>     file including the values for the parameters, separated by \"=\"\n")
    cat("<doublPerc>  expected percentage of doublets (1-100)\n\n")
    q()
}

OBJECT <- args[1]
IN_DIR <- args[2]
OUT_DIR <- args[3]
PARAMS <- args[4]
DOUBLET_PERC <- as.numeric(args[5])


params <- as.matrix(read.table(PARAMS, sep="="))
cores <- as.numeric(params[params[,1]=="cores", 2])
nfeats <- as.numeric(unlist(strsplit(params[params[,1]=="nfeats", 2], split=",")))

library("Seurat")
library("DoubletFinder")
library("ggplot2")
library("cowplot")
theme_set(theme_cowplot())

dir.create(OUT_DIR, showWarnings = FALSE)

object <- readRDS(OBJECT)

feature_method="vst"
for (nfeat in nfeats) {
    
    case <- paste("vst_top", nfeat, sep = "")
    dr <- paste("pca", case, sep = "_")
    in_dir <- file.path(IN_DIR, case)
    out_dir <- file.path(OUT_DIR, dr)
    dir.create(out_dir, showWarnings = FALSE)
    
    ExtractPCs(object, dr = dr, file = file.path(in_dir, "num_PCs.txt"))
    nPCs <- GetNPCs(in_dir)
    object <- DoubletFind(object, dr = dr, pcs = nPCs, expDoublPerc = DOUBLET_PERC, out_dir = out_dir)
    
    tsne_dr <- paste0("tsne_", dr)
    umap_dr <- paste0("umap_", dr)
    if (!(tsne_dr %in% names(object@reductions)))
        object <- RunTSNE(object, reduction = dr, dims = 1:nPCs, reduction.name = tsne_dr, reduction.key = tsne_dr)
    if (!(umap_dr %in% names(object@reductions)))
        object <- RunUMAP(object, reduction = dr, dims = 1:nPCs, reduction.name = umap_dr, reduction.key = umap_dr)
    DrPlot(object, dr = dr, dr_plot = tsne_dr, out_dir = out_dir, dr_type = "TSNE")
    DrPlot(object, dr = dr, dr_plot = umap_dr, out_dir = out_dir, dr_type = "UMAP")
}

saveRDS(object, file=OBJECT)


sessionInfo()
q()


