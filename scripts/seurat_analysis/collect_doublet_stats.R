
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
    cat("\nUsage: <cl_obj> <slot> <doubl_obj> <out.tsv>\n")
    cat("\n<cl_obj>      original object with clustering output\n")
    cat("<slot>       clustering slot in <cl_obj>\n")
    cat("<doubl_obj>  object containing doublet results from doubletfinder\n")
    cat("<out.tsv>    output statistics\n\n")
    q()
}

CL_OBJ <- args[1]
SLOT <- args[2]
DOUBL_OBJ <- args[3]
OUT_PREFIX <- args[4]

library("Seurat")

object <- readRDS(CL_OBJ)
cl_all <- object@meta.data[[SLOT]]
names(cl_all) <- colnames(object)

object <- readRDS(DOUBL_OBJ)
cells <- colnames(object)

dr <- gsub("_k.*", "", gsub("clusters_pca_", "", SLOT))
dr <- paste0("pca_", dr)
doublet_label_slot <- paste0("doublet_label_", dr)
doublet_score_slot <- paste0("doublet_score_", dr)

df <- object@meta.data[,match(c(doublet_label_slot, doublet_score_slot),colnames(object@meta.data))]
rownames(df) <- colnames(object)
colnames(df) <- c("doublet_label", "doublet_score")

OUT_TSV <- paste0(OUT_PREFIX, "_res.tsv")
write.table(df, file = OUT_TSV, sep = "\t", quote = FALSE)

cl <- cl_all[match(cells, names(cl_all))]
cl <- factor(cl, levels = sort(unique(cl)))

tab <- table(data.frame(clusters = cl, doublets = df[,1]))
doublet.frac <- tab[,"Doublet"] / rowSums(tab)
df <- data.frame(clusters = rownames(tab), doublets = tab[,"Doublet"], doublet.frac = doublet.frac)

OUT_TSV <- paste0(OUT_PREFIX, "_stats.tsv")
write.table(df, file = OUT_TSV, sep = "\t", row.names = FALSE, quote = FALSE)


q()

