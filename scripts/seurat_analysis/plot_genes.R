# merge the doublet score by doubletFinder and the known doublet cluster to define the final set of doublets
# the final number of doublet will match 10x guidelines, given the number of loaded cells

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
    cat("\nUsage: <dir> <obj> <genes>\n")
    q()
}

DIR <- args[1]
OBJ <- args[2]
GENES <- args[3]

library("Seurat")

obj <- readRDS(OBJ)
data <- read.table(GENES, sep = ";")[,1]

for (d in data) {
    v <- unlist(strsplit(d, split = "\t"))
    category <- v[1]
    cell_type <- v[2]
    genes <- v[3:length(v)]
    genes <- genes[genes %in% rownames(obj)]
    
    if (length(genes) > 0) {
        pdf(file.path(DIR, paste0("UMAP_",category,"_",cell_type,".pdf")))
        print(FeaturePlot(obj, reduction="umap", features = genes))
        dev.off()
    }
}


q()

