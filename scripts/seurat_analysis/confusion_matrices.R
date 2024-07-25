
DOUBLET_META <- "doublet_label_pca_vst_top5000" # Singlet and Doublet

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
    cat("\nUsage: <param> <stat_dir> <obj>\n")
    cat("\n<param>      file with clustering parameters\n")
    cat("<stat_dir>   output directory containing the statistics\n")
    cat("<obj>        object containing clustering solutions and other metadata\n")
    q()
}

PARAMS <- args[1]
STAT_DIR <- args[2]
OBJ <- args[3]

library("Seurat")

dir.create(STAT_DIR, recursive = TRUE, showWarnings = FALSE)

params <- as.matrix(read.table(PARAMS, sep="="))
disps <- as.numeric(unlist(strsplit(params[params[,1]=="disps", 2], split=",")))
nfeats <- as.numeric(unlist(strsplit(params[params[,1]=="nfeats", 2], split=",")))
k <- as.numeric(unlist(strsplit(params[params[,1]=="k", 2], split=",")))
res <- as.numeric(unlist(strsplit(params[params[,1]=="res", 2], split=",")))

object <- readRDS(OBJ)

# category / category_inferred confusion matrix
df <- table(object@meta.data[,c("category","category_inferred")])
TABLE <- file.path(STAT_DIR, "category-category_inferred.tsv")
write.table(df, file = TABLE, sep = "\t", quote = FALSE)

# category / doublet frac vector
for (meta in c("category", "Integrated_05")) {
    meta_lab <- sort(unique(object@meta.data[,meta]))
    v <- lapply(meta_lab, function(x) { 
                            y <- object@meta.data[,meta] == x 
                            z <- object@meta.data[,DOUBLET_META] == "Doublet"
                            sum(y & z) / sum(y)
                            } )
    v <- unlist(v)
    df <- data.frame(x = meta_lab, y = v)
    colnames(df) <- c(meta, "doublet_frac")
    TABLE <- file.path(STAT_DIR, paste0(meta, "-doublet_frac.tsv"))
    write.table(df, file = TABLE, sep = "\t", row.names = FALSE, quote = FALSE)
}

cases <- c()
for (ymin in disps)
    cases <- c(cases, paste("mean.var.plot_disp", ymin, sep=""))
for (nfeat in nfeats)
    cases <- c(cases, paste("vst_top", nfeat, sep=""))

for (case in cases) {
    dr <- paste("pca", case, sep="_")
    for (kk in k) {
        for (r in res) {
            cl_ident <- paste("clusters_", dr, "_k", kk, "_res", r, sep="")
            cl_labels <- object@meta.data[[cl_ident]]
            cl <- sort(unique(cl_labels))
            
            # cluster / meta matrices
            for (meta in c("category", "Integrated_05")) {
                df <- table(object@meta.data[,c(cl_ident,meta)])
                TABLE <- file.path(STAT_DIR, paste0(cl_ident, "-", meta, ".tsv"))
                write.table(df, file = TABLE, sep = "\t", quote = FALSE)
            }
            
            v <- lapply(cl, function(x) { 
                                     y <- object@meta.data[,cl_ident] == x 
                                     z <- object@meta.data[,DOUBLET_META] == "Doublet"
                                     sum(y & z) / sum(y)
                                     } )
            v <- unlist(v)
            df <- data.frame(x = cl, y = v)
            colnames(df) <- c(cl_ident, "doublet_frac")
            TABLE <- file.path(STAT_DIR, paste0(cl_ident, "-doublet_frac.tsv"))
            write.table(df, file = TABLE, sep = "\t", row.names = FALSE, quote = FALSE)
        }
    }
}



q()

