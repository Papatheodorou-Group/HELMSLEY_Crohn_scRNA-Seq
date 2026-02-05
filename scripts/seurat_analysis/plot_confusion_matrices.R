
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
    cat("\nUsage: <param> <stat_dir> <GCA_stats>\n")
    cat("\n<param>      file with clustering parameters\n")
    cat("<stat_dir>   output directory containing the statistics\n")
    cat("<GCA_stats>  input table containing the cell type - category mapping\n")
    q()
}

PARAMS <- args[1]
STAT_DIR <- args[2]
GCA_STATS <- args[3]

library("ComplexHeatmap")

params <- as.matrix(read.table(PARAMS, sep="="))
disps <- as.numeric(unlist(strsplit(params[params[,1]=="disps", 2], split=",")))
nfeats <- as.numeric(unlist(strsplit(params[params[,1]=="nfeats", 2], split=",")))
k <- as.numeric(unlist(strsplit(params[params[,1]=="k", 2], split=",")))
res <- as.numeric(unlist(strsplit(params[params[,1]=="res", 2], split=",")))

# category / category_inferred confusion matrix
TABLE <- file.path(STAT_DIR, "category-category_inferred.tsv")
df <- read.table(TABLE, sep = "\t", header = TRUE)

# plot confusion matrix 
M <- as.matrix(df)
M <- M/rowSums(M)
colnames(M) <- gsub("\\."," ",colnames(M))
H <- Heatmap(M, col = c("white","blue"), name = "frac",
             row_title = "scANVI category", column_title = "inferred category",
             cluster_rows = FALSE, cluster_columns = FALSE)
PDF <- file.path(STAT_DIR, "category-category_inferred.pdf")
pdf(PDF, width = 4, height = 4)
draw(H)
dev.off()

# category / doublet frac vector
TABLE <- file.path(STAT_DIR, "category-doublet_frac.tsv")
anno_cat <- read.table(TABLE, sep = "\t", header = TRUE)
TABLE <- file.path(STAT_DIR, "Integrated_05-doublet_frac.tsv")
anno_cell_type <- read.table(TABLE, sep = "\t", header = TRUE)

df_cat_celltype <- read.table(GCA_STATS, sep = "\t", header = TRUE)
orig_celltype <- df_cat_celltype$Integrated_05
df_cat_celltype$Integrated_05 <- gsub("[[:punct:]]", ".", df_cat_celltype$Integrated_05)
df_cat_celltype$Integrated_05 <- gsub("\\."," ",df_cat_celltype$Integrated_05)

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
            
            # cluster / meta matrices
            TABLE <- file.path(STAT_DIR, paste0(cl_ident, "-category.tsv"))
            df_cat <- read.table(TABLE, sep = "\t", header = TRUE)
            TABLE <- file.path(STAT_DIR, paste0(cl_ident, "-Integrated_05.tsv"))
            df_cell_type <- read.table(TABLE, sep = "\t", header = TRUE)
           
            TABLE <- file.path(STAT_DIR, paste0(cl_ident, "-doublet_frac.tsv"))
            anno_cl <- read.table(TABLE, sep = "\t", header = TRUE)
            
            M <- as.matrix(df_cat)
            M <- M/rowSums(M)
            colnames(M) <- gsub("\\."," ",colnames(M))
            x_anno <- HeatmapAnnotation(doublet = anno_cl[,2], which = "row", 
                                        show_legend = TRUE, show_annotation_name = FALSE)
            y_anno <- HeatmapAnnotation(doublet = anno_cat[,2], which = "column", 
                                        show_legend = TRUE, show_annotation_name = FALSE)
            H <- Heatmap(M, col = c("lemonchiffon","blue"), name = "frac",
                         row_title = cl_ident, column_title = "scANVI category",
                         cluster_rows = FALSE, cluster_columns = FALSE,
                         right_annotation = x_anno, bottom_annotation = y_anno)
            PDF <- file.path(STAT_DIR, paste0(cl_ident, "-category.pdf"))
            pdf(PDF, width = 4.5, height = 5)
            draw(H)
            dev.off()
            
            M <- as.matrix(df_cell_type)
            M <- M/rowSums(M)
            colnames(M) <- gsub("\\."," ",colnames(M))
            idx <- match(colnames(M), df_cat_celltype$Integrated_05)
            split <- df_cat_celltype$category[idx]
            colnames(M) <- orig_celltype[idx]
            x_anno <- HeatmapAnnotation(doublet = anno_cl[,2], which = "row", 
                                        show_legend = TRUE, show_annotation_name = FALSE)
            y_anno <- HeatmapAnnotation(doublet = anno_cell_type[,2], which = "column", 
                                        show_legend = TRUE, show_annotation_name = FALSE)
            H <- Heatmap(M, col = c("lemonchiffon","blue"), name = "frac", column_split = split,
                         row_title = cl_ident, column_title = sort(unique(split)),
                         cluster_rows = FALSE, cluster_columns = FALSE,
                         right_annotation = x_anno, bottom_annotation = y_anno)
            PDF <- file.path(STAT_DIR, paste0(cl_ident, "-cell_type.pdf"))
            pdf(PDF, width = 10, height = 7)
            draw(H)
            dev.off()
            
        }
    }
}



q()

