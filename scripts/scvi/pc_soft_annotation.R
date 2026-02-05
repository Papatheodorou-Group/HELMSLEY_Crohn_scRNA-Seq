args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
    cat("Usage: <soft_q_lab> <GCA_stats> <out_matrix>\n")
    q()
}
SOFT_Q_LAB <- args[1]
GCA_STATS <- args[2]
OUT_PREFIX <- args[3]

library("ggplot2")

data <- read.table(SOFT_Q_LAB, sep = "\t")
category <- read.table(GCA_STATS, sep = "\t", header = TRUE)
cats <- unique(category[,1])

colnames(data) <- data[1,]
rownames(data) <- data[,1]
data <- data[2:nrow(data),2:ncol(data)]

for (ca in cats) {

    idx <- match(category[category[,1] == ca, 2], colnames(data))
    
    if (length(idx) > 1) {
    
        data_cat <- data[,idx]
        data_cat <- apply(data_cat, 2, as.numeric)
        data_cat <- data_cat[rowSums(data_cat) > 0.5,,drop = FALSE]
        data_cat <- data_cat[,colSums(data_cat) > 20,drop = FALSE]

        if (nrow(data_cat) > 1 & ncol(data_cat) > 1) {
            entropy <- apply(data_cat, 1, function(x) -sum(log(x,base=length(x))*x))
            fuzzy_AND_top2 <- apply(data_cat, 1, function(x) min(sort(x)[c(length(x)-1,length(x))]))
            pc <- princomp(data_cat)
    
            df <- data.frame(x = c(pc$scores[,1]), y = c(pc$scores[,2]), e = entropy, f = fuzzy_AND_top2)

            df_all <- as.data.frame(matrix(NA, nrow=0, ncol = 5))
            colnames(df_all) <- c("x", "y", "e", "ct", "score")
            for (i in 1:ncol(data_cat)) {
                df_tmp <- cbind(df, ct = rep(colnames(data_cat)[i],nrow(data_cat)), score = data_cat[,i])
                df_all <- rbind(df_all, df_tmp)
            }
    
            g <- ggplot(data = df_all, aes(x = x, y = y, col = score)) + geom_point() + facet_wrap(~ ct, ncol = 3) + theme_classic()
            g_ent <- ggplot(data = df, aes(x = x, y = y, col = e)) + geom_point() + theme_classic()
            g_fuzzyANDtop2 <- ggplot(data = df, aes(x = x, y = y, col = f)) + geom_point() + theme_classic()
            
            pdf(paste0(OUT_PREFIX, "_", ca, "_score.pdf"), width = 9, height = 3*ceiling(ncol(data_cat)/3))
            print(g)
            dev.off()

            pdf(paste0(OUT_PREFIX, "_", ca, "_ent.pdf"), width = 4.5, height = 4)
            print(g_ent)
            dev.off()
            
            pdf(paste0(OUT_PREFIX, "_", ca, "_fuzzyANDtop2.pdf"), width = 4.5, height = 4)
            print(g_fuzzyANDtop2)
            dev.off()
        }
    }
}

q()

