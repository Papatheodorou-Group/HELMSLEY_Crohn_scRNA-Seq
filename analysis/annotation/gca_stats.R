df_cat <- read.table("GCA/adata/adult_pediatric_category.tsv", sep = "\t", header = TRUE)
df_ct <- read.table("GCA/adata/adult_pediatric_Integrated_05.tsv", sep = "\t", header = TRUE)
df <- data.frame(df_cat$labels_scanvi, df_ct$labels_scanvi)
df_out <- table(df)
tab <- as.data.frame(df_out)
colnames(tab) <- c("category", "Integrated_05", "count")
tab <- tab[tab$count > 0,]
tab <- tab[order(tab$category),]
tab
savehistory("gca_stats.R")
