
# NB: here assume that the cells contain the sample name in their ID

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
    cat("\nUsage: <in_meta> <sample_name> <out_meta> <meta_col>\n")
    q()
}

IN_META=args[1]
SAMPLE_NAME=args[2]
OUT_META=args[3]
META_COL=args[4]

meta <- unlist(strsplit(META_COL, split = ","))
SAMPLE_NAME <- gsub(",","\\|",SAMPLE_NAME) # select cells in any sample

df <- read.table(IN_META, sep = "\t", header = TRUE)

rownames(df) <- df[,1]
df <- df[grep(SAMPLE_NAME, rownames(df)),]
df <- df[,meta]

write.table(df, file = OUT_META, sep = "\t", quote = FALSE)

q()

