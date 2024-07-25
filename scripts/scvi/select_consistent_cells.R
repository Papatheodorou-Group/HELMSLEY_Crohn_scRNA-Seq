args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
    cat("Usage: <stat_q> <stat_r> <cells>\n")
    q()
}
STAT_Q <- args[1]
STAT_R <- args[2]
CELLS <- args[3]

META <- c("category","Integrated_05")

stat_q <- read.table(STAT_Q, sep = "\t", header = TRUE)
stat_r <- read.table(STAT_R, sep = "\t", header = TRUE)

ref <- unique(stat_r[,META])
idx <- which(apply(stat_q[,META], 1, function(x) { sum( apply( ref, 1, function(y) { sum(x == y) == 2 } )) == 1 } ) )

nms <- colnames(stat_q)
if (sum(nms == paste0(META[1], "_inferred")) == 0) {
    new_col <- apply(stat_q[,META], 1, function(x) ref[which(ref[,META[2]] == x[2]), META[1]] )
    stat_q <- cbind(stat_q, new_col)
    colnames(stat_q) <- c(nms, paste0(META[1], "_inferred"))
    write.table(stat_q, file = STAT_Q, sep = "\t", quote = FALSE, row.names = FALSE) 
}

cells <- stat_q[idx,"X"]
write(cells, file = CELLS)



q()


