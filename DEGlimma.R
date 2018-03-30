# ----------------------------
# common template for python script
# ----------------------------

#.libPaths("~/lib/R")
library(limma)
args <- commandArgs(TRUE)
eset <- read.table(args[1], header=TRUE, row.names=1,check.names=FALSE)
genes <- rownames(eset)

group_info.raw <- read.table(args[2], header=FALSE, row.names=1)
group_info <- as.vector(group_info.raw$V2)
names(group_info) <- rownames(group_info.raw)

deg_info <- read.table(args[3], header=FALSE)

if (args[4] == 'absolute_binary'){
	isSigned <- FALSE
	isBinary <- TRUE
} else if (args[4] == 'signed_binary'){
	isSigned <- TRUE
	isBinary <- TRUE
} else if (args[4] == 'absolute_tstats'){
	isSigned <- FALSE
	isBinary <- FALSE
} else if (args[4] == 'signed_tstats'){
	isSigned <- TRUE
	isBinary <- FALSE
}

for (i in 1:nrow(deg_info)){
	group1 <- deg_info$V1[i]
	group2 <- deg_info$V2[i]
	eset.sub <- eset[1:nrow(eset),group_info[names(eset)] == group1 | group_info[names(eset)] == group2]
	design <- model.matrix(~factor(group_info[names(eset.sub)]))
	colnames(design) <- c(group1,group2)
	fit <- lmFit(eset.sub, design)
	fit <- eBayes(fit)
	options(digits=2)
	table <- topTable(fit, coef=2, adjust="fdr", n=Inf)
	#result <- table[genes, 5]
	table <- table[order(row.names(table)),]
	# columns of table 1:logFC, 2:AvgExpr 3:t 4:P.value 5:adj.P.val 6:B
	result <- table[1:nrow(table), 5]
	#print(rownames(table)[0:10])
	if (isBinary == TRUE){
		result <- as.integer(table[1:nrow(table),5] < 0.05)
		if (isSigned == TRUE){
			result <- result * sign(table[1:nrow(table),1])
		}
	}
	names(result) <- rownames(table)
	if (i==1){
		result.all <- cbind(result)
	} else{
		result.all <- cbind(result.all, result)
	}
	colnames(result.all)[ncol(result.all)] <- paste(group1,",",group2, sep="")
}

#filename="./DEG_result.T1.txt"
filename=args[5]
write.table(t(c("gene", colnames(result.all))), file=filename, sep="\t", col.names=FALSE, row.names=FALSE, quote=FALSE)
write.table(result.all, file=filename, append=TRUE, sep="\t", col.names=FALSE, row.names=TRUE, quote=FALSE)
