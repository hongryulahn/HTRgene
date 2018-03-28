library(skmeans)
args<-commandArgs(TRUE)

fname=args[1]
k=strtoi(args[2])
seed=args[3]
outfile=args[4]
if(seed!="None"){
	set.seed(strtoi(seed))
}
gene_expr<-read.table(fname, header=TRUE, row.names=1, check.names=FALSE)
gene_exprmat=as.matrix(gene_expr)

cl<-skmeans(gene_exprmat, k, method="genetic", m=1, weights=1, control=list(start=gene_exprmat[1:k,]))
write.table(cl$cluster, file=outfile, sep="\t", quote=FALSE, col.names = FALSE)
#write.table(cl$prototypes, file=outprototype, sep="\t", quote=FALSE)
