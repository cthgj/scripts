# count.txt:
#pair		       Total_proteins	Nprot_GO1  Nprot_GO2	Nprot_common	ratio
#GO:0006893_GO:0060429	   10935   	   8	      376     	     0		  0

pval <- function(x) {
     r<-phyper(x[4], x[2], x[1] - x[2], x[3], low=T)	
     return(r)
}
ease <- function(x) {
     r<-phyper(x[4]-1, x[2], x[1] - x[2], x[3], low=T)	
     return(r)
}


counttxt <- read.table("count_under.txt", header=F, row=1,sep="\t")
pairs <- row.names(counttxt)
count <- data.matrix(counttxt)

#pv<-apply(count,1,pval)
pv<-apply(count,1,ease)
	
library(multtest)
mt <- mt.rawp2adjp(pv, proc="BH")
adj <- mt$adjp[order(mt$index),]
adjpv <- adj[,2]
table <- data.frame(pair=pairs, total=count[,1], c1 <- count[,2], c2 <- count[,3], common <- count[,4], ovr <- count[,5], pval=pv, adjpv=adjpv)
no <- order(table$pval)
write.table(table[no,], "hyper_under.txt", quote=FALSE, sep="\t", col=FALSE, row=FALSE)



