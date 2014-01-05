pval <- function(x) {
     r<-phyper(x[4], x[2], x[1] - x[2], x[3], low=T)	
     return(r)
}
ease <- function(x) {
     if(x[4]-1 == 0){r<-phyper(x[4], x[2], x[1] - x[2], x[3], low=T)}	
     else{r<-phyper(x[4]-1, x[2], x[1] - x[2], x[3], low=T)}	
     return(r)
}


counttxt <- read.table("calculate_probabilities1.out", header=F, row=1,sep="\t")
pairs <- row.names(counttxt)
count <- data.matrix(counttxt)

#pv<-apply(count,1,pval)
pv<-apply(count,1,ease)
	
library(multtest)
mt <- mt.rawp2adjp(pv, proc="BH")
adj <- mt$adjp[order(mt$index),]
adjpv <- adj[,2]
table <- data.frame(pair=pairs, total=count[,1], c1 <- count[,2], c2 <- count[,3], common <- count[,4], pval=pv, adjpv=adjpv)
no <- order(table$pval)
write.table(table[no,], "hyper_under_mine.txt", quote=FALSE, sep="\t", col=FALSE, row=FALSE)



