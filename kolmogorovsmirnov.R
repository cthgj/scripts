#!/usr/bin/env Rscript
args <- commandArgs(TRUE)
a<-read.table(args[1],sep="\t",col.names=c("Name","val"))
b<-read.table(args[2],sep="\t",col.names=c("Name","val"))
ks<-ks.test(a$val, b$val)
write(paste(format.pval(ks$p.value,digits=4),format.pval(ks$statistic,digits=4),sep="\t"),stdout()) 
