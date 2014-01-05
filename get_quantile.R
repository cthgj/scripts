#!/usr/bin/env Rscript
args <- commandArgs(TRUE)
a<-read.table(args[1])
write(quantile(a$V2,c(.80)),stdout())
