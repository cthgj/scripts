#!/usr/bin/env Rscript
args <- commandArgs(TRUE)
a<-read.table(args[1])
write(mean(a$V2),stdout())
