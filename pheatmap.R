#!/usr/bin/env Rscript
args <- commandArgs(TRUE)
df<-read.table(args[1],header=T,row.names=1)
mat<-data.matrix(df)
## Get the colors
library(RColorBrewer)
blues<-colorRampPalette(brewer.pal(9,"Blues"))(256)
library(pheatmap)
library(grid)
## This will make the names diagonal, see http://stackoverflow.com/a/15506652/1081936
draw_colnames_45 <- function (coln, ...) {
     m = length(coln)
     x = (1:m)/m - 1/2/m
     grid.text(coln, x = x, y = unit(0.96, "npc"), vjust = .5, hjust = 1, rot = 45, gp = gpar(...)) ## Was 'hjust=0' and 'rot=270'
}
assignInNamespace(x="draw_colnames", value="draw_colnames_45",ns=asNamespace("pheatmap"))
pdf(args[3],points=8,family="NewCenturySchoolbook"); 
pheatmap(mat,main=args[2],display_numbers=F,color=blues,fontsize=10,fontsize_col=5,fontsize_row=5);
dev.off()
