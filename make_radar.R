#!/usr/bin/env Rscript
args <- commandArgs(TRUE)
library(fmsb)
name<-args[3]
df1<-read.table(args[1],header=T,row.names=1)
png(filename=args[2], width = 800, height = 800, units = "px", pointsize = 18, bg = "white")
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "magenta4", "#0072B2", "lightslategrey", "#CC79A7")
# cbbPalette<-terrain.colors(nrow(df1)-2)
# cbbPalette<- c("dodgerblue2",
# 	     "#E31A1C", # red
#                 "green4",
#                 "#6A3D9A", # purple
#                 "#FF7F00", # orange
#                 "black",
# 		"gold1",
#                 "skyblue2",
# 		"#FB9A99", # lt pink
#                 "palegreen2",
#                 "#CAB2D6", # lt purple
#                 "#FDBF6F", # lt orange
#                 "gray70", 
# 		"khaki2",
#                 "maroon",
# 		"orchid1",
# 		"deeppink1",
# 		"blue1",
# 		"steelblue4",
#                 "darkturquoise",
# 		"green1",
# 		"yellow4",
# 		"yellow3",
#                 "darkorange4",
# 		"brown")
leg.txt<-c("Candidates",name,"Network average")
radarchart(df1,seg=4,maxmin=T,pcol=cbbPalette,plty=1,pty=16,plwd=2,cglty=1,cglwd=0.15,centerzero=F)
legend(0.75,1.28,leg.txt,pch =18, col=cbbPalette,bty="n")
dev.off()
