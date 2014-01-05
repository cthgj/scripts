#!/usr/bin/perl -w

use strict;

my @P=();




print<<HAhaHA;
library(limma)
modea<-scan("a.names",'character')
modeb<-scan("b.names",'character')
modeB<-scan("B.names",'character')
modei<-scan("i.names",'character')
modeP<-scan("P.names",'character')
HAhaHA

my @c=("blue","red","gold","green","magenta");
my @a=("a", "b", "B", "i","P");
for (my $k=0;$k<=$#a;$k++){
    for (my $kk=$k+1;$kk<=$#a;$kk++){
	for (my $kkk=$kk+1;$kkk<=$#a;$kkk++){
	    
	    print<<HOhoHO
png(filename="$a[$k]$a[$kk]$a[$kkk].png")
all <-sort ( unique( c(mode$a[$k],mode$a[$kk],mode$a[$kkk])))
Counts <- matrix(0, nrow=length(all), ncol=3)

colnames(Counts) <- c(paste("$a[$k] (",length(mode$a[$k]),")"), paste("$a[$kk] (",length(mode$a[$kk]),")"), paste("$a[$kkk] (",length(mode$a[$kkk]),")"))
cols<-c("$c[$k]","$c[$kk]","$c[$kkk]")
for (i in 1:length(all))
{
Counts[i,1] <- all[i] %in% mode$a[$k]
Counts[i,2] <- all[i] %in% mode$a[$kk]
Counts[i,3] <- all[i] %in% mode$a[$kkk]
}
vennDiagram( vennCounts(Counts),circle.col=cols )
HOhoHO
	}
    }
}


