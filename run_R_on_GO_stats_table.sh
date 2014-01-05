#!/bin/bash
n=`wc -l $1 | gawk '{print $1}'`;
species=$2;
k=$(($n / 8))
echo "splitting..." >&2;
split -l $k $1 $species;
c=0;


for i in $(ls $species*);
do
    c=$((c+1));
    out=$i.out; 
    names[$c]=$i;
    cat <<EOF >> $i.R
pval <- function(x) {
     r<-phyper(x[4], x[2], x[1] - x[2], x[3], low=T)	
     return(r)
}
counttxt <- read.table("$i", header=F, row=1,sep="\t")
pairs <- row.names(counttxt)
count <- data.matrix(counttxt)

pv<-apply(count,1,pval)
	
table <- data.frame(pair=pairs, total=count[,1], c1 <- count[,2], c2 <- count[,3], common <- count[,4], ovr <- count[,5], pval=pv)
write.table(table, "$out", quote=FALSE, sep="\t", col=FALSE, row=FALSE)

EOF
    
done;
for i in $(seq 1 $((c-(c/2)-1)));
do
    echo -n "nice R CMD BATCH ${names[$i]}.R && "
done;
echo -n "nice R CMD BATCH  ${names[((c-(c/2)))]}.R &"
echo;
for i in $(seq $((c-(c/2)+1)) $((c-1)));
do
    echo -n "nice R CMD BATCH ${names[$i]}.R && "
done;
echo -n "nice R CMD BATCH ${names[$c]}.R &"


echo;
