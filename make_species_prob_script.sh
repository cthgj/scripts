#!/bin/sh
## This silly little script creates the scripts that will calculate probabilities for different species.
## Relevant data will be placed in the ./local directory. See research/pogo/README. The first argument
## is the species name and the second (if present) is the output dir.

## Directory to save/search for output files
bdir="local";

species=$1;

inter=0;
## If there are no arguments, echo usage info
if [ -z $1 ];
then

    echo "This script takes a species' name as a first argument. The second argument is an optional"
    echo "output dir name. If \"inter\" is given a second argument, it runs the interactome version"
    echo "and interprets the THIRD argument as the output_dir";
    echo "";
    echo "${0##*/} <species> [OPTIONAL:inter] [OPTIONAL:output_dir]";
    exit;
fi
## If a second argument is passed, then that is
## either 'inter', in which case we should run the 
## interactome version
sec=$2;
if [ $2 ]
then
    if [ "$2" = "inter" ]
    then
 	inter=1;
    else
    	bdir=$2;
    fi
fi
if [ $3 ];
then
    bdir=$3;
fi

rscript="pvall <- function(x) {xxx     r<-phyper(x[4], x[2], x[1] - x[2], x[3], low=T)	xxx     return(r)xxx}xxxpvalh <- function(x) {xxx     r<-phyper(x[4], x[2], x[1] - x[2], x[3], low=F)	xxx     return(r)xxx}xxxcounttxt <- read.table(\\\"$bdir/FILENAME1\\\", header=F, nrows=LENGTH,sep=\\\"\\\t\\\", row=1)xxxpairs <- row.names(counttxt)xxxcount <- data.matrix(counttxt)xxxpvl<-apply(count,1,pvall)xxxpvh<-apply(count,1,pvalh)xxxpvlhm<-p.adjust(pvl, method =\\\"holm\\\" , n = LENGTH)   xxxpvlho<-p.adjust(pvl, method =\\\"hochberg\\\" , n = LENGTH)         xxxpvlbh<-p.adjust(pvl, method =\\\"BH\\\" , n = LENGTH)   xxxpvlby<-p.adjust(pvl, method =\\\"BY\\\" , n = LENGTH)xxxpvhhm<-p.adjust(pvh, method =\\\"holm\\\" , n = LENGTH)  xxxpvhho<-p.adjust(pvh, method =\\\"hochberg\\\" , n = LENGTH)  xxxpvhbh<-p.adjust(pvh, method =\\\"BH\\\" , n = LENGTH) xxxpvhby<-p.adjust(pvh, method =\\\"BY\\\" , n = LENGTH) xxxtable <- data.frame(pair=pairs, total=count[,1], c1 <- count[,2], c2 <- count[,3], common <- count[,4], pvall=pvl, pvalh=pvh, pvlc1=pvlhm, pvlc2=pvlho , pvlc3=pvlbh , pvlc4=pvlby , pvhc1=pvhhm, pvhc2=pvhho , pvhc3=pvhbh , pvhc4=pvhby) xxxwrite.table(table, append=FALSE, \\\"$bdir/FILENAME2\\\", quote=FALSE, sep=\\\"\\\t\\\", col=FALSE, row=FALSE)xxx";

if [ $inter = 0 ]
then  
    echo  "#!/bin/bash";
    
    ## Collect gopair stats
    echo "nice cross_onto_probabilities.pl -vG data/GO.terms_alt_ids -a data/gene_association.$species -g data/all_three.genealogy > $bdir/$species.gostats 2>$species.error  &&"; 
    
    ## Condense to avoid repeat calculations
    echo "nice perl -e 'my %k; while(<>){chomp; @a=split(/\\\t/); shift(@a); \$b=join(\"_\",@a); \$k{\$b}++;} map{@a=split(/_/,\$_); \$\"=\"\\\t\"; print \"\$_\\\t@a\\\n\"}keys(%k)' $bdir/$species.gostats > $bdir/$species.gostats.cond  &&" ;     
    
    ## Count the lines of the uncondensed file (for R's multi-testing)
    echo "lines=\`wc -l $bdir/$species.gostats | gawk '{print \$1}'\`;";
    ## prepare the R script
    echo "echo \"$rscript\" | sed -e \"s/FILENAME1/$species.gostats.cond/g\" -e \"s/FILENAME2/$species.prob.cond/g\" -e \"s/LENGTH/\$lines/g\"  |perl -ne 's/xxx/\\\n/g; print' > $species.R  &&"
    echo "nice R CMD BATCH --no-restore --no-save $species.R &&";
    echo "uncondense_condensed_probs.pl $bdir/$species.prob.cond $bdir/$species.gostats > $bdir/$species.prob &&"
    echo "echo \"All done\" >> $species.error ";
    

## If we are running the interactome version
else

    echo  "#!/bin/bash";
    
    ## Collect gopair stats
    echo "nice cross_onto_probabilities.pl -ivG data/GO.terms_alt_ids -a data/gene_association.$species -g data/all_three.genealogy -n $species.gr -m data/$species.map > $bdir/$species.inter.gostats 2>$species.inter.error  &&"; 
    
    ## Condense to avoid repeat calculations
    echo "nice perl -e 'my %k; while(<>){chomp; @a=split(/\\\t/); shift(@a); pop(@a); \$b=join(\"_\",@a); \$k{\$b}++;} map{@a=split(/_/,\$_); \$\"=\"\\\t\"; print \"\$_\\\t@a\\\n\"}keys(%k)' $bdir/$species.inter.gostats > $bdir/$species.inter.gostats.cond  &&" ; 
    
    ## Count the lines of the uncondensed file (for R's multi-testing)
    echo "lines=\`wc -l $bdir/$species.inter.gostats | gawk '{print \$1}'\` &&";

    ## prepare the R script
    echo "echo \"$rscript\" | sed -e \"s/FILENAME1/$species.inter.gostats.cond/g\" -e \"s/FILENAME2/$species.inter.prob.cond/g\" -e \"s/LENGTH/\$lines/g\"  |perl -ne 's/xxx/\\\n/g; print' > $species.inter.R &&" 
	## Run it
    echo "nice R CMD BATCH --no-restore --no-save $species.inter.R &&";
    echo "uncondense_condensed_probs.pl -i $bdir/$species.inter.prob.cond $bdir/$species.inter.gostats  > $bdir/$species.inter.prob &&"
    echo "echo \"All done\" >> $species.inter.error";


    # ## Check if R is still running for this species    
    # echo "while [ 1 ]; do";
    # echo "user=\`whoami\`;" ; 
    # echo "running=\`ps aux | grep \$user |grep inter| gawk '{print \$2,\$11}' | grep R | wc -l\`;"; 
    # ## If it isn't, concatente the output files into one, and then delete them.
    # echo "if [ \$running -eq 0 ]";
    # echo "then";
    # echo "echo \"All done\" >> $species.inter.error ";
    # echo "exit;";
    # echo "else";
    # echo "sleep 5 ;";
    # echo "fi"; 
    # echo "done";
    
fi
