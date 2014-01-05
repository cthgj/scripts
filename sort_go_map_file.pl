#!/usr/bin/perl -w
############################################################################
## This script is usuful for dealing with huge files of the following type:
##       GO pair                    stats
##   GO:0008711_GO:2000960	13753_0_0_0_0
##   GO:0008711_GO:2000198	13753_8_10_4_0
##   GO:0008711_GO:0042250	13753_0_0_0_0
##   GO:0006634_GO:0008711	13753_0_1_0_0
##
## And turn it into a file like :
##
##   stats                      LIST of gopairs with these stats
## 13753_0_0_0_0	        GO:0008711_GO:2000960 	GO:0008711_GO:0042250 
## 13753_0_1_0_0                GO:0006634_GO:0008711
############################################################################

my %index;

while(<>){
    chomp;
    $index{$_}=$index{$_}||tell();
    print "ii : $index{$_}\n";
    open (OUT,">>out");   
    seek(OUT,$index{$_},0);
    print OUT "hihihoho $index{$_}\n";
    close(OUT);
    
}

