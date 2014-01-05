#!/usr/bin/perl -w

use strict;
my ($f,$t,$lf,$lt)=("bob","bob","bob","bob");
while(<>){
    next if /^From/;
    /^(.+?)\s+(.+?)\s/; 
    $f=$1;
    $t=$2; 
   
    open(L,">l"); 
    print L "$2\n";
#    system("uniprot_fix_ids.pl l >kk"); 
    my $a=`uniprot_fix_ids.pl l`;
    print "$f\t$a"; 
    $lt=$t;
    $lf=$f;

}
