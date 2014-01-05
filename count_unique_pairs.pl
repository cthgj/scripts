#!/usr/bin/perl -w
use strict;
use Getopt::Std;
my (%opts);
getopts('vpts:c:',\%opts);
if ($opts{c}){
    my $ll=($opts{c} * ($opts{c}-1))/2;
    print "TT : $ll\n";
    exit;
}
die("Need a geneology file (\$ARGV[0])\n") unless $ARGV[0];

my $infile=$ARGV[0]||"biological_process.genealogy";
open(A,"$infile");
my @GOs;
while(<A>){
    chomp;
    my @a=split(/\t/);
    push @GOs,$a[0];
}

if($opts{p}){
    my $l=0;
    for (my $n=0; $n<=$#GOs;$n++){
	for (my $k=$n+1; $k<=$#GOs;$k++){
	    $l++;
	    print STDOUT "$GOs[$n] $GOs[$k]\n";
	}
	printf(STDERR "$n of " . ($#GOs+1) ."\r") if $n % 100 == 0 && $opts{v};
    }
    print STDERR scalar(@GOs) . " terms => $l total pairs\n";
}
else{
    print scalar(@GOs) . " terms => " . (($#GOs-1) * $#GOs)/2 . " unique pairs\n"

}
# ## not that Elegant!! crap...
# if ($opts{t}){
#     for(my $a=0; $a<scalar(@GOs);$a++){
# 	$tot=$tot+scalar(@GOs)-$a;
#     }
#     my $real_tot=$tot-scalar(@GOs);
#     print STDERR "$k GO terms => $real_tot unique pairs\n" ;
# }
