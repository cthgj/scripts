#!/usr/bin/perl -w
## This script will take a condensed prob file:
## num_annot_to_onto\tnum_annot_togo1\tnum_annot_togo2\tnum_annot_toboth\tprob(s)
## e.g.: 25736_227_83_5	25736	227	83	5	0.94595922589822
##
## and the uncodensed gostats file:
## GO:0001664_GO:0021782	25736	227	83	5
## and create an uncondensed prob file:
## 
## GO_pair\tgo1\tgo2\tstats\tprobs
use strict;
use Getopt::Std;

my %opts;
getopts('vi',\%opts);
unless($ARGV[1]){
    print STDERR "USAGE: $0 [-i] condensed.prob uncondensed.gostats\n\n\t-i for interactome, assumes that the last column is the number of self interactions\n";
    exit;
}
my $inter=$opts{i}||undef;

#####################################
# Open condensed probabilities file #
#####################################
open(A,"$ARGV[0]"); 
my %k; 
while(<A>){
    print STDERR "$.\r"  if $opts{v};
    chomp; 
    my @a=split(/\t/); 
    my $id=shift(@a); 
    $k{$id}=join("\t",@a);
} 
close(A);
print STDERR "Done\n" if $opts{v}; 
#####################################
# Open uncondensed gostats file     #
#####################################
open(B,"$ARGV[1]") or die "Could not open argv1 : $ARGV[1] : $!\n";
my $same="";
$inter ? 
    print "#GO\tTotal\tgo1\tgo2\t#go1\t#go2\tboth\tpval-low\tpval-high\tpval-low(holm)\tpval-low(hochb.)\tpval-low(BH)\tpval-low(BY)\tpval-high(holm)\tpval-high(hochb.)\tpval-high(BH)\tpval-high(BY)\tsame\n" :
    print "#GO\tTotal\tgo1\tgo2\t#go1\t#go2\tgo1\tgo2\tboth\tpval-low\tpval-high\tpval-low(holm)\tpval-low(hochb.)\tpval-low(BH)\tpval-low(BY)\tpval-high(holm)\tpval-high(hochb.)\tpval-high(BH)\tpval-high(BY)\n";

while(<B>){
    chomp; 
    my @a=split(/\t/);
    my $go=shift(@a); 
    my @b=sort {$b lt $a} split("_",$go);
    $same="\t" . pop(@a) if $inter; 
    my $id=join("_",@a); 
    print "$go\t$b[0]\t$b[1]\t$k{$id}\t$same\n";
}
close(B);


