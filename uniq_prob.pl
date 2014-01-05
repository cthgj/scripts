#!/usr/bin/perl -w
use strict;

my %k; 
my %l; 
while(<>){
    chomp; 
    my @a=split(/\t/);
    $k{$a[0] . "_" .  $a[1]} = $k{$a[1] . "_" .  $a[0]}=$a[2];
}

foreach my $gp(keys(%k)){
    my @b=split(/_/,$gp);
    unless ((defined($l{$b[1] . "_" .  $b[0]}))||
	    (defined($l{$gp})) ){
	print "$b[0]\t$b[1]\t$k{$gp}\n" ;
    } 
    $l{$b[1] . "_" .  $b[0]}=$l{$gp}=1;
}

