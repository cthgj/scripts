#!/usr/local/bin/perl

## This script should be run on the output of this command:
## parseblast blast.out | gawk '$14=="+"{print $12,$13"\t",$18,$19}' | sort -nk1
use strict;

my %pos;

while(<>)
{
   
    my @a = split;
#   my $mm= $. + $a[0]; 
#    $pos{$mm}=[@a];#push @{$pos{$a[0]}}, $a[1], $a[2], $a[3];
    push @{$pos{$a[0]}}, [@a];
}


foreach my $position (keys(%pos))
{
    my $start = 

    print "@{@{$pos{$position}}[0]}\n";




}
