#!/usr/bin/perl

use warnings;
use strict;

my %dupes;
my (%k,%l);
while(<>){
    chomp;
    my @a=split(/\t/,$_);
    push @{$k{$a[0]}},$a[1];
    push @{$l{$a[1]}},$a[0];
}
