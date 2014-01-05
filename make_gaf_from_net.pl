#!/usr/bin/perl -w

use strict;
unless($ARGV[0]){
    print STDERR "USAGE: $0 NETWORK MAP GAF > subGAF\n"; exit;
}
open(N, $ARGV[0])||die("Need a net file \$ARGV[0]\n");
my (%prots);
while(<N>){
    chomp;
    my @a=split(/\t/);
    $prots{$a[0]}=$prots{$a[1]}=1;
}
close(N);
open(M, $ARGV[1])||die("Need a map file \$ARGV[1]\n");
while(<M>){
    chomp;
    my @a=split(/\t/);
    next unless defined($prots{$a[0]});
    $prots{$a[1]}++;
}
close(M);
open(G, $ARGV[2])||die("Need a map file \$ARGV[1]\n");
while(<G>){
    next if /^!/;
    my @a=split(/\t/);
    next unless defined($prots{$a[1]});
    print;			
}
