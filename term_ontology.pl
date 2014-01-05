#!/usr/bin/perl -w

use strict;
my $termfile="/home/terdon/research/GO/GO.terms_alt_ids";
unless (-e $termfile){
    $termfile="/cobelix/chapple/backup/research/GO/GO.terms_alt_ids";
}
my $term=$ARGV[0];

open(A,$termfile)||die("Could not open $termfile:$!\n");
while(<A>){
    chomp;
    next unless /$term/;
    /\t([CPF])\t/;
    print "$1\n";
}
