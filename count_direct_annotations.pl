#!/usr/bin/perl -w

use strict;
use Getopt::Std;
my (%names,%opts,%gos);
getopts('hs:',\%opts);
my $names=$ARGV[0]||die("Need a list of ids\n");
my $gaf=$ARGV[1]||die("Need a GAF file\n");
my $subonto=$opts{s} || "P";

open(A,"$names")||die("cannot open $names: $!\n");
while(<A>){
    chomp;
    $names{$_}++;
}
close(A);

open(A,"$gaf")||die("cannot open $gaf: $!\n");
while(<A>){
   	next if /^!/;
	chomp;
	my @tmparray=split(/\t/);
	next unless $tmparray[8] eq $subonto;
	next unless $tmparray[3]=~/^$/;  ## skip the NOT annotations
	next unless defined($names{$tmparray[1]});
	$gos{$tmparray[1]}{$tmparray[4]}++; ## $tmparray[4] == GO
}

map{
    print "$_\t",scalar(keys(%{$gos{$_}})),"\n";
}keys(%gos);
