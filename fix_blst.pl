#!/usr/local/bin/perl -w

## This script will merge different blast outfiles into one. For example
## if you have sent the same query file to more than one process in an effort 
## to get it done faster, or if the damn thing crashed or whetever, and you 
## have ended up with duplicate entries for the same query in the different 
## outfliles, just cat the outfiles through this script and you get a normal 
## blast outfile with just one of the entries for each query 

use strict;
my %hsps;
my $ll=0;
my $query;
my @lines;
my $blast;

while(<>)
{
    $. == 1 && do { /^(.*?)\s/; $blast = $_;};

    if(/Query=/)
    {
	$hsps{$query} = [@lines] unless $ll == 0;
	/Query=\s(.*)/;
	$ll = 1;
	$ll = 0 if $1 =~ /^\s*$/; #/Query=\s+$/;
	next if $ll == 0;
	$query = $1; 

	@lines = ();
    }
    push @lines, $_;

}
 
$hsps{$query} = [@lines] unless $ll == 0;

my @keys = keys(%hsps);

print "$blast\n\n";
foreach my $key (@keys)
{
    map {print "$_"}@{$hsps{$key}}
}

