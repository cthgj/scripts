#!/usr/bin/perl -w


# This script will parse blast .out files and return those queries with NO HSPs 

use strict;

my $query;
my %queries = ();


while(<>)
{
    if (/^Query=/)
    {
	/Query=\s(.*)/;
	$query = $1;
    }
    if (/No\shits\sfound/)
    {
	print STDOUT "$query\n";
    }
    else { next;}
}


#my @keys = keys(%queries);

#map {print STDOUT "$_\n"} @keys;
