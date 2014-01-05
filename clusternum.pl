#! /usr/bin/perl -w

# This script takes the list of erpin hits and the clusters they belong to
# and returns the copy number of each cluster.


use strict;
my %clusters;



open (FILE, "$ARGV[0]");
while (<FILE>)
{
    /^.*?\s(.*)$/og;
    my @a = split(/\s/, $1);
    map{ $clusters{$_}++ }@a;
}
close(FILE);

my @b = keys(%clusters);
map{ print STDOUT "$clusters{$_} $_\n"} @b;
