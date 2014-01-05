#! /usr/bin/perl -w

# This script will take as input a file with erpin hitcounts for each cluster and 3 files with SECISearch 
# hitcounts for each cluster (one file per pattern)

#   FORMAT example : tmp.pl ego4_060903.fa.erpin.clustercount ../std/cluster_morethan2_std ../non_std/cluster_morethan2_non_std ../twil/cluster_morethan2_twil


use strict;
my %erpin = ();
my %secis = ();

open(FILE,"$ARGV[0]");
while(<FILE>)
{
    /^(\d)?\s(.*?)$/og;
    $erpin{$2} = $1;
    $secis{$2} = 0;
}
close(FILE);
open(FILE, "$ARGV[1]");
while(<FILE>)
{
    /^(\d)?\t(.*?)$/og;
    $secis{$2} = $1;
}
close(FILE);

open(FILE, "$ARGV[2]");
while(<FILE>)
{
    /^(\d)?\t(.*?)$/og;
    $secis{$2} = $1;
}
close(FILE);

open(FILE, "$ARGV[3]");
while(<FILE>)
{
    /^(\d)?\t(.*?)$/og;
    $secis{$2} = $1;
}
close(FILE);


my @a = keys(%erpin);
my @b;


foreach my $key (sort { $erpin{$b} <=> $erpin{$a} } keys %erpin)
{
  push @b, "$key \terpin : $erpin{$key}  SECISearch : $secis{$key}\n";
}
print STDOUT @b;
