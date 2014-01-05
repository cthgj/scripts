#!/usr/bin/perl -w

### This script will take a list of IDs ($ARGV[0]) and look through the gff output of SECISearch ($ARGV[1]) 
### printing a gff file of IDs present in both files


my %secis = ();
my %shared = ();
open(FILE, "$ARGV[0]");
while(<FILE>)
{
    /^(.*?)$/;
    $secis{$1}++
}
close(FILE);

open(FIL,"$ARGV[1]");
while(<FIL>)
{
    /^(.*?)\s/;
    $shared{$1} = $_;
}
close(FILE);

foreach my $a(keys(%secis))
{
    print STDOUT "$shared{$a}" if $shared{$a};
}
