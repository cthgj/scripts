#!/usr/bin/perl -w

use strict;
my %ego;
my %tog;
my @array;
open(SEQ,"/home/ug/cchapple/research/seq/EGO/ego3_020103_seq_clusters_nopine.tbl") || die "0 cannot open seq: $!";

my $prevtog =  'TOG248582'; ### The first TOG# of the file.
   
while (<SEQ>)
{
    /^(.*?)-(.*?\s.*?)$/;
    if ($1 eq $prevtog)
    {
	push @array, "$1 "   . $2 . "\n";  ### array contains TOG# . species . sequence
    }
    else
    {
	$tog{$prevtog} = "@array";  ### $tog{TOG#} = all names and sequences of that tog.
	$prevtog = $1;
	@array= ();
    }
}
close (SEQ);

open(FILE, "/home/ug/cchapple/research/seq/EGO/known_sel_clusters")|| die "1 cannot open clusters: $!";

my @sels = <FILE>;

open(SEL, ">/home/ug/cchapple/research/seq/EGO/sel_clusters.tbl");

foreach my $sel (@sels)
{
    chomp $sel;
    print SEL "$tog{$sel}\n";
}
close(SEL);
exit();

	
