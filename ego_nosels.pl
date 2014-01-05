#!/usr/bin/perl -w

# This script produces a tbl file of all EGO sequences which do not belong
# to a cluster which produced a significant hit with a known human selenoprotein


use strict;




open(SELS, "/home/ug/cchapple/research/seq/EGO/known_sel_clusters") || die "cannot open file $_";

my @sels = <SELS>;
close(SELS);

open (FILE, "/home/ug/cchapple/research/seq/EGO/ego3_020103_seq_clusters_nopine.tbl") || die "cannot open file $_";

while (<FILE>)
{
    foreach my $sel(@sels)
    {
	unless (/^$sel/)
	{
	    print;
	}
    }
    
}

close(FILE);
