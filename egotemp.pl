#!/usr/bin/perl -w

use strict;

my %ego = ();


&ego();

open(FILE, "$ARGV[0]");
while(<FILE>){
    /^(.*?)$/;
    print STDOUT "$1 $ego{$1}\n";
    
}


sub ego   ### Get full gene sequence
{
    {
	open(SEQ,"/home/ug/cchapple/research/seq/4ego/ego4_060903.seq") || die "0 cannot open seq: $!";
	while (<SEQ>)
	{
	    /^(.*?)\s(.*?)$/;
	    $ego{$1} = $2;     ### $ego{seq name} = sequence
	}
	close (SEQ);
	print STDERR "*** get full seq done\n";
    }
}
