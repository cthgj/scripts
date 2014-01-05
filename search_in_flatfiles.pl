#!/usr/bin/perl -w 

use strict;
my $flatfile=$ARGV[0]||die("no ARGV[0] : $ARGV[0]\n");
my $pattern=$ARGV[1]||die("no ARGV[1] : $ARGV[1]\n");
my $id=$ARGV[2]||die("no ARGV[2] : $ARGV[2]\n");
	

open(A,"$flatfile");
my $c=0;
my @syns;
while(<A>){
    
    if(/$id/){
	$c=1;
	next;
    }
    next unless $c>0;
    if(/^DR\s+(.+)/){
	my $a=$1;
	my @b=split(/\;/,$a);
	push @syns,@b;
	if(/($pattern.+?)[\.\;\,]?\s/){
	    print "$id\t$1\n" unless $c>1;
	    $c++;
	}
    }
    if ((/^ID/) && ($c>1)){
#	print "@syns\n";
	exit(0)
    }
}
print "$id\tmissing\n";

