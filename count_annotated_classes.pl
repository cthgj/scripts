#!/usr/bin/perl -w

use Getopt::Std;
use strict;

for my $file (@ARGV) {
    my ($classes, $ca, $cm)=(0,0,0);
    open(A, $file);
    while (<A>) {
	next if /#/;
	next if /^\s*$/;
	$classes++ if /^\[CLASS/ ;
	if (/^CA/){
	    $ca++ unless /_unknown/;
	}
	if (/^CM/){
	    $cm++ unless /_unknown/;
	}
    }    
    close(A);
    my ($perA, $perM)=(0,0);
    if ($classes>0) {
	$perA=($ca * 100)/$classes;
	$perM=($cm * 100)/$classes;
    }
    
    print "$file:\n\t$classes Classes\n\t  $ca (";
    printf("%.2f",$perA);
    print "%) Classes have 1ary annot.\n\t  $cm (";
    printf("%.2f",$perM);
    print "%) Classes have 2ary annot.\n";
}

