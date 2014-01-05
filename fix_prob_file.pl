#!/usr/bin/perl -w

## This script will check that probs in a network prob file (eg HSHQ.gr.prob) are the
## same as thos in the data directory files, eg ~/research/testing/new/data/gostats/human/annotations/00023.prob

## Usage:  $0 <prob file> <prob dir>
use strict;
my $stats_dir=$ARGV[1]||die("Need a stats dir as argv[1], eg ~/research/testing/new/data/gostats/human/annotations/\n");
my (%f_probs,%prob);
open(A, "$ARGV[0]")||die("Need a network prob file argv[0]\n");
while(<A>){
    chomp;
    my @a=split(/\s+/);
    #my @g=split(/_/,$a[0]);
    $f_probs{$a[0]}=$a[1];
}
close(A);

foreach my $gopair (keys(%f_probs)){
    unless(defined($prob{$gopair})){
	$gopair=~/^GO:(.....)/||die("cannot match gopair : $gopair\n");
	my $file=$1 . ".prob";
	open(A,"$stats_dir/$file")|| die("cannot open $stats_dir/$file : $!\n$gopair\n");
	while(<A>){
	    chomp;
	    if(/$gopair/){
		my @tt=split(/\t/);
		$prob{$gopair}=$tt[6];	
		print STDOUT "$gopair\t$prob{$gopair}";
		if ($f_probs{$gopair} != $tt[6]){
		    print STDOUT "\tBAD\n" ;

		}
		else{print STDOUT "\tGOOD\n";}

	    }
	    else{next}
	}
	print STDERR "$gopair not found in $file\n" unless defined($prob{$gopair});
    }
    
}

