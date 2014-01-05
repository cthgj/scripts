#!/usr/bin/env perl

my $target=$ARGV[0];
my $AnnotStats=$ARGV[1];
my $InterStats=$ARGV[2];

open(A,"$target")||die("Could not open file $target : $!\n");
while (<A>) {
    if(/id=.InterTab/){
	print $InterStats;
    }
    elsif (/id=.AnnotTab/) {
	print $AnnotStats;
    }
    else{print}
}
