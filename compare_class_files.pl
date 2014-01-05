#!/usr/bin/perl -w

my $f1=$ARGV[0];
my $f2=$ARGV[1];

my %classes;
my $a=0;
foreach my $file ($f1, $f2){
    open(A,"$file");
    my $class;
    while(<A>){
	## Deal with different formats:
	$a++ if /Final\s*classes/i;
	next unless $a>0;
	if(/^\s*(\d+)\s*$/){
	    $class=$1;
	}
	else{
	    my @b=split(/\s+/);
	    map{$classes{$_}{$class}++}@b;
	}



    }
    $a=0;


}
