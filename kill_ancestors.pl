#!/usr/bin/perl -w
use strict;

my $geneology_file=$ARGV[0];
my (%ancestors,%offspring);

&load_ancestors();

while(<>){
    my @a=split(/\t/);
    unless( (exists($offspring{$a[0]}{$a[1]}) ) ||
	    (exists($offspring{$a[1]}{$a[0]}) )	){
	print;
    }


}


sub load_ancestors{
    open (ANC,"$geneology_file")|| die("cannot open $geneology_file:$!\n");
    while(<ANC>){
	chomp;
	my @terms=split(/\t/);
	my $child=shift(@terms);
	$ancestors{$child}=[@terms];
	map{$offspring{$_}{$child}++}@terms;
    }
    close(ANC);
}
