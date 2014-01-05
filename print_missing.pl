#!/usr/bin/perl -w
use strict;
use Getopt::Std;
my %opts;
getopts('hf', \%opts);
my $fast=$opts{f}||undef;
&usage() unless $ARGV[1];
&usage() if $opts{h};
my %k;

open(F1,"$ARGV[0]")||die("Cannot open $ARGV[0]: $!\n");
while(<F1>){
    chomp;
    $k{$_}++;
}
close(F1);

open(F2,"$ARGV[1]")||die("Cannot open $ARGV[1]: $!\n");
while(<F2>){
    if($fast){
	my @a=split(/\s+/,$_);
	$k{$a[0]}++ if defined($k{$a[0]});
    }
    else{
	foreach my $key(keys(%k)){
	    if (/^$key\W/ || /\W$key$/ ||  /\W$key\W/ ){
		$k{$key}++ ;
	    }
	}
    }
}
close(F2);
map{print "$_\n" unless $k{$_}>1}keys(%k);


sub usage{
    print STDERR <<EndOfHelp;

  USAGE: print_missing.pl FILE1 FILE2

This script will print all the lines of FILE1 that are not
found in FILE2. FILE one must be one search pattern per line,
the search pattern need only be contained within one of the 
lines of FILE2.

The -f flag assumes that the search pattern should be the first
characters of each line of FILE2 until the first non-word char.

EndOfHelp
    exit(0);
}
