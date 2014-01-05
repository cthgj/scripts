#!/usr/bin/perl -w
use strict;
use Getopt::Std;
use Switch;
my %opts;
getopts('',\%opts);

my @files=@ARGV;
    my %candidates;
foreach my $file (@ARGV){
    my %fcands;
	$file=~/.+\.([bimpPcd])\./;
	my $mode=$1;
    open(A,"$file")||die("Cannot open $file : $!\n");
    while(<A>){
	my @a=split(/\t/);

    ## Check what type of file this is
	switch($mode){
	    case /[bipP]/{
		$fcands{$a[0]}++;
	    }
	    case /[cdm]/{
		my @b=split(/\s+/,$a[0]);
		$fcands{$b[0]}++;
		@b=split(/\s+/,$a[2]);
		$fcands{$b[0]}++;
	    }
	    
	} ## end switch
    } ## end while(<A>)
	map{$candidates{$_}++}keys(%fcands);
} ## end foreach my $file

map{print "$_\t$candidates{$_}\n"}keys(%candidates);

