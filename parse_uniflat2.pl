#!/usr/bin/perl
use strict;
use warnings;
use Data::Dumper;
use Getopt::Std;
use SWISS::Entry;
my %opts;
getopts('hIfFen:t:i:',\%opts) || do { print "Invalid option\n"; exit(1); };
&usage() if $opts{h};
&print_ID_types if $opts{I};
my $toID=$opts{t}||"AC";
my $fromID=$opts{f}||'ID';
my $flat=$ARGV[0]||die("Need a flat file\n");
my $names = $opts{n}||undef;
my $IDlist;
my %entries;
my %ids;
my $return_fasta=$opts{F}||undef;
my $looking_for;
my $return_entry=$opts{e}|| undef;


## get desired ids
if($names){
    map{
	$ids{$_}++;
	$looking_for++;
    }split(/\s+/,$names);
}
else{
    $IDlist=$opts{i}||die("Need a list of desired IDs. Either as a text file (-i) or as a space separated list (-n)\n");
## Read desired IDs. Expecting a text file, one ID per line
    open(ID,$IDlist)||die("Cannot open ID list $IDlist :$!\n");
    while(<ID>){
	chomp;
	$ids{$_}++;
	$looking_for++;
    }
}

# Change the line termination string so we read an entire entry at a time
local $/ = "\n//\n";

if($flat=~/\.gz$/){
    open(A,"zcat $flat |")|| die("cannot open flat file $flat : $!\n");
}
else{
    open(A,"$flat")|| die("cannot open flat file $flat : $!\n");
}

while (<A>) {
    exit(0) if $looking_for==0;
#    printf STDERR "($looking_for,$.)\r";
    # Read in all the entries and fill %entries 
    my $entry = SWISS::Entry->fromText($_);
    my $protID=${$entry->IDs->elements}[0];
    print "pp : $protID\n";
	
}
