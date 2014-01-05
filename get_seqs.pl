#!/usr/bin/perl -w

# This script takes as input a list of FASTA ID lines
# and read a sequence file containing these sequences from STDIN. It will print out those seqs
# of the seq file whose ids were in the id file 

# USAGE  zcat/cat FASTA input file | get_seqs.pl name.file > outfile

use strict;

my %seqs = ();
my %gis = ();
my %names = ();

my @names;

my $match = 0;
my $name;

my @seq;


&usage unless $ARGV[0];

open(NAMES,"$ARGV[0]")|| die("can't open $ARGV[0] : $!\n");  ### file with FASTA ID lines, one per line 

while(<NAMES>)
{
 ### get list of spc names
    /^(.*)\n/;
    push @names, $1;
}

close(NAMES);

#print STDERR "names : $names[2]\n"; die;
while(<STDIN>)
{
    
    if (/^>/)
    {
#print STDERR "m : $match\n";
	### if we have reached the end of the last desired sequence ==> match == 1
	### put sequence list into hash and set @seq to 0

	if ($match == 1)
	{
	    $seqs{$name} = [@seq]; ### $seqs{gi} = @sequence
	    @seq = ();
	}

	$match = 0;
	
	foreach my $n (@names)   ### foreach species name
	{ 
	    if ($_ =~ /$n/)      ### if FASTA ID line matches name
	    {
	#	print STDOUT "name : $n\n";
	#	print STDOUT "\$_ : $_\n";
		print STDERR "*";
		$match = 1;      ### set match to 1
		$name = $_;      ### and name to ID line
	    }
	    else {next; }
	}
    }
    else 
    {
	if ($match != 0)   ### if this sequence is desired
	{
	    push @seq, $_;
	    next;
	    
	}
	next if $match == 0;  ### next if it is  not
    }
}


my @sequences = keys(%seqs);
#print STDERR "seqs : @sequences\n"; die();
map { print STDOUT "$_"; print STDOUT @{$seqs{$_}};  } @sequences;  ### produce FASTA output


sub usage
{
    print STDERR "USAGE : zcat/cat FASTA input file | get_seqs.pl name.file > outfile\n";
    exit(0);
}
