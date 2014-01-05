#!/usr/bin/perl -w

# This script will return the longest HSP... dunno really, read the sourcecode
# was written specifically for tetraodon selnprtn search. Go through the code before use...


use strict;
use Getopt::Std;
my %tetra = ();
my $tk = 0.1;
my $sk = 0.1;
my $lasttetra;
my %sel = ();
my %seq = ();
my %hum = ();
my %tet = ();
my %tetseq = ();
my %opts;
getopts('ts',\%opts) || do { print "Invalid option,-s for sels and -t for tetra"; exit(1); };

my $tet  = $opts{t} || 0;
my $sp   = $opts{s} || 0;


open(FILE,"$ARGV[0]"); ### Open parsed(parseblast.pl -P) blast output
while(<FILE>)
{
    next if /^\#/;  ## skip #s and tetra sequences
    if (/^F/)  ### if tetra sequence
    {
	/^(.*?)_.*?\s(.*?)$/;
	$lasttetra = $2;
	my $la = length($2);
	if ($tetra{$1})
	{
	    $tk = length($tetra{$1});
	    if ($la > $tk)
	    {
		$tetra{$1} = $2; 
	    }
	    else {next;}
	}
	else
	{
	    
	    $tetra{$1} = $2; 
	    next;
	}
    }
    else
    {
	/^(.*?)\s+(.*?)$/;
	$sel{$lasttetra} = $1;
	$seq{$lasttetra} = $2;
	
	if ($sp)
	{
	    my $sa = length($2);
	    if ($hum{$1})
	    {
		$sk = length($hum{$1});
		if ($sa > $sk)
		{
		    $hum{$1} = $2; 
		}
		else {next;}
	    }
	    else
	    {
		$hum{$1} = $2; 
		next;
	    }
	}
    }
}

close(FILE);

if($tet)
{
    my @keys = keys(%tetra);
    
    foreach my $key (@keys)
    {
	my $la = length($key);
	my $lo = length($sel{$tetra{$key}});
	
	my $le = $la - $lo;
	
	print STDOUT "$key\t$tetra{$key}\n$sel{$tetra{$key}}" ;
	print STDOUT " " x $le;
	print STDOUT "\t$seq{$tetra{$key}}\n";
    }
}
elsif ($sp)
{
    my @keys = keys(%hum);
    
    foreach my $key (@keys)
    {
	my $la = length($key);
	my $lo = length($tet{$hum{$key}});
	
	my $le = $lo - $la;
	
	print STDOUT "$key";
	print STDOUT " " x $le;
	print STDOUT "\t$hum{$key}\n$tet{$hum{$key}}"; 
	print STDOUT "\t$tetseq{$hum{$key}}\n***\n"; 
    }
}

