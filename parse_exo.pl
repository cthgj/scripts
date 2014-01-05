#!/usr/bin/perl -w

use Getopt::Std;
use strict;


my %opts;
getopts('hq:',\%opts);

&usage() if $opts{h};
&usage() unless $ARGV[0];

my (%qlength, %tlength, %id);

my $quer_name = $opts{q} || undef;

my ($query, $target);
my %hits;
my $prevq = 0;
my (@hsp, @targets);
my%hsps;

while (<>)
{
    if (/Query:/) 
    { 
	/Query:\s(.*)/; 
	$query = $1; 
    }
    if (/Target:/)
    {
	/Target:\s(.*)/; 
	$target = $1;
	if ($query eq $prevq)
	{
	    push @targets, $target ; 
	}
	else
	{
	    $prevq = $query;
	    $hits{$query} = [@targets];
	    @targets = ();
	}
    }
    if( /^[^\s]/) 
    { 
	my $n = $query . "##" . $target if $query && $target;
	if( /vulgar/) 
	{
	    $hsps{$n} = [@hsp];
# print STDOUT "Query : $query\nTarget : $target\n@hsp\n";
	    @hsp = ();
	}
      	if (/^RES/) {/[-+]\s\d*?\s(\d*?)\s(\d*?)\s([\.\d]*)/; $qlength{$query} = $1; $tlength{$target} = $2; $id{$n} = $3; }
    }
    else
    {
	push @hsp, $_;
    }
}

if ($quer_name)
{
    my $q = $quer_name;

    foreach my $hit (@{$hits{$q}}) { my $k = $q . "##" . $hit; print STDOUT "aa :$q($qlength{$q}) : $hit($tlength{$hit}), $id{$k}% \n@{$hsps{$k}}\n" } 
		     
}
else{print "shit!\n";}

sub usage
{   
    open(HELP, "| more") ;
    print HELP <<EndOfHelp;

USAGE:   parse_exo.pl [q] exonerate.out 
    
parse_exo will take an exonerate out file as input and either :

	a) Return all HSPs, slightly trimming the file
        b) Return only the HSPs for the query name given by the -q flag

EndOfHelp
close(HELP);
    exit(1);

}
