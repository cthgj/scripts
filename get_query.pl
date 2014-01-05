#!/usr/bin/perl -w

# This script will take one or more query names and return their entries form a blast outfile
# get_query.pl outfile -n name

use strict;
use Getopt::Std;
my %opts;
&getopts('n:f:',\%opts);
my $file = $opts{f} || undef;
my @queries;
my $c=0;
#my $blastfile  = $ARGV[0];
my $wanted_query = $opts{n} || undef;

my $qq = 0;  # counter
my %q;
my $name;
my $got=0;

my $num;


if($file)
{
    open(F,"$file");
    while(<F>)
    {
	chomp;
	$q{$_}++;
	$num++;

    }
    @queries = keys(%q);
}
else
{
    $queries[0] = $wanted_query;
    map{$q{$_}++; $num++} @queries;
     
}

#map{ print "$c : $_\n"; $c++;}@queries;

my $kk=0;

while(<>)
{
   
	if(/^\s?.?BLAST/)
	{
	    print STDOUT  "$_\n\n" if $kk==0; 
	    $kk=1;
	    
	    if($qq == 1)
	    {
		$got++;
		print STDOUT  "$_\n\n";
	    }
	    $qq = 0;
	}
	elsif(/Query=/)
	{
	    /Query=\s?([^\s]*)/ || die("cannot match query : $_\n");
	    $name = $1;
	    if (defined($q{$name})){$qq=1}
	    else{$qq=0}
	}
	else
	{
	    exit (0) if $got == $num;
	    next unless $qq==1;

	}
	printf(STDERR "$.\r");
	exit (0) if $got == $num;
	print STDOUT if $qq==1;
	
    }

