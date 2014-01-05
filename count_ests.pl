#!/usr/bin/perl -w

# This script runs on the output of ~scaste/bin/Add_Cys_pos2queryID.awk. It prints out
# those queries which have alignments with at least a user specified (-m) minimum of different subjects.
# If no min specified, the default (>=5) is used. 
# -s for silent mode...
# The output format is  Query_name<TAB>space separated list of subject_names



use strict;
use Getopt::Std;
my %peptides = (); ### query
my %ests     = (); ### sbjct

my $query;
my $sbjct;


my %opts;
getopts('sm:',\%opts);

my $est_number = $opts{m} || 'nnn';
my $s = $opts{s} || 0;
#print STDERR "$est_number\n"; #die();
if ($est_number eq 'nnn') 
{
    print STDERR "No number specified, using >4...\n" unless $s; $est_number = 5
}
else
{
    print STDERR "Selecting queries supported by >$est_number ESTs\n" unless $s;;
}

while(<>)
{
    next if /^\s+$/;  ### skip blank lines
    
    if (/^Query=/)
    {
	/Query=\s(.*?)\s/;
	$query = $1;
	next;
    }
    
    if (/^>/)
    {
	#s/\s(.)/_$1/g;
	/>(.*)$/;
	$sbjct = $1;
	
	push @{$ests{$query}} , $sbjct
    }

}
my @keys = keys(%ests);

    map {
	my $a = @{$ests{$_}};
	if ($a >= $est_number) {
	    print STDOUT "$_\t@{$ests{$_}}\n" 
	    }
    } @keys;
