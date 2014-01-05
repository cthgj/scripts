#!/usr/bin/perl -w

# simple parseblast script. Will return query names and the sequence they hit against. 
# If -n switch then returns queries with NO hits.
# If -m switch then returns one file per query with query and subject names
# If -i switch then returns only IDENTICAL or near identical (100% coverage >98%id) hits


### NOT very good, remember to make sure it works before using

use Getopt::Std;
use strict;


my %opts;
getopts('vnmi',\%opts);

my %nohits = ();
my %queries = ();
my @sbjcts;
my $query = undef;
my ($q,$subj,$sb,$id,$qlength,$slength,$coverage);
my $n = $opts{n} || undef;
my $m = $opts{m} || undef;
my $i = $opts{i} || undef;
my $verbose = $opts{v} || undef;

while(<>)
{
#    next unless (/Query=/ || />/ || /Identities/ || /Length/ || /letters\)/);#/\s\d{1,3}\s\s\d/);  #only parse relevant lines, Query and list of subjects, and stats
    if (/Query=/)
    { 
	print STDERR "-" if $verbose;
	/Query=\s(.*?)\s/;
	$q = $1;                   # Query name = $q
#	print STDOUT "q : $1\n";
	$nohits{$q}++;      # provisionaly, put $nohits{$q}++ until we see if it has a hit
	if ($query)         # if this is not the first query
	{
	    $queries{$query} = [@sbjcts];   # $queries{$query} = list of subjects it hit against
	    if ($m)
	    {
		unless ($n)
		{
		    my $new_name = $query;
		    $new_name =~ s/\|/_/g;
		    my $query1 = $query;
		    # special case where query name ends in |frame number (which number we leave in the filename)
		    $query1 =~ s/\|\d$//;
		    open(TMP,">$new_name");
		    print TMP "$query1\n";
		    map {print TMP "$_\n"} @{$queries{$query}}; 
		    close(TMP);
		    @sbjcts = ();     ## Set vars back to 0 to go through list again
		    $query = undef;   ##
		}
	    }
	    else
	    {
		print STDOUT "\n$query  :\n@{$queries{$query}}\n" unless $n;
		@sbjcts = ();     ## Set vars back to 0 to go through list again
		$query = undef;   ##
	    }
	}
    }
    elsif(/(\d+)\s+letters/){
	$qlength=$1;
    }
    elsif(/Length/){
	$slength=$1;
    }
    elsif(/Identities\s*=\s*(\d+)\/(\d+)\s*\((\d+)%/){
	$id=$3;
	$coverage=$1;
	$query = $q;
	if($i){
	    my $a=$qlength-$coverage;
	    print "aa : $a\n"; die();
	    $id==100 && push @sbjcts, $1 unless grep {$_ eq $1} @sbjcts;   # push subject name => @sbjcts unless already there
	}
	else{
	    push @sbjcts, $1 unless grep {$_ eq $1} @sbjcts;   # push subject name => @sbjcts unless already there
	}
	delete($nohits{$query});  # this one has  a hit ==> remove from %nohits 
    }
### i know, this is silly...
    elsif(/>/)
    {
	/^>(.*?\s)/; # get subject name
	$subj=$1;	
	print STDERR  "." if $verbose;
    }
}
if($n)
{
    my @keys = keys(%nohits);
    map {print STDOUT "$_\n";}@keys;
}


