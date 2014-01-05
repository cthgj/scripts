#!/usr/bin/perl -w

#
# This script takes the EGO database in pseudo tbl format (one line per sequence)
# writes each cluster into a tmp file, runs tblastx of the file against itself and then
# runs SECISearch on it.

use strict;
use Getopt::Std;


my %opts;
getopts('sbc',\%opts) || do { print "Invalid option"; exit(1); };

my $cluster;
my $id;
my $seq;
my $prevcluster;
my $counter = 0;

my $SECISearch = '/home/ug/cchapple/scripts/SECISearch.pl';
my $secisout   = 'SECISearch\.out';
my $secislog   = 'SECISearch\.log';

my $s = $opts{s} || 0;
my $b = $opts{b} || 0;
my $c = $opts{c} || 0;




open (FILE, "/home/ug/cchapple/research/selenoproteins/EGO/blast/non_sel/non_sel.tbl")|| die "cannot open file $_";

my $fline = <FILE>;


&get_cluster($fline);
$prevcluster = $cluster;

while (<FILE>)
{
   	 
    if (&get_cluster($_) eq $prevcluster)
    {

	
#	print STDERR "\$prevcluster is $prevcluster\n"; 
#	print STDERR "\$cluster is  $cluster\n";


	open (TMP , ">>$cluster\.tmp")|| die "cannot open $cluster\.tmp $!";

	print TMP "$cluster\-$id $seq\n";
	$prevcluster = $cluster;
	close (TMP) || die "cannot close $cluster\.tmp $!";
	print STDERR "just closed $cluster\.tmp\n"; 
    }
    elsif($c)
    {

	system("cp $prevcluster\.tmp $prevcluster\.tbl");
	unlink ("$prevcluster\.tmp");
	$prevcluster = $cluster;


    }
       
    else
    {
	&blastit();

        &secisearch();
#       print STDERR "else \$cluster is now $cluster\n";
#	print STDERR "\$prevcluster is now $prevcluster\n"; 
	unlink ("$prevcluster\.tmp");
	unlink ("$prevcluster\.ftmp");
	$prevcluster = $cluster;

#	print STDERR "\$prevcluster is finally $prevcluster\n"; 
	
    }
    
}


###################################################


sub get_cluster
{
    $_[0] =~ m/^(.*?)-(.*?)\s(.*?)$/;
    $cluster = $1;
    $id = $2;
    $seq = $3;
    return $cluster;
}

sub blastit
{
    unless ($b)
    {
	
	
	
	system("TblToFasta $prevcluster\.tmp > $prevcluster\.ftmp"); 
	system("wu-formatdb -p F -t $prevcluster\.ftmp -i $prevcluster\.ftmp");
	    print STDERR "********************\nwu-tblastx $prevcluster\.ftmp $prevcluster\.ftmp > $prevcluster.out\n**********************";
	    
	    
	    system("wu-tblastx $prevcluster\.ftmp $prevcluster\.ftmp > $prevcluster.out");#||  die "cannot run tblastx: $!";
	unlink("$prevcluster\.ftmp\.xnt");
	unlink("$prevcluster\.ftmp\.xns"); 
	unlink("$prevcluster\.ftmp\.xnd");
    }
}

sub secisearch
{
    unless ($s)
    {
	system("$SECISearch -sdlI -f $prevcluster\.tmp >>$secisout")|| die "cannot run $SECISearch: $!";
    }
}
