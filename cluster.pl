#!/usr/bin/perl -w

#   This script will take the .secis output of SECISearch and create separate tbl files for each 
#   cluster which returned >2 SECISes. The file cluster_morethan2_PATNAME is created by the 
#   parse_egosecis.sh script. 
#       
#           Correct usage is "cluster.pl [SECISout.secis] [cluster_morethan2_PATNAME]
#

use strict;

my @clusname = ();

my $secisout = $ARGV[0];   # get secisout filename
chomp $secisout;

my $clusfile = $ARGV[1];   # get cluster name file
chomp $clusfile;

$clusfile =~ /cluster_morethan2_(.*?)$/;
my $patname = $1;



open (CLUSTERS, "$clusfile")|| die "cannot open $clusfile: $!";

my @clusters= <CLUSTERS>;  # get cluster names
chomp @clusters;

foreach my $a (@clusters)
{
    $a =~ /\d\s?(.*?)\s?$/g;
    my $lala = $1;
    push @clusname, $lala ;
}



system("FastaToTbl $secisout > tmptbl");


open (FILE, "tmptbl")|| die "cannot open tmptbl: $!";

while (<FILE>)
{
    foreach my $cluster (@clusname)
    {
	if (/$cluster/)
	{
	    open (CLUS, ">>$cluster\_$patname\.tbl")|| die "cannot open $cluster\.tbl: $!";
	    print CLUS;
	    close (CLUS);
	}
    }
}
close(FILE);
unlink ("tmptbl");
