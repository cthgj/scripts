#!/usr/local/bin/perl -w




use strict;
use Getopt::Std;

my %opts;
getopts('l:',\%opts) || do { print "Invalid option\n"; exit(1); };

unless ($ARGV[0])
{
    print STDERR "****\nNeed an exonerate file!!\n****\n";
    &usage();
}
my $filename = $ARGV[0];


## Set number of nts downstream to look for SECIS in
my $LENGTH = $opts{l} || 1000;


my ($query, $target);

open(FILE, "$ARGV[0]");

my @file = (<FILE>);
my $size = @file;



for(my $i=0; $i<$size; $i++)
{
    my $line = $file[$i];

   #  if ( $line =~ /^Query/) 
#     { 
# 	$line =~ /Query:\s(.*)/; 
# 	$query = $1; 
# 	next;
	
#     }
#     if ( $line =~ /^Target:/)
#     {
# 	$line =~ /Target:\s(.*)/; 
# 	$target = $1;
# 	next;
#     }
## Only the gff lines have "\t"
    if ( $line =~ /^[^\s].*\t.*gene/ )
    {
	$line =~ /^(.*?)\t.*?gene\t(\d*?)\t(\d+).*?([+-]).*?sequence\s(.*?)\s/ || die();
	$target = $1;
	my $start = $2;
	my $end = $3;
	my $strand = $4;
	$query = $5;

	my $new_end = $end + $LENGTH;
	if ($strand eq '-') 
	{
	    $start = $start - $LENGTH;
	    $end = $start;
	}
	

	open(GFF, ">$target.gff");
	print GFF "$target\t.\t.\t$end\t$new_end\t.\t$strand\t.\t.\t.\n";
	close(GFF);
	
	system("getfasta ../../Ptetraurelia_V2.fasta ../../index $target > temp.fa; SSgff.chr temp.fa $target.gff > check.fa" ); 
	my @patterns = ("s","n","t");
	
	foreach my $pat (@patterns)
	{
#	    print STDOUT "SECISearch -vp $pat check.fa > $query.secis\n"; 
	    system("SECISearch.pl -sp $pat check.fa >> $query.secis");
	}

	unlink "$target.gff";
#	print "t : $target q : $query s : $start e : $end n : $new_end\n";
#	print "$line";
    }
    
}
 
