#!/usr/bin/perl -w

#
#  This script is trying to produce two things. A file, clusters.tmp, which is each cluster
#  and the sequences in it and the STDOUT which is each sequence and the clusters it belongs to.
#             USAGE :   manyclusters.pl egoX_xxxxxx.orth > outfile
#

use strict;

my $c = 0;
my %T = ();
my @array;
my @names;
my $prevcluster = '1';

unless ($ARGV[0])
{
    print STDERR "\nUSAGE :   manyclusters.pl egoX_xxxxxx.orth > outfile\n\n";
exit();
}

open (FILE, "$ARGV[0]")|| die "cannot open $ARGV[0] $_";

while(<FILE>)
{
    /(.*?)\s/;
    push @names, $1;  ### @names will have TOGs and names
}
close(FILE);

open(CLUS, ">clusters");

foreach my $name (@names)
{

   
    if ($name =~ /^>(TOG.*?)$/o)
    {	

	print CLUS "$prevcluster\t@array\n" unless $prevcluster eq '1';

	my $cluster = $1;
	$prevcluster = $cluster;
	@array=();
	$c++;
	print STDERR "_" unless $c % 1000;  ;
    } 
    else
    { 
	chomp $name;
	push @array, $name ;
    }
}

print CLUS "$prevcluster\t@array\n";
print STDERR "\n";

close(CLUS);
open(CLUS,"clusters");
###########################
$c = 0;
while (<CLUS>)
{
    my  @F = split /\s+/og, $_;

    for (my $n=1; $n < scalar(@F); $n++) 
    {
	defined($T{$F[$n]}) || (@{ $T{$F[$n]}} = ()); 
	push @{ $T{$F[$n]} }, $F[0]; 
	$c++;
	print STDERR "." unless $c % 1000;  
    };
 
}
$c = 0;
print STDERR "\n";

	foreach my $k (sort keys %T) 
	{
	    $c++;
	    print STDERR "+" unless $c % 1000; ;
	    print STDOUT "$k: @{ $T{$k} }\n"; 
	}; 
print STDERR "\n";



