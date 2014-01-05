#!/usr/bin/perl -w

use strict;
use Getopt::Std;
$ARGV[0] =~ /^(.*)\.exo/;
my %opts;
getopts('bf:',\%opts);
my $n;
my $fastafile = $opts{f} || die("need a fasta file (-f)\n");
my $name;
my $best = $opts{b} || 0;
my $query;
my $k=0;
while(<>)
{
  
    if(/Query:\s*(.*?)\s/)
    {
	$query=$1;
	$query =~ s/\|/_/g;
    }
    if (/exonerate:protein2genome/ && /\sexon\s/)
    {
	
	/^((.+?)_.+?)\s/ || die("$1 : $2 $_\n");
#	die("$1 : $2 ---  $_\n");
	$n=$2;
	$name = $1;
	print STDERR "name : $name\n";
	
	unless (-r "$name.faa")
	{
	    print STDERR "retrieveseqs.pl -vfn $fastafile $name >$name.faa\n";
	    system("retrieveseqs.pl -vfn $fastafile $name >$name.faa");
	}
	open(F,">$name.gff");
	/^(.*\.\s\w+)/|| die();
	my $gffline = $1;
	if($gffline =~ /.\t-\t/)
	{
	    $gffline =~ /\t(\d+)\t(\d+)/ || die();
	    my $a = $1;
	    my $b = $2;
	    my $c=$a+1;
	    my $d=$b+1;
	    $gffline =~ s/$a/$d/;
	    $gffline =~ s/$b/$c/;
	    print F "$gffline\t.\n"; 
	    print STDOUT "SSgff $name.faa $name.gff >> $query.dna\n";
	    system("SSgff $name.faa $name.gff >> $query.dna");
	    
	}
	elsif($gffline =~ /.\t\+\t/)
	{
	    print F "$gffline\t.\n"; 
	    print STDOUT "SSgff $name.faa $name.gff >> $query.dna\n";
	    system("SSgff $name.faa $name.gff >> $query.dna");
	}
	else{die("gff parse problem : $!\n");}
    


	
    }
}
