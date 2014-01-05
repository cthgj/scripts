#!/usr/bin/perl -w
## This script will take a .map file, eg:
##
## 1433B_HUMAN	P33992 MCM5_HUMAN
## 1433E_HUMAN	P62258 1433E_HUMAN
##
## and a fasta file of this format:
##
## >sp|P33992|MCM5_HUMAN DNA replication licensing factor MCM5 OS=Homo sapiens GN=MCM5 PE=1 SV=5
##
## It will add to the appropriate FASTA header, any IDs in the map file to 
## that are not there and should be. Eg:
##
## >sp|P33992|MCM5_HUMAN|1433B_HUMAN DNA replication licensing factor MCM5 OS=Homo sapiens GN=MCM5 PE=1 SV=5
##
## Useful to map names from obsolete network iles to sequences

use strict;
&usage() unless $ARGV[1];
my $mapfile=$ARGV[0];
my $fasta=$ARGV[1];

my %map;

open(A,"$mapfile")||die("Cannot open $mapfile: $!\n");
while(<A>){
chomp;
my @a=split(/\s+/,$_);
push @{$map{$a[1]}{NEW}},$a[2];
push @{$map{$a[1]}{OLD}},$a[0];
}
close(A);

open(F,"$fasta")||die("Cannot open $fasta: $!\n");
while(<F>){
    if(/>/){
	chomp;
	/(>.+?)\s+(.+)/;
	my ($name, $description)=($1,$2);
	$name=~s/(>..\|)//;
	my $header=$1;
	my @a=split(/\|/,$name);
	my $netname;
	# if(defined($map{$a[1]}{ACC})){ $netname=$map{$a[1]}{ID};}
	# elsif(defined($map{$a[2]}{ACC})){ $netname=$map{$a[2]}{ID};}
	# else{die("PROBLEM: $_\n")}
	
	if(defined($map{$a[0]})){
	    foreach my $old (@{$map{$a[0]}{OLD}}){
		unless(/\|$old/){
		   $name="$old|$name"
		}
	    }
	    print "$header$name\t$description\n";
	}
	else{print "$_\n"}
    }
    else{print;}
}

sub usage{
    print <<EndOfHelp;
Usage: add_misssing_names.pl map fasta

This script will take a .map file, eg:

1433B_HUMAN	P33992 MCM5_HUMAN
1433E_HUMAN	P62258 1433E_HUMAN

and a fasta file of this format:

>sp|P33992|MCM5_HUMAN DNA replication licensing factor MCM5 OS=Homo sapiens GN=MCM5 PE=1 SV=5

It will add to the appropriate FASTA header, any IDs in the map file to 
that are not there and should be. Eg:

>sp|P33992|MCM5_HUMAN|1433B_HUMAN DNA replication licensing factor MCM5 OS=Homo sapiens GN=MCM5 PE=1 SV=5

Useful to map names from obsolete network files to sequences
EndOfHelp
exit();
}
