#! /usr/bin/perl 

## This script will take a list of names and a FASTA file and return those names
## that are not found as part of the ID of any seq in the fasta file

## The FASTA file is expected to have various IDs separated by "|"
## They will be split on the "|" and that list will be searched.

## eg.:
## >sp|A0JLT2|MED19_HUMAN|Q8IZD1_HUMAN	Mediator of RNA polymerase II transcription subunit 19 OS=Homo sapiens GN=MED19 PE=1 SV=2

## The ids extracted will be (note the absence of "sp":
## A0JLT,MED19_HUMAN,Q8IZD1_HUMAN


## USAGE find_missing_names_in_fasta_file.pl names fasta

my %k; 
open(A,"$ARGV[0]");
while(<A>){
    chomp; 
	$k{$_}++
}
open(B,"$ARGV[1]");
while(<B>){
    next unless />/; 
    @a=split(/\s+/,$_); 
    $a[0]=~s/>..\|//; 
    @b=split(/\|/,$a[0]); 
    for($i=0; $i<scalar(@b); $i++){
	if (defined($k{$b[$i]})){
	    $k{$b[$i]}++; 
	}
    }
}
map{print "$_\n" unless $k{$_}>1}keys(%k);
