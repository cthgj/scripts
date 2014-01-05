#!/usr/bin/perl -w
# Script to convert embl ID files to FASTA
my $a=0; 

while(<>){
    if(/^ID/){
	/^ID\s+(.*?)\s/;
	print ">$1 ";
	$a=0;
    } 
    if(/^DE/ && $a==0){
	/^DE\s(.*)/; 
	$a=1; 
	print "$1\n";
    } 
    if(/^\s/){
	s/\s+//g;
	s/\d+/\n/g; print;
    }
}
