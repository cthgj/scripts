#!/usr/bin/perl -w
## this script should be run on the output of get_hsp_coords.pl.
## it will extract each target sequence from a multifasta file (-f), then
## use fastasubseq to extract the HSP region and then align the relevant protein 
## to it using genewise. The query is either extracted from a multifasta file (-p).
## or if there is only one query passed as an input fil (-P). 
## Genewise results in <blast query name>.genewise

use strict;
use Getopt::Std;
my %opts;
getopts('f:p:P:l:',\%opts);

my $fasta = $opts{f} || die("need a fasta file to extract the sequences from\n");
my $protein_file = $opts{p} || undef; ## multifasta of protein queries
my $protein = $opts{P} || undef;      ## FASTA file of protein query
my $length = $opts{l} || undef;       ## length of target sequence to extract for genewise run
my $q;
my @a; 
my %queries; 
my $l;

while(<>){
    @a=(); 
    /^(.*?)\s(.*?)\s(\d+)\s(\d+)\s([+-])/ || die("Match problem : $_");
    $q=$1; 
    push @a, $2, $3, $4, $5;
    push @{$queries{$q}}, [@a];
 }
my @k=keys(%queries);

if($protein_file)
{    
    foreach my $qu (@k){ 
	foreach my $s(@{$queries{$qu}}){
	    if($length){ $l = $length;}
	    else{$l=${$s}[2]-${$s}[1] +100;}
	    if(${$s}[3] eq "+") {
		print STDERR "aaa retrieveseqs.pl -vfn $fasta ${$s}[0] > a; fastasubseq  -s ${$s}[1] -l $l -f a > b; FastaToTbl $protein_file | grep \"$qu\" | TblToFasta > o; genewise -cdna o b >> \"$qu.genewise\"\n";
		system("retrieveseqs.pl -vfn $fasta ${$s}[0] > a; fastasubseq  -s ${$s}[1] -l $l -f a > b; FastaToTbl $protein_file | grep \"$qu\" | TblToFasta > o; genewise -cdna o b >> \"$qu.genewise\"");
	    } 
	    else{
		print STDERR "bbb retrieveseqs.pl -vfn $fasta ${$s}[0] > a; fastasubseq  -s ${$s}[1] -l $l -f a > b; FastaToTbl $protein_file | grep \"$qu\" | TblToFasta > o; genewise -cdna -trev o b >> \"$qu.genewise\"\n";
		system("retrieveseqs.pl -vfn $fasta ${$s}[0] > a; fastasubseq  -s ${$s}[1] -l $l -f a > b; FastaToTbl $protein_file | grep \"$qu\" | TblToFasta > o; genewise -cdna -trev o b >> \"$qu.genewise\"");
	    }
	}
    }
}
elsif($protein)
{
    foreach my $qu (@k){ 
	foreach my $s(@{$queries{$qu}}){
	    if($length){ $l = $length;}
	    else{$l=${$s}[2]-${$s}[1] +100;}
	    if(${$s}[3] eq "+") {
		print STDERR "aaa retrieveseqs.pl -vfn $fasta ${$s}[0] > a; fastasubseq  -s ${$s}[1] -l $l -f a > b; genewise -cdna $protein b >> \"$qu.genewise\"\n";
		system("retrieveseqs.pl -vfn $fasta ${$s}[0] > a; fastasubseq  -s ${$s}[1] -l $l -f a > b; genewise -cdna $protein b >> \"$qu.genewise\"");
	    } 
	    else{
		print STDERR "bbb retrieveseqs.pl -vfn $fasta ${$s}[0] > a; fastasubseq  -s ${$s}[1] -l $l -f a > b; genewise -cdna -trev $protein b >> \"$qu.genewise\"\n";
		system("retrieveseqs.pl -vfn $fasta ${$s}[0] > a; fastasubseq  -s ${$s}[1] -l $l -f a > b; genewise -cdna -trev $protein b >> \"$qu.genewise\"");
	    }
	}
    }
}
else{die("need a protein file\n");}
