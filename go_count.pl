#!/usr/bin/perl -w
## This script will count all the occurrences (implicit and explicit) of the given GO 
## term in a gene_association file.
##
## Usage: go_count.pl -o <ontology GAF file> -g <geneology file> GO:term
use strict;
use Getopt::Std;

my (%opts,%ancestors,%want,%seen);
getopts('g:o:',\%opts);

my $subonto="P";
my $gaf_file=$opts{o}||"gene_association.goa_human";
my $geneology_file=$opts{g}||"biological_process.geneology";
my $prot_name="Bob";
my $go=$ARGV[0]||die("Need a go term");
my $count=0;
my @children;
open(A,"$geneology_file")||die("cannot open $geneology_file: $!\n");
while(<A>){
    next unless /$go/;
    chomp;
    my @terms=split(/\t/);
    $want{$terms[0]}++;
}
close(A);
open(A,"$gaf_file")||die("cannot open $gaf_file: $!\n");
while(<A>){
    next if /^!/; 
    chomp;
    my @tmparray=split(/\t/);
    next unless defined($want{$tmparray[4]});
    next unless $tmparray[3]=~/^$/;  ## skip the NOT annotations
    my $nn=$tmparray[2]; ## $nn == protein name
    $count++ unless defined($seen{$nn});
    $seen{$nn}++;
}
close(A);
print "$count occurrences of $go\n";
