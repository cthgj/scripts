#!/usr/bin/perl -w
## This script will parse one of benoit's class files search for each of the proteins
## given. Protein names can be given in a file (one per line) or with the optional
## -p flag, as a comma separated list in the command line. If no file or proteins are
## given, the gold standard set is used
use strict;
use Getopt::Std;


my %opts;
getopts('p:f:m:M:ah',\%opts);
&usage() if $opts{h};
unless ($ARGV[0])
{
    &usage("Need a class file!!");
}
my($name,$annot,$num,$class);
my (%singletons,%prot);
my @gold;
my $cc=0;

my $min=$opts{m}||0;
my $max=$opts{M}||100000;

if(defined($opts{p})){
    @gold=split(/,\s*/, $opts{p});
    $cc=1;
} 
if(defined($opts{f})){
    open(A,"$opts{f}")||die("cannot open $opts{f} for reading : $!\n");
    while(<A>){
	chomp;
	push @gold,$_;
    }
    $cc=1;
}
if ($cc==0){
    @gold=("ARRB1", "ARRB2", "SLC4A1", "CFTR", "CHMP1A", "CRYAA", "CRYAB", "CYCS", "DAB1", "DAB2", "APPL1", "APPL2", "EP15_HUMAN", "EPS15L1", "EPN1", "ERCC2", "FGF1", "FGF2", "GAPDH", "GPI", "GPHN", "GK", "GPX4", "GSK3B", "ALAD", "HIP1", "XRCC6", "XRCC5", "LGALS1", "LGALS2", "LGALS3", "LONP1", "ABCB1", "NRP1", "PGK1", "PCBD1", "PLK1", "PLOD3", "PMS2", "PRDX6", "SMC3", "TSG101", "TMSB4X", "TYMP");
}
while(<>){
    if(/\[CLASS:\s*(\d+)/){
	$class=$1;
    }
    elsif(/^CA\s*(.+)/){
	$annot=$1;
    }
    elsif(/^P\#\s*(\d+)/){
	$num=$1;
    }
    elsif(/^PN\s*(.+)/){
	my @prots=split(/,\s+/, $1);
	map{
	    $prot{ANNOT}{$_}=$annot;
	    $prot{NUM}{$_}=$num;
	    $prot{CLASS}{$_}=$class;
	}@prots;
    }
}
print "max : $max, min: $min\n";
map{
    if ((defined($prot{CLASS}{$_})) && ($prot{NUM}{$_}>=$min) && ($prot{NUM}{$_}<=$max)){
	print "\nPROTEIN $_\n";
	print "ANNOTATION $prot{ANNOT}{$_}\n";
	print "NUMBER $prot{NUM}{$_}\n";
	print "CLASS $prot{CLASS}{$_}\n------------------------\n";
    }
    else{
	print "$_ was not present in the original network\n" if $opts{a};
    }
   
}@gold;


sub usage{
    
    print  <<EndOfHelp;

USAGE:
       parse_class.pl [-p protein list] [-f protein file] <CLASS FILE>

This script will parse a class file and search for each of the proteins
given. Protein names can be given in a file (one per line) or with the 
optional -p flag, as a comma separated list in the command line. If both
are given, the -p proteins will be added to the file. If no file or 
proteins are given, the gold standard set is used. 

COMMAND-LINE OPTIONS:
    -a : Print output for ALL proteins given, even if they were not present
         in the original network.
    -f : File containing the names of the proteins of interest, one per line.
   
    -h : Print this help and exit.
    -m : Min class size. Only return info for those classes with >than this
         number of members
    -M : Max class size. Only return info for those classes with <than this
         number of members

    -p : Comma separated list of proteins of interest.

EndOfHelp

    exit();
}
