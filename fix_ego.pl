#!/usr/bin/perl -w


#  This script will take the EGO database files and return TOG# seqnames #

use strict;

my %seq;
my %orth;
my %tog;
my $togno;
my @array;


open(ORTH, "$ARGV[0]") || die "1 cannot open orth: $!";  ### Open orth file ARGV[1]

while (<ORTH>)
{ 
    next if /^\s+$/;
    if(/^\>(TOG\d{5,7})/)
    {
	$tog{$togno} = "@array" if $togno;
	@array=();   ### Grab new TOG and set array to 0
	$togno = $1; ###

    } 
    else
    {
	/^(.*?)\|(.*?)\s/;
	my $c = "$1" . '_' . "$2";
	push @array, $c; ### Push the togno and ORTH line into the array
    }
    $tog{$togno} = "@array"; ### now $tog{togno} will return the TOG# and description of the sequence
   
    
}
close(ORTH);

my @a = keys(%tog);
map{ print STDOUT "$_ $tog{$_}\n"; } @a;

sub get_seq
{
    open(SEQ,"$ARGV[1]") || die "0 cannot open seq: $!";  ### Open seq file, ARGV[0]
    while (<SEQ>)
    {
	/^(.*?)\s(.*?)$/;
	$seq{$1} = $2; ### $seq{seqname} = sequence
    }
    close(SEQ);   
    print STDERR "########################\n       Got SEQ\n########################\n";
}

sub fasta  ### Was necessary in previous EGO version not so now. Keep just in case
{
    
    my @temp = %tog;
    open (EGO, ">/home/ug/cchapple/research/seq/EGO/ego4_060903.tmp")|| die "2 cannot open ego: $!";
    print EGO "@temp";
    close(EGO);
    print STDERR "########################\nOpened EGO\n########################\n";
    open (EGO, "/home/ug/cchapple/research/seq/EGO/ego4_060903.tmp")|| die "3 cannot open ego: $!";
    open (FA, ">/home/ug/cchapple/research/seq/EGO/ego4_060903.fa")|| die "4 cannot open FA: $!";
    while (<EGO>)
    {
	if (/TOG\d{5,7}\sTOG\d{5,7}/)
	{
	    s/(TOG\d{5,7})\sTOG\d{5,7}/$1/;
	}
	my $la = $_;
	chomp $la;
	$la =~ s/^\s+//;  ### Do all this to remove spaces from some of the strings
	
	/^.?TOG\d{5,7}\s(.*?)\s/;
	
	print STDERR "EGO \$1 is $1\n";
	
	print FA ">$la\n$seq{$1}\n";
    }
    
    
}
