#!/usr/bin/perl -w


use Getopt::Std;
use strict;


my $has_secis;
my $is_gene=0;
my @gene;
my @last_gene;
my $gene_type = "hh";
while(<>)
{
    next if /^\s+$/;
    if(/^\s?\#/)
    {
	$is_gene = 0;	
    }
    if(/>/)
    {
#	s/>\s/>/;
	@gene = ();
	push @gene, $_;
	$is_gene =1;
	next;
    }
    if ($is_gene == 1)
    {
	push @gene, $_ unless /^\s+$/;
    }
    if(/AA\-/)
    {
	if(/\sSECIS/)
	{
	   map{print "$_"} @gene unless $gene_type =~ /SECIS/;
	}
	/^\s*([^\s]+?)\s/;
	$gene_type=$1;
	
    }
    

