#!/usr/bin/perl -w

use strict;
use Getopt::Std;

my $c = 0; ## these are both just counters to know where I am inside the file
my $d = 0;

my %opts;
getopts('anfAl:e:',\%opts) || do { print "Invalid option\n"; exit(1); };

my $n          = $opts{n} || 0; ## print names
my $a          = $opts{a} || 1; ## print all alignments (default)
my $EVALUE     = $opts{e} || 1; ## maximum evalue
my $list       = $opts{l} || 0; ## print only output for list of names
my $fasta      = $opts{f} || 0; ## return fasta output
my $aln        = $opts{A} || 0; ## return aligned FASTA output

$a = 0 if $n;
$a = 0 if $fasta;
$fasta = 1 if $aln;


my($evalue,$name);
my (@alignment, @names, @keys);
my (%hit, %ev);
while(<>)
{
    if (/Alignments/){
	$c=1;
	next;
    }
    if(/Histogram/){
	$d=1;
	next;
    }
    next unless $c==1 && $d==0;;
    if (/^[\S]/)
    {
	if($name)
	{
	 
	 $hit{$name} = [@alignment];
	 @alignment = ();
	 
	}
	
	/^(.*?):.*E\s=\s(.*)\n/ || die("cannot match : $_\n");
	$name = $1;
	$evalue = $2;
	if ($evalue =~ /^\s?e/)
	{
	    $evalue = "1" . $evalue;
	}
	push @alignment, $_; 
	$ev{$name} = $evalue; 
	
	if($n)
	{
	    print STDOUT "$name\n" if 	$ev{$name} <= $EVALUE;
	}
	next;
    }
    push @alignment, $_ unless /Alignme/; 
    
}
@keys = keys(%hit) unless $list;
my @hits;

map{ push @hits, $_ if $ev{$_} <= $EVALUE} @keys;  ## get the good hits

&list_output() if ($list) ;

&all_alignments_output if($a);

&fasta_output if($fasta);


sub list_output
{
    open(LIST, "$list");

    while(<LIST>)
    {
	/^(.*)\n/;
	push @keys, $1;
    }
    close(LIST);
}

sub all_alignments_output
{
    foreach my $gi (@hits)
    {
	print STDOUT "@{$hit{$gi}}\n";
	
    }
    
}

sub fasta_output
{
    foreach my $gi (@hits)
    {
	my $seq;
	my $len = @{$hit{$gi}};
	print STDOUT ">@{$hit{$gi}}[0]";
	for(my $n=1; $n<$len; $n++)
	{
	    next if ${$hit{$gi}}[$n] =~ /^\s+$/;
	    if(${$hit{$gi}}[$n]=~ /\s\d+\s+([^\d]*)\s/)
	    {
		$seq = $seq . $1;
		$seq =~ s/-//g unless $aln;
	    }
	}
	my $l = length($seq);
	my $k = 0;
	while($k<$l)  
	{
	    print STDOUT substr($seq,$k,60) . "\n"; 
	    $k=$k+60;
	}
    }
    
}    
