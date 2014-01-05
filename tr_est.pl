#!/usr/bin/perl -w

# This script takes a .aln file (with frame info) and a FASTA file of the subject sequences.
# It will return one multi FASTA (called subjects.pep.fa, cant STDOUT because I'm using transeq
# to translate) with the translated subjects in the right frame.


# USAGE : tr_est FILE.aln FILE.fa outfile_name

use Getopt::Std;
use strict;

my ($fr, $sbjct);
my %frame = ();
my %seq   = ();
my @sequence;
my $gi = undef;
my %opts;

getopts('on::',\%opts);
my $outfile = $opts{o} || die("You must specify an output file name with the -o flag\n");
my $namelist = $opts{n}

open(HA,"$ARGV[0]") || die("$ARGV[0] : $!\n");

# Get the frame of each target seq from
# the .aln file
my $n;
while(<HA>)
{
    if(/Frame/ || /^>/)
    {
	if (/>/)
	{
	    />(.*?)\s/;
	    $n=$1;

	}
	else
	{
	    /Frame\s=\s(.*?)\s/;
	    $frame{$n} = $1;
	}
    }
}

open(SEQ,"$ARGV[1]");

my $fname;
my %fnames = ();
# Get sequences
while(<SEQ>)
{
    next if /^\s*$/;
    if (/^\#\#/) ## special case used on pseudo fasta I made (never mind)
    {
	/\#\#(.*)/ || die("shit : $.\n");
	$fname = $1; 
    }
    elsif (/>/) 
    {
	if ($gi)
	{
	    $seq{$gi} = [@sequence];
	    @sequence = ();
	    $fnames{$gi} = $fname || "undef";
	   # print STDOUT "fname : $fname\n";
	}

	/>(.*?)\s/;
	$gi = $1;

    }
    else
    {
	push @sequence, $_;
    }
}
close(SEQ);

# deal with last sequence
$seq{$gi} = [@sequence];
$fnames{$gi} = $fname || "undef";
#print STDOUT "fname : $fname\n";

#print STDERR ".";

my @names = keys(%seq);
my $hit;

#while ((my $ka,my $v) = each %fnames) { $k{$v}++; } my @a = keys(%k); map{print STDOUT "$_\n"} @a; die();    #my $k = @a; print STDERR "k : $k\n";die();

# Print output
while ((my $k,my $v) = each %frame )
{

#    print STDERR "$k";
    open(TMP, ">$$");
    my $ho = $k;
    $k =~ s/\|/_/g;
    ($hit) = grep {$_ =~ /$ho/;} @names;

    ## Get out file name and revert to -o value if not found
    my $outfile_name = $fnames{$hit} . "\.nfa" || $outfile;
    
    print STDERR "fnames{$hit} : $fnames{$hit}\n";
    print STDERR "hit : $hit --- outfil : $outfile_name\n";
    print STDERR "seq{$hit} : @{$seq{$hit}}\n"; 
    
    if (defined($seq{$hit})) {
    print TMP ">$k\n";
    map { print TMP "$_";} @{$seq{$hit}}; 
}
    close(TMP); 
    if ($v =~ /\+/) { $v =~  /\+(\d)/; $v = $1;}  ### no need to do same with -, transeq accepts it
    
    print  "transeq -frame $v -outseq $$.pep $$ : $k\n";
    system("transeq -frame $v -outseq $$.pep $$; cat $$.pep >> $outfile_name ; rm $$ $$.pep; sed 's/_/|/g' $outfile_name > $$; mv $$ $outfile_name" );
    print STDERR "#######\n";

#    print "$k => $v\n";
}
