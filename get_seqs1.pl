#!/usr/bin/perl -w

# This script takes as input a list of FASTA ID lines
# and read a sequence file containing these sequences from STDIN. It will print out those seqs
# of the seq file whose ids were in the id file 
# will work on rel small FASTA files... retrieveseq.pl better


# USAGE  zcat/cat FASTA input file | get_seqs.pl name.file > outfile

#


use strict;
use Getopt::Std;
my %seqs = ();
my %gis = ();
my %names = ();



my $match = 0;
my $name;

my @seq;


my %fasta = ();
my %opts;
getopts('f',\%opts);
my $fastaseq = $opts{f} || undef;

#&usage unless $ARGV[0];
#my @names = qw(
#gi_14826977_dbj_AU199676.1_AU199676_1  gi_275187_gb_M75834.1_M75834_1        gi_31237652_dbj_AU112553.2_AU112553_3 gi_14831523_dbj_AU202073.1_AU202073_2  gi_30730122_gb_CB388412.1_CB388412_6  gi_31246912_dbj_BJ125791.2_BJ125791_3 gi_18246813_dbj_BJ104143.1_BJ104143_1  gi_30734090_gb_CB392379.1_CB392379_2  gi_4719645_gb_AI641170.1_AI641170_4 gi_18281409_dbj_BJ121276.1_BJ121276_3  gi_30743290_gb_CB401563.1_CB401563_2  gi_7742956_dbj_AV413780.1_AV413780_3
#);

my @names = "gi_275187_gb_M75834.1_M75834_1";

map
{
    open(NAMES,"$_");
    my @n = <NAMES>;
    $names{$_}= [@n];   # $names{name} = list of names
    close(NAMES);
} @names;

print STDERR "got names\n";

# get sequences


my @sseq;
my $seq_name;
my @nombres = keys(%names);


&fasta if ($fastaseq);
&tbl unless ($fastaseq);



####
# Functions
####




sub fasta
{
    while(<>)
    {
	if(/^>/)
	{
	    
	    if ($seq_name) ## if NOT first seq
	    {
		$seq_name =~ s/\|/_/g;
		$fasta{$seq_name} = [@sseq];
		@sseq = ();
	    }
	    />(.*?)\s/;
	    $seq_name = $1;
	}
	else
	{
	    push @sseq, $_;
	}
    }




### Get last sequence 
    $seq_name =~ s/\|/_/g;
    $fasta{$seq_name} = [@sseq];
    @sseq = ();
    print STDERR "\n***********\ngot hash\n***********\n";

    my @keys = keys(%fasta);

     ### each $names{nombre} is a list of names
    
    foreach my $nombre (@nombres)  # foreach list of sequence names
    {
	chomp $nombre; 
	
	map {
	    my $curr = $_;  # curr = partial name in file 
	    chomp $curr;
	    chop $curr;
	    $curr =~ s/\|/_/g;  
	    my @fin_name = grep {$_ =~ /$curr/;} @keys;  # find key which matches partial name
	    open(NAM,">>$nombre\.fa");
	    print NAM ">$_";
	    map {print NAM  "$_";} @{$fasta{$fin_name[0]}};
	    close(NAM);
	} @{$names{$nombre}};
    }
    
}

sub tbl
{
    my %tbl = ();
    while(<>)
    {
	/^(.*?)\s(.*)/;
	my $tmp = $1;
	my $seq = $2;
	$tmp =~ s/\|/_/g;
	$tbl{$tmp} = $seq;

    }
    print STDERR "\n***********\ngot hash\n***********\n";
    my @keys = keys(%tbl);
    foreach my $nombre (@nombres)  # foreach list of sequence names
    {
	chomp $nombre; 
	
	map {
	    my $curr = $_;  # curr = partial name in file 
	    chomp $curr;
	    chop $curr;
	    $curr =~ s/\|/_/g;  
	    my @fin_name = grep {$_ =~ /$curr/;} @keys;
	    open(NAM,">>$nombre\.fa");
	    print NAM ">$_";
	    my $l = length($tbl{$fin_name[0]});
	    my $i = 0;
	    while($i<$l)
	    {
		print NAM substr($tbl{$fin_name[0]},$i,60) . "\n"; 
		$i=$i+60;
	    }
	}  @{$names{$nombre}};
    }
}

