#! /usr/bin/perl -w

# This script will parse geneid output, select those sequences with >1  exon and a SECIS element. In order to get these sequences easily,
# the script will rerun geneid on these sequences. Assuming, that is that you have a directory with all the fasta sequences...





use strict;
use Getopt::Std;
my $genes = 0;
my $sequence;
my $exons = 0;
my @sequences;
my %seq =();


my %opts;
getopts('gh',\%opts) || do { print "Invalid option, options are -h and -g\n"; exit(1); };

my $outfile = $ARGV[1] || 'STDOUT';
my $error   = $ARGV[2] || 'STDERROR';

#print STDERR "\$ARGV[1] = $ARGV[1]\n \$ARGV[2] : $ARGV[2]\n\$ARGV[0] : $ARGV[0]\n";  die();



my $geneid = $opts{g} || 0;
my $help = $opts{h} || 0;


&usage if $help;
    &usage unless $ARGV[0];
open(FILE,"$ARGV[0]") || die ("cannot open file : $!\n");

while(<FILE>)
{   	
#    s/CONTIG/contig/;
    if (/Sequence\s(.*?)\s/)
    {
#	print STDERR "@";
	$sequence = $1;
	next;
    }
    if (/2\sgenes/)
    {
#	print STDERR "!";
	$genes = 1;
	next;
    }
    if (/(\d)\sexons/)
    {
	$exons = 1 if $1 > 1;
#	print STDERR "*";
	next;
    }

    if(($genes == 1) && ($exons == 1))
    {
#	push @sequences, $sequence;
	$seq{$sequence}++;
	$genes = 0;
	$sequence = '';
	$exons = 0;
    }
    else 
    {
	$genes = 0;
	$sequence = '';
	$exons = 0;
    }
  #  print STDERR "*****\$genes : $genes\n\$sequence : $sequence\n\$exons = $exons\n*****";

}
close(FILE);

@sequences = keys(%seq);

&geneid if $geneid;
map{ print STDOUT "$_\n"} @sequences unless $geneid;


sub geneid{
    map
    {    
	system( "~scaste/geneid/bin/geneid -vP /home/ug/cchapple/research/selenoproteins/tetraodon/SECISearch/potential/tetraodon.param.3.No_SECIS.cDNA_full_length.param -R /home/ug/cchapple/research/selenoproteins/tetraodon/SECISearch/potential/gff/twil/$_\.gff /home/ug/cchapple/research/selenoproteins/tetraodon/SECISearch/potential/fasta/twil/$_ >> $outfile 2>> $error");

print STDERR "***********************************************\n~scaste/geneid/bin/geneid -vP /home/ug/cchapple/research/selenoproteins/tetraodon/SECISearch/potential/tetraodon.param.3.No_SECIS.cDNA_full_length.param -R /home/ug/cchapple/research/selenoproteins/tetraodon/SECISearch/potential/gff/twil/$_\.gff /home/ug/cchapple/research/selenoproteins/tetraodon/SECISearch/potential/fasta/twil/$_ >> $outfile 2>> $error\n***********************************************\n";
#	die();

    } @sequences;
   
}




sub usage
{
   print STDOUT <<"EOF";
USAGE:
   parsegenid.pl <geneid.out> <new_geneid.out> <new_geneid.error>
EOF

exit;
}
