#!/usr/bin/perl -w
# This script should be run on output of AlignCys2ESTs_frame.awk. It will open a list of seq names (ARGV[0], one 
# name per line), look through the .aln file to see what frame their hsp was in and then return a list of the 
# names and their frames (default) to STDOUTor trasnlate the appropriate frame into $seq_name.pep.fa (-f). It 
# assumes the presence of FASTA sequence files of each name in the list, of the format name.fa
############################-NOTE-#############################
#                                                           #
# This is a very specific script, relying on certain files  #
# read/modify the source before using.                      # 
#                                                           #
############################-NOTE-#############################
# Usage : frame.pl [-f] list_of_names input_file.aln



use strict;
use Getopt::Std;


my $name;
my %frame = ();
my %pep = ();
my %dels = ();

my %opts;
getopts('f',\%opts) || die("bad option\n"); 

my $f = $opts{f};

open(K, "$ARGV[1]");

while(<K>) ## go trhough .aln file to get frames
{
    next unless ( /Frame/ || /Query=/);  ## if want subjects use />/
    if (/Query=/)
    {
	/Query=\s(.*?)\s+\(/;
	$name = $1;
    }
    else
    {
	/Frame\s=\s(..)/;
	my $fr = $1;
	if ($fr =~ /\+/) { $fr =~ /\+(\d)/; $frame{$name} = $1;}
	else { $frame{$name} = $fr;}  ### frame{est} = est's frame 
    }
    
}
close(K);
my @frames = keys(%frame);

&get_fasta();  ## get fasta seq of each est


open(NAMES, "$ARGV[0]") || die("error cannot open $ARGV[0] : $!\n");

my @ests = <NAMES>;
my @peps = keys(%pep);    

if ($f)
{
    foreach my $pep (@peps)
    {
	foreach my $nombre (@{$pep{$pep}})
	{
	   chomp $nombre;
	  # print STDERR "transeq -frame $frame{$nombre} \"$nombre\.fa\" -outseq tmp.fa\n";
	   system("transeq -frame $frame{$nombre} \"$nombre\.fa\" -outseq tmp.fa") && die("$nombre : $!\n");
	   
	   ### Put gi names back to original to avoid confusion (- back to |) and cat into file
	   system("cat tmp.fa | sed \'s/-/|/g\' >> $pep.pep.fa "); 
	   $dels{$nombre}++
	   	 
	}
	

    }
}
else
{
    map { print STDOUT "$_ $frame{$_}\n";} @frames;
}




unlink("tmp.fa");
my @dk = keys(%dels);

map {
    unlink("$_.fa")
    } @dk;


### Get FASTA sequence foreach set of ESTs
sub get_fasta
{
    open(PEP,"wormpep.fa_all_nematode_ests_5ormore.strict");
    while(<PEP>)
    {
	/^(.*?)_.*?\s(.*)$/;
	my $nam = $1; 
	my $g = $2; 
	$g =~ s/_[^_]*?_[+-]\s?/ /g; ### the [^_] class is necessary for certain sequence names
	my @a = split/\s/, $g;       ### @a : list of EST names for this pep
	$pep{$nam} = [@a];           ### pep{peptide_name} = list of ests
	map
	{
	    
	    system ("FastaToTbl ../wormpep.fa_all_nematode_ests_5ormore.strict.est.fa | grep \"$_\" | sed \'s/|/-/g\' |  TblToFasta > \"$_\".fa") && die("$_ : $!\n");  ### pass through sed to change gi anmes (| confuses things)
	} @a; 
    }
    return %pep;
}


