#!/usr/bin/perl 
use Getopt::Std;
use strict;
my %opts;
getopts('sb:h',\%opts);


sub usage{

    $0=~/.+?([^\/]+)$/;
    my $p_name=$1;
    print STDERR <<EOF;
DESCRIPTION:
  This script will search a psicquic file for a given list of interactions. It takes 4 files as input. A .blastmap (optional), a .psi2id , a list of the interactions in question and a psiqcuic file (compressed or uncompressed).

USAGE:
  $p_name [-b blastmap] <psi2id> <psiqcuic> <interaction_list>

OPTIONS:
    -s  :  For each interaction, only print whether it was found or not and the line number 
           of the psicquic file line(s) that corresponds to the interaction. Without
           this option, the entire psicquic line will be printed for each interaction.
    -h  :  Print this message and exit.

EOF
exit;
}
usage() unless ($#ARGV==2) ;
usage() if $opts{h};
my $blastmap=$opts{b}||undef;

my $psi2id=$ARGV[0]||usage();
my $psifile=$ARGV[1]||usage();
my $inter_file=$ARGV[2]||usage();


## Parse blastmap
my %bmap;
my $fh;
if ($blastmap) {
    open($fh, "<", "$blastmap")||die ("Could not open $blastmap for reading : $!\n");
    while (<$fh>) {
	chomp;
	my @a=split(/\s+/); 
	$bmap{$a[0]}=$a[1];
    }
close($fh);
}

## Load missing interactions
my (%inters, %out, %prots, %orig);
open($fh, "<", $inter_file)||die ("Could not open $inter_file for reading : $!\n");
while (<$fh>) {
    chomp; 
    my @a=split(/\s+/);
    my $pair=join("\t", sort {$b lt $a} @a);
    $prots{$a[0]}=$prots{$a[1]}=1;
    $inters{$pair}++;
    $out{$pair}{STATUS}="MISSING";
}
close($fh);

## Load psi2id info
my %map; 
open($fh, "<",  "$psi2id")||die ("Could not open $psi2id for reading : $!\n");
while(<$fh>){
    chomp; 
    my @b=split(/\t/); 
    $map{$b[0]}=$b[1];# if defined($prots{$b[1]});
    $map{$b[0]}=$b[1];# if defined($prots{$bmap{$b[1]}});
} 
close($fh);


## Look for the missing interactions in the psicquic file
if ($psifile=~/\.gz$/ || $psifile=~/\.Z$/) {
    open($fh, "-|", "zcat $psifile")
	||die ("Could2 not open $psifile for reading : $!\n");     
}
else {
    open($fh, "<",  $psifile)
	||die ("Could1 not open $psifile for reading : $!\n");    
}


while(<$fh>){
    chomp;
    my @a=split(/\t/); 
    ## Translate psifile name to network name
    my $p1=$map{$a[0]}; 
    my $p2=$map{$a[1]};
    my $o1=$p1;
    my $o2=$p2;
#    print "PP $p1 : $p2\n";
    ## If this prot has been mapped to another in the blast step
    ## get the mapped name
    defined($bmap{$p1}) && do {	$p1=$bmap{$p1}};
    defined($bmap{$p2}) && do {$p2=$bmap{$p2}};
    my $pair=join("\t", sort {$b lt $a} $p1,$p2);
    $orig{$pair}{$p1}=$o1;
    $orig{$pair}{$p2}=$o2;
    
    if ($opts{s}) {
	if (defined $inters{$pair}){
	    $out{$pair}="FOUND ($.)";
	}
    }
    else {
	print "$p1 ($orig{$p1}) $p2 ($orig{$p2})\t$_\n" if defined $inters{$pair};	
    }

}

if ($opts{s}) {
    map{
	my @a=split(/\t/,$_);
	# $bmap{$a[0]}=$a[0] unless defined($bmap{$a[0]});
	# $bmap{$a[1]}=$a[1] unless defined($bmap{$a[1]});
#	my @b=split(/_/, $out{$pair}{ORIG});
	
	print "$a[0] ($orig{$_}{$a[0]})\t$a[1] ($orig{$_}{$a[1]}) : $out{$_}\n"

    }keys(%inters);
}
