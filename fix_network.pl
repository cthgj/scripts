#!/usr/bin/perl -w

# This script will consolidate a network file. If two synonyms of a given
# protein are present in the network, it will keep the most used synonym
# of the two, and consolidate their interactions.  

print STDERR "not working yet\n"; exit();
use strict;
use Getopt::Std;
use IO::File;
my %k;
my $sfile=$ARGV[0]||die("need a annotation file (GAF)");
my $nfile=$ARGV[1]||die("need a network file");

my %is_synonym;
my $fh = IO::File->new("< $nfile")|| die("Cannot open $nfile : $!\n"); 
while(<$fh>){
    chomp;
    next if /^\d+$/;
    my ($a,$b)=split(/\t/,$_);
    $k{$a}++;
    $k{$b}++;
}
close($fh);
$fh = IO::File->new("< $sfile")|| die("Cannot open $sfile : $!\n"); 
## read gaf file
my %oo;
while(<$fh>){
    my $silivar="!";
    next if /!/;
    next if /^\s+$/;
    my @tmparray=split(/\t/);
    my @names=split(/\|/,$tmparray[10]);
    push @names,$tmparray[1],$tmparray[2];
    for(my $n=0;$n<scalar(@names);$n++){
	if(defined($k{$names[$n]})){

	}
	# for(my $k=scalar(@tmparray)-1;$k>-1;$k--){
	# 	  $is_synonym{$tmparray[$n]}=$tmparray[$k];
	# }
    }
    
}
close($fh);



$fh = IO::File->new("< $nfile")|| die("Cannot open $nfile : $!\n"); 
while(<$fh>){
    chomp;
    next if /^\d+$/;
    my ($a,$b)=split(/\t/,$_);
    $k{$a}++;
    $k{$b}++;
}
close($fh);

