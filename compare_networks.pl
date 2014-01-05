#!/usr/bin/perl -w
use strict;
use Getopt::Std;
my %opts;
getopts('md', \%opts);

my $map_file=$opts{M}||undef;
my $net1=$ARGV[0] ||&usage();
my $net2=$ARGV[1] ||&usage();
my (%map,%interactions, %prots);
my $missing=$opts{m}||1;
my $dupes=$opts{d}||undef;

$missing=undef if $dupes;

&usage() unless $missing || $dupes;


if ($map_file) {
    %map=%{get_map()};
}

open(B,"$net1")||die("Cannot open $net1: $!\n");
my @net1_names;
while(<B>){
    chomp;
    next if /^\d+$/;
    my @fields=split(/\t/,$_);
    $fields[0]=$map{NAME}{$fields[0]} if defined($map{NAME}{$fields[0]});
    $fields[1]=$map{NAME}{$fields[1]} if defined($map{NAME}{$fields[1]});
    my $pair=join("\t", sort(@fields));
    $interactions{NET1}{$pair}++;
    $prots{NET1}{$fields[0]}=$prots{NET1}{$fields[1]}=1;
}
close(B);

open(C,"$net2")||die("Cannot open $net2: $!\n");
my (@MISSING, @DUPES, %dprots);
while(<C>){
    chomp;
    next if /^\d+$/;
    my @fields=split(/\t/,$_);
    $fields[0]=$map{NAME}{$fields[0]} if defined($map{NAME}{$fields[0]});
    $fields[1]=$map{NAME}{$fields[1]} if defined($map{NAME}{$fields[1]});
    my $pair=join("\t", sort(@fields));
    $interactions{NET2}{$pair}++;
    $prots{NET2}{$fields[0]}=$prots{NET2}{$fields[1]}=1;
    unless (defined($interactions{NET1}{$pair})){
	    push @MISSING,$pair;
    }
    if (defined($prots{NET1}{$fields[0]})){
	    $dprots{$fields[0]}=1;
    }
    if (defined($prots{NET1}{$fields[1]})){
	    $dprots{$fields[1]}=1;
    }

}
close(C);

my @aa=keys(\%{$interactions{NET1}});
my @bb=keys(\%{$interactions{NET2}});
my $pairs1=$#aa+1;
my $pairs2=$#bb+1;

my @pairs=(@aa,@bb);
my (@net1, @net2);
foreach my $pair (@pairs) {
    if (defined($interactions{NET1}{$pair}) && defined($interactions{NET2}{$pair})) {
		push @DUPES,$pair;
    }
    else {
	if (defined($interactions{NET1}{$pair})){
	    push @net1, $pair;
	}
	elsif (defined($interactions{NET2}{$pair})){	    
	    push @net2, $pair;
	}
	else {
	    die("Something weird here...\n");
	}
    }
}
my ($s1,$s2)="";

## This is such a waste of time, but I like clean formatting :)
my @p1s=keys(%{$prots{NET1}});
my @p2s=keys(%{$prots{NET2}});
my @dps=keys(%dprots);
my @string1s=($ARGV[0], $ARGV[1], "Unique interactions", "Only in $ARGV[0]", "Only in $ARGV[1]", "Common interactions", "Prots in $ARGV[0]", "Prots in $ARGV[1]", "Common Prots");
my @string2s=($pairs1,$pairs2, $#net1+1+$#net1+1+$#DUPES+1,$#net1+1,$#net2+1, $#DUPES+1, $#p1s+1, $#p2s+1, $#dps+1);

my $longest=0;
map{
    $longest=length($_) if $longest < length($_);
}@string1s;

for (my $i=0; $i<=$#string1s; $i++) {
    my $s;
    my $xx=length($string1s[$i]);
    while (length($string1s[$i])<$longest) {
	$string1s[$i].=" ";
#	$xx=length($string1s[$i]);
    }
    print "$string1s[$i]\t:\t$string2s[$i]\n";
}

# $s1="\t" if length($ARGV[0])<20;
# $s2="\t" if length($ARGV[1])<20;
# print "$ARGV[0]\t\t$s1:\t$pairs1\n";
# print "$ARGV[1]\t\t$s2:\t$pairs2\n";
# print "Unique interactions\t\t:\t" , $#net1+1+$#net1+1+$#DUPES+1 , "\n";
# print "Only in $ARGV[0]\t\t\t:\t" , $#net1+1 , "\n";
# print "Only in $ARGV[1]\t:\t" , $#net2+1 , "\n";
# print "Common interactions\t\t:\t" , $#DUPES+1 , "\n";



sub get_map{
    my %map;
    open(A,"$map_file")||die("Cannot open $map_file: $!\n");
    while (<A>) {
	chomp;
	my @fields=split(/\t/,$_);
	##what kind of map is it?
	if ($#fields==6) {
	    $map{ACC}{$fields[0]}=$fields[2];
	    $map{ACC}{$fields[1]}=$fields[3];
	    $map{NAME}{$fields[2]}=$fields[0];
	    $map{NAME}{$fields[3]}=$fields[1];
	} elsif ($#fields==1) {
	    $map{ACC}{$fields[0]}=$fields[1];
	    $map{NAME}{$fields[1]}=$fields[0];
	
	} elsif ($#fields==2) {
	    ## Oldname\tNEWNAME\tAC
	    $map{ACC}{$fields[0]}=$fields[2];
	    $map{ACC}{$fields[1]}=$fields[2];
	    $map{NAME}{$fields[0]}=$fields[1];
	    $map{NAME}{$fields[2]}=$fields[1];
	    
	} else {
	    die("unkown map file format\n");
	}
    }
    close(A);
    return(\%map);
}

sub usage{
    print STDERR <<EndOfHelp;

  USAGE: compare_networks.pl [-M map] net1 net2

This script takes a map file (of various formats) and two networks. 

   -m : return interactions present in net2 and missing in net1
   -d : return interactions present in both net1 and net2

EndOfHelp
      exit(0);
}
