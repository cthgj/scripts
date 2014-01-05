#!/usr/bin/perl -w

use strict;
my $onto="$ARGV[0]" || 'P';
my $n=$ARGV[1]||die("need a list of names as ARGV1\n");
my %k;
my %codes =(
    "EXP" => "Experimental",
    "IDA" => "Experimental",
    "IPI" => "Experimental",
    "IMP" => "Experimental",
    "IGI" => "Experimental",
    "IEP" => "Experimental",
    "ISS" => "Computational",
    "ISO" => "Computational",
    "ISA" => "Computational",
    "ISM" => "Computational",
    "IGC" => "Computational",
    "IBA" => "Computational",
    "IBD" => "Computational",
    "IKR" => "Computational",
    "IRD" => "Computational",
    "RCA" => "Computational",
    "TAS" => "Experimental",
    "NAS" => "Experimental",
    "IC" => "Computational",
    "ND" => "Computational",
    "IEA" => "Computational"
    );

my $c=0;
my %found;

my @names;
open(N,"$n")||die("Could not open $n:$!\n");
while(<N>){
    chomp;
    push @names,$_;
}

open(A, "gawk -F\"\t\" '{if(\$9==\"$onto\" && \$4==\"\"){print \$2,\$7}}' $ARGV[2] |")||die("Can't open $ARGV[2]:$!\n");
while(<A>){
    $c++; 
    /(.+)\s+(.+)/;
    $k{$1}{$codes{$2}}++;
#    $found{$1}++;
}    
my $none=0;
my ($exp,$com)=(0,0);
foreach my $prot (@names){
    if (defined($k{$prot}{"Experimental"})){
	$exp++
    }
    elsif( defined($k{$prot}{"Computational"})){
	$com++
    }
    else{$none++}
}
print "Experimental:\t";
print $exp*100/$#names . " ($exp)\n";
print "Computational:\t";
print $com*100/$#names . " ($com)\nNONE: " . $none*100/$#names . " ($none)\n";


# map{print "$_ $k{$_}\t" . $k{$_}*100/$c . "\n"}keys(%k);
# #my %k; my $c=0; while(<>){$c++; /(.+)\s+(.+)/; $k{$2}++} map{print "$_ $k{$_}\t" . $k{$_}*100/$c . "\n"}keys(%k);
# my $a=$c-scalar(keys(%found));
# print "Missing : $a\t" . $a*100/$c   . "\n";
