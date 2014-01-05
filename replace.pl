#!/usr/bin/perl -w

use strict;
use Getopt::Std;
my (%opts,%from);
getopts('f',\%opts);


unless($ARGV[0]){
    print STDERR "USAGE: $0 <replacements> <input_file> <output_file>\n\n";
    exit();
}

die("Need an output file\n")unless $ARGV[2];

open(A,"$ARGV[0]");
open(B,"$ARGV[1]");
open(C,">$ARGV[2]")||die("Need an output file\n");
while(<A>){
    chomp;
    my @a=split(/\t/);
#    $a[1]=$a[0] if $a[1]=~/^\s*$/;
    $a[1]=$a[0] unless $a[1];
    die("$a[0] is already defined as $from{$a[0]}\n") if defined($from{$a[0]});
    $from{$a[0]}=$a[1];
}
while(<B>){
 if($opts{f}){
     s/(ENSG.+?)\s/$from{$1} /g;
 }

 else{
     foreach my $pat (keys(%from)){
	 s/$pat/$from{$pat}/g;
     }
 }
print C;

}
