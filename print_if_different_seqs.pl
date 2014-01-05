#!/usr/bin/perl -w
use strict;

my %k; 
my (@a,@b,@a1,@b1);
open(A,"$ARGV[0]")||die("bad filename: $ARGV[0]\n"); ; 
while(<A>){
    chomp;
    $k{$_}{NAME}=$_;
} 
close(A);
open(F,"$ARGV[1]");
while(<F>) {
    @a=split(/\t/,$_); 
    my $seq=$a[1]; 
    $a[0]=~s/^..\|//; 
    @b=split(/\|/,$a[0]); 
    foreach my $name (@b){
	if (defined($k{$name}{NAME})){
	    $k{$name}{SEQ}=$seq;
	}
    }
}
close(F);
open(B,"$ARGV[2]")||die("bad filename: $ARGV[2]\n"); 
while(<B>){
    @a1=split(/\t/,$_); 
    my $seq1=$a1[1]; 
    $a1[0]=~s/^..\|//; 
    @b1=split(/\|/,$a1[0]); 
        foreach my $name (@b1){
	    if (defined($k{$name}{NAME})){
		if (defined($k{$name}{SEQ})){
		    print $k{$name}{NAME} . "\n" unless $seq1 eq $k{$name}{SEQ}
		}
		else{
		    print "$k{$name}{NAME}\n";
		}
	    }
	}
} 
