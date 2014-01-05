#!/usr/bin/perl -w
use strict;

open(R,"ARGV[0]");
open(H,"ARGV[1]");
my %found;
while(H){
    @a=split;
    $found{$a[0]."xxx".$a[1]}++;
}

my @aa;
while(<R>){
    if((/^hhyper/) || (/^lhyper/)){
	push @a;
    }
    elsif(/^hypers/){
	/\"(.+?)\",\"(.+?)\"/||die("cannot match\n");
	unless(defined($found{$a[0]."xxx".$a[1]}) || 
	       defined($a[1]."xxx".$a[0])){
	    push @a;
	}
    }
    elsif(/write/){
	print "@a$_";
	@a=();
    }
