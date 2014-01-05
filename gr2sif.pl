#!/usr/bin/perl -w
use strict;


while(<>){
    my @a=split(/\s+/,$_);
    print "$a[0]\tpp\t$a[1]\n";
}
