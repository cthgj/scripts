#!/usr/bin/perl -w

use strict;

my($subonto);
local $/ = "\n\n";
my $a=0;
while(<>){
    next if $.==1;
    print;die();
    if (/namespace:\s*(.+)/){
	$subonto=$1;
    }
    



}
