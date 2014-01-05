#!/usr/bin/perl -w

use strict;

my %record;
foreach my $file (@ARGV){
    open(A,"$file");
    my $id="";
    my @lines=();
    while(<A>){
	push @lines,$_;
	if (/^ID\s+(.+?)\s/){
	    $id=$1;
	}
	elsif(/^\/\//){
	    $record{$id}=[@lines] unless defined($record{$id});
	    @lines=();
	}
	
    }
    close(A);
}

map{print @{$record{$_}}}keys(%record);
