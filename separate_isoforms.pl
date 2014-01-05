#!/usr/bin/perl 
use strict;

my (%found,%info);
while(<>){
    if (/^\!/){next};
    my ($a,$pid,$pname,$go,@rst)=split(/\t/);
    $rst[scalar(@rst)-1]=~/:(.+)/;
    my $iso=$1;
    push  @{$info{$pname}{ISOFO}}, $iso unless defined($found{$iso});
    $found{$iso}++;
    push  @{$info{$pname}{ANNOT}{$iso}},$go;
	
}
foreach my $pname (keys(%info)){
    print "########$pname##########\n";
	map{
	    my @a=sort(@{$info{$pname}{ANNOT}{$_}});
	    print "\n$_ : @a ";
	    }@{$info{$pname}{ISOFO}};
    print "\n";
}
