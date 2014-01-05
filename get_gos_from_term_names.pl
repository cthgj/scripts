#!/usr/bin/perl -w

my %k; 
while(<>){
    s/\s*GO/GO/;
    chomp;

    @a=split(/\s:\s/);
    
    $k{$a[1]}=$a[0];
}
open(A,"umbilical_Cons_D.txt.50maj.FM"); 
while(<A>){
    if(/^CA\s+(.+)/){
	@a=split(/\s+/,$1); 
	foreach my $j (@a){
	    $ka=$j;
	    $j=~s/_/ /g;
	    die("\$k{$j}\n") unless defined($k{$j});
	    s/$ka/$ka($k{$j})/g;
	}
    }
print;
}

