#!/usr/bin/perl -w

## This file will take two term files, the "wrong" and "new" and a 
## class anotation file and correct any of the "wrong" terms in the 
## class file with the "new" terms

## USAGE fix_terms.pl -w <WRONG> -n <NEW> <ANNOTATED_CLASS_FILE>

use strict;
use Getopt::Std;

my (%wrong_term,%new_term,%wrong_GO,%new_GO,%opts);

getopts('w:n:',\%opts);

open(WRONG,$opts{w});
while(<WRONG>){
    next if /!/;
    /^(.+?)\t(.+?)[\t\n]/;
    my ($a,$b)=($1,$2);
    $b=~s/ /_/g;
    $wrong_GO{$b}=$a;
    $wrong_term{$a}=$b;    
    print STDERR "0bb : $b,$a\n" if /0071844/;
}
close(WRONG);

open(NEW,$opts{n});
while(<NEW>){
    next if /!/;
    chomp;
    my @a=split(/\t/);
    my @b=split(/\s/, $a[$#a-2]);
#    my @b=grep(/GO:\d+/, @a);
#    print STDERR "1bb : $a[$#a-2] : @a : @b : $b[1]\n$_" if /0071844/;
#    print "a : @b\n$a[$#a-2]\n"; die();
    $a[$#a-1]=~s/ /_/g;
    $new_term{$a[0]}=$a[$#a-1];   
    $new_GO{$a[$#a-1]}=$a[0];
    map{
	$a[$#a-1]=~s/ /_/g;
#	if(/0071844/){print STDERR "bb $_ ::\$new_GO{$a[$#a-1]}=$b[0] \$new_term{$_}=$a[$#a-1]\n";}
	$new_GO{$a[$#a-1]}=$a[0];
	$new_term{$_}=$a[$#a-1];    
    }@b;
#    /^(.+)\t(.+?)$/;
#    print "xx $1 :: $2\n";
    # my ($a,$b)=($1,$2);
    # $b=~s/ /_/g;
    # $new_GO{$b}=$a;
    # $new_term{$a}=$b;    
}
close(NEW);

while(<>){
    s/\t\n/\n/; 
    s/\t\s/\t/g;
    if((/^CA\s+(.+)/) ||(/^CM\t(.+?)\s+\d/)){	
	#print STDERR "======================================================\na1 $_"; 
	my @processes=split(/\s+/,$1);
	foreach my $a (@processes){
	    #print STDERR "a1AA  a:$a  w:$wrong_GO{$a}  n:$new_GO{$a} : \$new_term{$wrong_GO{$a}}=$new_term{$wrong_GO{$a}} \n$_\n" ;
	    unless(defined($new_GO{$a})){
		s/$a/$new_term{$wrong_GO{$a}}/;
	    }
	}
    }
	print;

}
