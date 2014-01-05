#!/usr/bin/perl -w

use strict;

&usage unless $ARGV[2];
my @terms;
my (%new_terms,%old_terms);

open(A,$ARGV[0]);
while(<A>){
    next if /^!/;
    chomp;
    my @aa=split(/\t/);
    my @a=@aa;
    my @GOs=($aa[0]);
    if(/obs$/){pop(@a);}
    if($aa[1] =~/GO:/){
	push @GOs,split(/\s+/,$aa[1]);
    }
    $a[$#a-1]=~s/ /_/g;
    ## $old_terms{GOs}{GO:0008150} : biological_process_unknown
    ## $old_terms{TERMS}{biological_process_unknown} : GO:0008150
    $old_terms{TERMS}{$a[$#a-1]}=$a[0];
    map{$old_terms{GOs}{$_}=$a[$#a-1];}@GOs;
}
close(A);

open(A,$ARGV[1]);
while(<A>){
    next if /^!/;
    chomp;
    my @aa=split(/\t/);
    my @a=@aa;
    my @GOs=($aa[0]);
    if(/obs$/){pop(@a);}
    if($aa[1] =~/GO:/){
	push @GOs,split(/\s+/,$aa[1]);
    }
    $a[$#a-1]=~s/ /_/g;
    $new_terms{TERMS}{$a[$#a-1]}=$a[0];
    map{$new_terms{GOs}{$_}=$a[$#a-1];}@GOs;
}
close(A);

open(A,$ARGV[2]);
while(<A>){
    if(/^CA\s+(.+)\s*$/){
	my @terms=split(/\s+/,$1);
	foreach my $term (@terms){
	    my $go=$old_terms{TERMS}{$term}||die("BAD term : -$term-\n$old_terms{TERMS}{$term} : $new_terms{TERMS}{$term}");
	    if($new_terms{GOs}{$go} ne $old_terms{GOs}{$go}){
#		print STDERR "OLD ($go) : $old_terms{GOs}{$go}\n";
#		print STDERR "NEW ($go) : $new_terms{GOs}{$go}\n";
		s/$term/$new_terms{GOs}{$go}/;
	    }
	}
    }
    elsif(/^CM\s+(.+)\s+on\s*\d/){
	my @terms=split(/\s+/,$1);
	foreach my $term (@terms){
	    my $go=$old_terms{TERMS}{$term}||die("BAD term : -$term-\n$old_terms{TERMS}{$term} : $new_terms{TERMS}{$term}");
	    if($new_terms{GOs}{$go} ne $old_terms{GOs}{$go}){
		s/$term/$new_terms{GOs}{$go}/;
	    }
	}
    }
	print;
    
}


# map{print "$_ : $old_terms{GOs}{$_}\n"}keys(%{$old_terms{GOs}});
# print "\n-------------------\n";
# map{print "$_ : $old_terms{TERMS}{$_}\n"}keys(%{$old_terms{TERMS}});
# die();

sub usage{
    print STDERR "USAGE : $0 <old terms file> <new terms file> <class file>\n";
    exit();
}
