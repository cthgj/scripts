#!/usr/bin/perl -w

use strict;
use Getopt::Std;
my %opts;

getopts('g',\%opts);
my $fname=$ARGV[0] . "_grouped";

my (%assoc,%k); 
while(<>){
    /^(.+?)\s(.+)/; 
    push @{$k{$1}}, $2
}
 if ($opts{g}){
     open (A,">$fname");
     map{print A "$_\t@{$k{$_}}\n"}keys(%k);
}


foreach my $gene_name (keys %k){
    my @terms=sort(@{$k{$gene_name}});
 #   print"TERNMS : @terms\n";
#    next if (scalar(@terms==1));
    for(my $n=0; $n<scalar(@terms)-1; $n++){
	for(my $k=0; $k<scalar(@terms)-1; $k++){
	    next if $terms[$n] eq $terms[$k];
	    
	    my $aa=$terms[$n] . "--" . $terms[$k];
	    my $bb=$terms[$k] . "--" . $terms[$n];
	    my $gterm;
	    
	    if(defined($assoc{$aa}) && (!defined($assoc{$bb})) ){
		$gterm=$aa;
	    }
	    elsif(!defined($assoc{$aa}) && (!defined($assoc{$bb})) ){
		$gterm=$aa;
	    }
	    elsif(!defined($assoc{$aa}) && (defined($assoc{$bb})) ){
		$gterm=$bb;
	    }
	    else{die("kkkkkn $aa $bb");}

	    die("crap\n") if(defined($assoc{$bb}) && defined($assoc{$aa}));
	
	    my $kkname=$gene_name . "--" . $gterm;
	    
	    $assoc{$gterm}{NUM}++ unless defined($assoc{$gterm}{GENES}{$kkname});
	    push @{$assoc{$gterm}{NAMES}},$gene_name unless defined($assoc{$gterm}{GENES}{$kkname});
	    $assoc{$gterm}{GENES}{$kkname}++;
	}
    }
}

map{print "$_ : $assoc{$_}{NUM} (@{$assoc{$_}{NAMES}})\n"}keys %assoc;
