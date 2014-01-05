#!/usr/bin/perl -w


use Switch;
use strict;

my (%ancestors,%offspring);
my $geneology_file=$ARGV[0];

open (ANC,"$geneology_file")|| die("cannot open $geneology_file:$!\n");
while(<ANC>){
    chomp;
    my @terms=split(/\t/);
    my $child=shift(@terms);
    $ancestors{$child}=[@terms];
    map{$offspring{$_}{$child}++}@terms;
    
}
close(ANC);

map{
    &get_direct_papas($_);
}keys(%ancestors);

sub get_direct_papas{
    my $term=shift;
    my @papas=@{$ancestors{$term}};
    my %count;
    print "tt : $term : \n";
    for(my $n=0; $n<=$#papas; $n++){
	$count{$papas[$n]}++;
	for(my $k=0; $k<=$#papas; $k++){
	    #print "$papas[$k] $papas[$n] $offspring{$papas[$k]}{$papas[$n]}\n";
	    if(defined($offspring{$papas[$n]}{$papas[$k]})){
		$count{$papas[$n]}++;
		print "xx $papas[$n] $papas[$k]\n";
	    }
	}
    }
    map{print "$_ : $count{$_}\n" if $count{$_}==1  }keys(%count);
    die();
	
    # for(my $n=0; $n<$#papas; $n++){
    # 	for(my $k=$n; $k<=$#papas; $k++){
    # 	    if(defined($offspring{$papas[$n]}{$papas[$k]})){
    # 		die("$papas[$n] : $papas[$k]\n")
    # 	    }
    # 	}
    # }
    
}
sub get_direct_papas_papas{

}
