#!/usr/bin/perl -w

use strict;
my @clusnames;
my @mclus;

my $pat = 'std';


    open (CLUSTERS, "cluster_morethan2_$pat")|| die "cannot open cluster_morethan2_$pat: $!";
    
    my @clusses= <CLUSTERS>;                  ### get all cluster names

    foreach my $a (@clusses)
    {
	my @e = split /\s+/go, $a;            ### @e = list of INDIVIDUAL clusters
	map { push @clusnames, $1 if $_ =~ /(TOG\d{6})/go;} @e;   ### @clusnames = clusters_found > 2
    }
close(CLUSTERS);

 





   
CLUS:    foreach my $clus (@clusnames){
    print STDOUT "=====================$clus=================================\n";
	
    open(FILE,"TOG_files/$clus\_$pat\.tbl")|| next;
    my @lines = <FILE>;
    $lines[0] =~ /^(.*?)_/;
    my $prevspec = $1;
    my @nms;
    close(FILE);

    open(FILE,"TOG_files/$clus\_$pat\.tbl")|| die "cannot open TOGfile: $!";
   
    print STDOUT "\$prevspec is now $prevspec\n";
 
    while(<FILE>){
	print STDOUT "----------$.--------------\n";
	next if $_ =~ /$prevspec/;
	print STDOUT "###\n";
	push @nms, $clus;
	next CLUS;}
	



    
	    
    open(MCLUS, ">mclus");
    map{print MCLUS "$_\n"} @mclus;
    close(MCLUS);
	
    close(FILE);
    
    
}
