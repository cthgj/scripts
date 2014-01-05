#!/usr/bin/perl -w

use strict;

my $sec_ac="/home/terdon/research/data/sec_ac.txt";
my $acindex="/home/terdon/research/data/acindex.txt";

my (%accs,%ids,%names);

open(A1,"$sec_ac")|| die("Could not open $sec_ac :$!\n");
while(<A1>){
    next unless /^[\w\d]+\s+[\w\d]+$/;
    chomp;
    my @a=split(/\s+/,$_);
    ## $accs{old_acc}=new_acc
    $accs{$a[0]}=$a[1];
    ## $accs{new_acc}=new_acc
    $accs{$a[1]}=$a[1];
}
close(A1);
open(A2,"$acindex")|| die("Could not open $acindex :$!\n");
my $c=0;
while(<A2>){
    /______/ && do{$c=1; next};
    next if /^\s+$/;
    next unless $c==1;
    chomp;
    my @a=split(/\s+/,$_);
    ## ids{acc}=id, eg $id{A0A111}=PAFA_MYCTU
    $ids{$a[0]}=$a[1];
    ## $accs{PAFA_MYCTU}=A0A111
    $accs{$a[1]}=$a[0] unless defined($accs{$a[1]});
    
}
close(A2);
open(A3,"$ARGV[0]")|| die("need a list of accessions\n");
while(<A3>){
    chomp;
    #    print "$_\t
    ## These are the names I am interested in
    $names{$_}++;
}
close(A3);
foreach my $name (keys(%names)){
    unless((defined($ids{$name})) || defined($accs{$name})){
	print "$name\tmissing\n";
	next;
    }
    unless (defined($accs{$name})){
	if(defined($ids{$name})){
	    $accs{$name}=$accs{$ids{$name}};
	}
	else{$accs{$name}=$name;}
    }
    unless (defined($ids{$name})){
	$ids{$name}="no_id\n";
    }
    print "$name\t$ids{$name}\t$accs{$name}\n";
    
    


}


