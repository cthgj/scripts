#!/usr/bin/perl -w

use strict;
use Getopt::Std;
use Term::ANSIColor; 

my %opts;
getopts('chv', \%opts);
my $verbose=$opts{v}||undef;
my $flat=$opts{f}||undef;
my $color=$opts{c}||undef;
my $sec_ac="/home/terdon/research/data/sec_ac.txt";
my $acindex="/home/terdon/research/data/acindex.txt";
my $del="/home/terdon/research/data/deleted";
my $id;
my (%primary_accs,%del,%name,%accs,%ids,%names,%acc);

open(A1,"$sec_ac")|| die("Could not open $sec_ac :$!\n");
while(<A1>){
    next unless /^[\w\d]+\s+[\w\d]+$/;
    chomp;
    my @a=split(/\s+/,$_);
    ## $accs{old_acc}=new_acc
    push @{$primary_accs{$a[0]}}, $a[1];
    ## $accs{new_acc}=new_acc
    push @{$primary_accs{$a[1]}}, $a[1];
}
close(A1);
open(A2,"$acindex")|| die("Could not open $acindex :$!\n");
my $c=0;
while(<A2>){
    /______/ && do{$c=1; next};
    next if /^\s+$/;
    next unless $c==1;
    chomp;
    my @a=split(/\s\s/,$_);
    my @b=split(/,/,$a[1]);
    # if (defined($accs{$a[0]})){
    # 	print "$a[0] : @{$accs{$a[0]}}\n";
    # }    

## ids{acc}=id, eg $ids{A0A111}=PAFA_MYCTU
    push @{$ids{$a[1]}},$a[0];
    my @la=split(/,/,$a[1]);
## name{ID}=acc, eg $name{PAFA_MYCTU}=A0A111
    push @{$name{$a[0]}},$a[1];
    ## $accs{PAFA_MYCTU}=A0A111
    # if(defined($primary_accs{$a[0]})){
    # 	map{print "$a[0] : $_,"}@{$primary_accs{$a[0]}};
    # 	print "\n";
    # 	$accs{$a[0]}=0;
    # }
    # push @{$accs{$a[1]}},$a[0];
    
}
close(A2);
open(A3,"$del")|| die("Could not open $del :$!\n");
$c=0;
while(<A3>){
    /______/ && do{$c=1; next};
    next if /^\s+$/;
    next unless $c==1;
    chomp;
    $del{$_}++;
    
}
close(A3);
## Read flat files
if($flat){
    if($flat =~ /\.gz$/) 
    {
	open(F,"zcat $flat |")||die("Need a flat file : $!\n");
    } 
    else{
	open(F,"$flat")||die("Need a flat file :$!\n");
    }
    my @lines;
    my $o;
    while(<F>){
	my @acc_list=();
	$.==1 && do{
	    die("This is not a flat file\n") unless /^ID/;
	};
	chomp;
	push @lines,$_;
	if(/^ID\s+(.+?)\s/){
	    $id=$1;
	}
	elsif(/^AC\s+(.+)/){
	    @acc_list=split(/;\s*/,$1);
	    map{ push @{$acc{$_}},$acc_list[0] }@acc_list;
	    $acc{$id}=[@acc_list];
	    $ids{$acc_list[0]}=$id;
	    
	}
	
	
    }

}
close(F);
open(A3,"$ARGV[0]")|| die("need a list of accessions\n");
while(<A3>){
    chomp;
    #    print "$_\t
    ## These are the names and accessions I am interested in
    my ($id,$acc)=split(/\t/,$_);
    $names{$_}++;
    &update_name($id,$acc);
    
}
close(A3);


sub update_name{
    my $id=shift;
    my $acc=shift;
    my $prob=0;
    my ($newname,$newacc)=($id,$acc);
    my $a;
    
    ## ids{acc}=id, eg $id{A0A111}=PAFA_MYCTU
    ## if we have an id for this accession in $acindex
    if(defined($ids{$id})){
	## If we have found >1 id for this accession
	if(${$ids{$id}}[1]){
	    ## foreach accession 
	    foreach my $accession (@{$ids{$id}}){
		my $last_acc=${$primary_accs{$accession}}[0] || die("No primary accession for $accession\n");
		## Check that the 1ary acc for this acc is always the same
		foreach my $p_acc (@{$primary_accs{$accession}}){
		  $p_acc ne $last_acc && do {
		      $prob=1;  
		      print "PROBLEM: $id : $accession : $p_acc\n";
		  };
		  $last_acc=$p_acc;
		}
	    }
	    ## Unless there were multiple 1ary accessions, acc is now
	    ## the accession read from $seq_ac
	    if($prob==0){
		$newacc=${$primary_accs{${$ids{$id}}[0]}}[0];
	    }
	}
	## If we only have 1 id for this accession
	## then that is the new accession
	else{
	    $newacc=${$ids{$id}}[0];
	}
    }
    
    else{
	## If there is no id for this accession in $acindex
	## but there is a name for this accession
	if(defined(${$name{$acc}}[0])){#die("$id : $acc\n");
	    $newname=${$name{$acc}}[0];
	}
	else{
	   # die("No info for $id:$acc\n");
	}
## If there is a primary accession for this $acc
	if(defined(${$primary_accs{$acc}}[0]) ){
	    $newacc=${$primary_accs{$acc}}[0];
	}
    }
    if(defined($del{$acc})){
	$newname=$newacc="deleted";
    }
    $a="$id\t$newname\t$newacc\n";
    if($id ne $newname){
	$color && do {$a=~s/(\t$newname\t)/color("bold blue").$1.color("reset")/ige	unless $newname=~/deleted/}
    }
   if($acc ne $newacc){
      $color && do { $a=~s/(\t$newacc)/color("bold red").$1.color("reset")/ige 	unless $newacc=~/deleted/;}
    }
    print "$a";
 }
