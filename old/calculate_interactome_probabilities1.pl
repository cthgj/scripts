#!/usr/bin/perl 
## Get the numbers necessary to calculate the probabilities of association of each GO pair
## based on protein interactions and NOT annotations

use strict;
use Getopt::Std;
use IO::File;
# use GO::Basic;
# use GO::Parser;
# use GO::Model::Association;
my %opts;
getopts('a:g:n:s:o:uh',\%opts) || do { print "Invalid option, try '$0 -h' for more information\n"; exit(1); };
my $gaf_file=$opts{a}||die("Need a gene annotation (GAF) file, -a\n");
my $gen_file=$opts{g}||die("Need a geneaology file, -g\n ");
my $net_file=$opts{n}||die("Need a network file, -n\n ");
my $synfile=$opts{s}||die("Need a synonyms file, -s\n ");
my $subonto=$opts{o}||'P';
my %synonyms;
my %papas=&load_geneology($gen_file);
## parse GAF
my %proteins=&load_network($net_file); 
my $suffix='HUMAN';
my (%synonyms,%gos, %count,%gocount);
my $under_only=$opts{u}||undef;


open(GAF,"$gaf_file")||die("cannot open $gaf_file: $!\n");
while(<GAF>){
    next if /^!/;
    print STDERR "Gaf : $.\r";
    chomp;
    my @tmparray=split(/\t/);
    unless($subonto eq 'all'){next unless $tmparray[8] eq $subonto;}
    next unless $tmparray[3]=~/^$/;  ## skip the NOT annotation
    my $name=&get_name($tmparray[1],"NET",1);
    next unless defined($proteins{$name});
    $gos{$name}{$tmparray[4]}++ if defined($proteins{$name});## $tmparray[1]== name, $tmparray[4] == GO
}
close(GAF);
my (%seen,%go_count,%good_genes);


## Get ancestors
foreach my $prot (keys(%proteins)){
    foreach my $ggo (keys(%{$gos{$prot}})){
	$go_count{$ggo}++ unless $seen{$ggo}{$prot};
	$seen{$ggo}{$prot}++;
	map{
	    $go_count{$_}++ unless $seen{$_}{$prot};
	    $seen{$_}{$prot}++;
	}@{$papas{$ggo}}
    }
}
my %pairs;
my $pcount=scalar(keys(%proteins));
print STDERR "\nPROT : $pcount\n";
my @GOs=keys(%seen);
for (my $n=0; $n<=$#GOs;$n++){
    print STDERR "$n of $#GOs\r";
    next if $GOs[$n] eq 'GO:0008150' ;
    next if $GOs[$n] eq 'GO:0005575' ;
    for (my $k=$n+1; $k<scalar(@GOs);$k++){
	next if $GOs[$k] eq 'GO:0008150' ;
	next if $GOs[$k] eq 'GO:0005575' ;
	my @b=sort {$b lt $a} ($GOs[$n],$GOs[$k]);
	my $pair = join("_",@b);
	foreach my $prot (keys(%{$seen{$GOs[$n]}})){
	    $pairs{$pair}++ if defined($seen{$GOs[$k]}{$prot})
	}
    }
}

foreach my $pair (keys(%pairs)){
    my @gos=split(/_/,$pair);
    my $o=$pairs{$pair}*$pcount/$go_count{$gos[0]}/$go_count{$gos[1]};
    if($under_only) {
	print "$pair\t$pcount\t$go_count{$gos[0]}\t$go_count{$gos[1]}\t$pairs{$pair}\t$o\n" if $o < 1;
    }
    else{
	print "$pair\t$pcount\t$go_count{$gos[0]}\t$go_count{$gos[1]}\t$pairs{$pair}\t$o\n";
    }
}


sub load_geneology{
    my $file=shift;
    my %ancestors;
    open(A,"$file")||die("Cannot open $file : $!\n");
    while(<A>){
	chomp;
	my @terms=split(/\t/);
	my %seen;
	@terms = grep { ! $seen{$_} ++ } @terms;
	my $child=shift(@terms);
	$ancestors{$child}=[@terms];
    }
    print SDTERR "GEN : $.\r";

    return(%ancestors);
}
###############################################################
sub load_network {
    my %proteins;
    my $network_file=shift;
    my $A = IO::File->new("< $network_file")|| die("Cannot open $network_file : $!\n");
    print STDERR "Loading Network...";
    while(<$A>){
	next if /^\d+$/;
	next if /^!/;
	next if /^\s*$/;
	chomp;
	my ($bait,$target)=split(/\t/,$_);
	$bait=&get_name($bait,"NET",2);
	$target=&get_name($target,"NET",3);
	## %proteins will hold all the info on all the 
	## network's proteins
	$proteins{$bait}{INTERACTORS}{$target}++;
	$proteins{$target}{INTERACTORS}{$bait}++;
    }
    print STDERR "Done\n";
    close(A);
    return(%proteins)
}
############################################################
sub get_name{
    my $name=$_[0];
    my $old_name=$name;
    my $want=$_[1]||"null";
    my $call=$_[2]||undef;
    
    if ($synonyms{LOADED}==0){
	if(-e $synfile){
	    open(S,$synfile);
	    while(<S>){
		chomp;
		my @a=split(/\t/);
		map{s/\s//g}@a;
		$synonyms{NET}{$a[1]}=$a[0];
		$synonyms{GAF}{$a[1]}=$a[1];
		$synonyms{GAF}{$a[0]}=$a[1];
		$synonyms{NET}{$a[0]}=$a[0];
	    }
	    close(S);
 	$synonyms{LOADED}=1;

	}
    }
    $want eq 'NET' ? 
	return($synonyms{NET}{$name}) :
	return ($synonyms{GAF}{$name});

   
}
