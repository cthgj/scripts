#!/usr/bin/perl -w
use strict;
use Getopt::Std;


my %opts;
getopts('hvo:g:G:',\%opts) || do { print "Invalid option, try '$0 -h' for more information\n"; exit(1); };

my %go_pairs;
my $geneology_file=$opts{g}||"/home/terdon/research/GO/biological_process.geneology";
my $gaf_file=$opts{G}||"./data/gene_association.human";
my $verbose=$opts{v}||undef;
my $subonto=$opts{o}||"P";



open(A, $ARGV[0]);
while(<A>){
    chomp;
    next if /GO:0008150/;
    my @a=split(/\s/);
    $go_pairs{$a[0]}{$a[1]}=$a[2];
}
my %ancestors=&load_ancestors($geneology_file);
my %proteins=&parse_gaf_file($gaf_file);

##Now, count

my $c=0;
my $b=scalar(keys(%proteins));
foreach my $prot (keys(%proteins)){
	printf STDERR ("$c of $b\r") if $verbose;
    foreach my $go1 (keys(%go_pairs)){
	foreach my $go2 (keys(%{$go_pairs{$go1}})){
	    if((defined($proteins{$prot}{GOs}{$go1})) &&
	       (defined($proteins{$prot}{GOs}{$go2}))){
		print "$prot\t$go1\t$go2\t$go_pairs{$go1}{$go2}\n"; 
	    }
	}
    }
    $c++;
}




sub load_ancestors{
    print STDERR "Loading ancestors...\n" if $verbose;
    open (ANC,"$_[0]")|| die("cannot open $_[0]:$!\n");
    my (%ancestors,%offspring);
    while(<ANC>){
	chomp;
	my @terms=split(/\t/);
	my %seen;
	@terms = grep { ! $seen{$_} ++ } @terms;
	my $child=shift(@terms);
	$ancestors{$child}=[@terms];
	map{$offspring{$_}{$child}++}@terms;
    }
    close(ANC);
    return(%ancestors)
}
sub parse_gaf_file{    
    print STDERR "Parsing GAF file...\n" if $verbose;
    open(A,"$_[0]")|| die("Cannot open $_[0]:$!\n");
    my %proteins;
    while(<A>){
	next if /^!/;
	chomp;
	my @tmparray=split(/\t/,$_);
	next unless $tmparray[8] eq $subonto;
	next if $tmparray[3] !~/^$/;
	my $name=$tmparray[1];
	$proteins{$name}{GOs}{$tmparray[4]}="DIR";
	map{$proteins{$name}{GOs}{$_}="ANC"}@{$ancestors{$tmparray[4]}};
	
    }
    close(A);
    return(%proteins);
}
