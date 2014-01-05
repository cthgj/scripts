#!/usr/bin/perl -w

use strict;
use Getopt::Std;
my %opts;
getopts('a:n:g:',\%opts) || do { print "Invalid option\n"; exit(1); };
my $genealogy_file=$opts{g} || die("Need a genealogy file, -g (eg ~/research/GO/all_three_Feb11-2013.genealogy: $!\n");
my $alt_file=$opts{a}||"/home/terdon/research/GO/GO.terms_alt_ids";

my %obs;
open(ALT,"$alt_file")||die("Could not open $alt_file: $!\n");
while(<ALT>){
    chomp;
    next unless /^GO:/;
    die("crap: $_\n") unless /^(GO:.+?)\s/;
    my $t=$1;
    $obs{$t}++ if /\tobs$/;
    
}

my ($a,$b,$longest_branch,$term_num)=&load_ancestors($genealogy_file);
#my (%ancestors,%offspring)=(%$a,%$b);
my %ancestors=%$a;
my %offspring=%$b;

my @terms;
if($opts{n}){
    @terms=split(/\s/,$opts{n});
}
elsif ($ARGV[1]) {
    open(A,"$ARGV[1]")||die("Need a list of terms either as a space separated list (-n) or a second argument\n");
    while(<A>){
	chomp;
	next unless /^..:/;
	push @terms, $_;
    }

}
## Read from STDIN
else {
  while(<>){
	chomp;
	next unless /^..:/;
	push @terms, $_;
    }  
}

foreach my $term (@terms){
    next if defined($obs{$term});
    my $spec;
    if (defined($ancestors{$term})){
	$spec=&calculate_specificity($term);
	print "$term\t$spec\n";
    }
    else{
	print "$term\t0\n";
    }
}


sub load_ancestors{
    my $genealogy_file=shift;
    my (%ancestors,%offspring);
    my $longest_branch=0;
    my $term_num=0;
    open (ANC,"$genealogy_file")|| die("cannot open $genealogy_file:$!\n");
    while(<ANC>){
	chomp;
	$term_num++;
	my @terms=split(/\t/);
	my $child=shift(@terms);
	$ancestors{$child}=\@terms;
	scalar(@terms)>$longest_branch ? $longest_branch=scalar(@terms) : $longest_branch=$longest_branch;
	map{$offspring{$_}{$child}++}@terms;
    }
    close(ANC);
    return(\%ancestors,\%offspring,$longest_branch,$term_num);
}

sub calculate_specificity{
    my $term=shift;
    my $anc=@{$ancestors{$term}}||1;
    my $kids=keys(%{$offspring{$term}})||1;
#    print STDOUT "log($kids/$term_num*$anc)/$term_num*$longest_branch\n";
    my $spec= -(log($kids/($term_num*$anc))/log($term_num*$longest_branch));
    return $spec;

}
