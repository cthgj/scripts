#!/usr/bin/perl -w
## USAGE : termstogoids.pl <terms file> <terms>

use strict;
use Getopt::Std;
my %opts;
getopts('gmt:',\%opts);
my $terms_file=$opts{t}||"$ENV{HOME}/research/GO/GO.terms_alt_ids";
my $return_terms=$opts{g}||undef;
my $print_missing=$opts{m}||undef;
# my $return_gos;
# my $return_terms==1 ? (my $return_gos =1) : (my $return_gos =0);

# print STDERR "$return_terms:$return_gos\n";die();


my (%terms_to_GOs,%GOs_to_terms);
die("USAGE : termstogoids.pl [-t <terms file>] <terms>\n\nBy default, it expects a list of term descriptions (eg system_process) and will return a list of GOs. For the reverse, use the -g flag\n") unless $ARGV[0];
my @terms;
if(-e $ARGV[0]){
    open(A,$ARGV[0]);
    while(<A>){
	chomp;
	s/\s+$//;
	s/\s/_/g;
	push @terms, $_;
    }
    close(A);
}
else{$terms[0]=$ARGV[0];}
open(T,"$terms_file")|| die("Cannot open $terms_file : $!\n");
while(<T>){
    next if /^\!/ || /^\s+$/; 
    chomp;
    s/^\s*//;
    my @a=split(/\t/);
    push my @gos, $a[0];
    ## If we have alternative IDs
    if($a[1]=~/./){
	map{push @gos,$_}split(/\s+/,$a[1]);
    }
    $a[2] =~s/\s/_/g;
    map{
	#print "xxxx \$terms_to_GOs{$a[2]}=$_\n\$GOs_to_terms{$_}=$a[2]\n";
	$terms_to_GOs{$a[2]}=$_;
	$GOs_to_terms{$_}=$a[2];
    }@gos;
}
close(T);

if ($print_missing) {
	map{print "$_\n" unless defined($GOs_to_terms{$_})}@terms;
}
else {
    if($return_terms){
	map{s/ /_/g;print "$_ : $GOs_to_terms{$_}\n";}@terms;
    }
    else{
	map{s/ /_/g;print "$_ : $terms_to_GOs{$_}\n";}@terms;
    }
}





