#!/usr/bin/perl -w
# Author: cchapple 
# Created: 03 Apr 2012
use warnings;
use strict;
use Getopt::Std;
require "MY_SUBS.pl";

my (%ints, %opts);
getopts('vh:', \%opts);
our $verbose=$opts{v}||undef;
usage() if $opts{h};
usage() unless $#ARGV>0;
my (@dupes,@new);

## Open 1st net file
open(my $fh, "<", $ARGV[0]) or die "Could not open file $ARGV[0]: $!\n";
while(<$fh>){
    chomp;
    my @a=split(/\t/);
    my $p=join("\t", sort (@a));
    $ints{$p}++;
}
my $total_old_net=scalar keys(%ints);
## Open 2nd net file
open(my $fh1, "<", $ARGV[1]) or die "Could not open file $ARGV[1]: $!\n";
while (<$fh1>) {
    chomp;
    my @a=split(/\t/);
    my $p=join("\t", sort (@a));
    $ints{$p}++;
    if ($ints{$p}>1) {
        $verbose && do {
            push @dupes, $p;
        } ;
    }
    else {
        $verbose && do {
            push @new, $p;
        } ;
    }  
      
}
my %prots;
map{
    my @a=split(/\t/);
    $prots{$a[0]}=$prots{$a[1]}=1;
    print "$_\n"
}keys(%ints);

my $aa=scalar(@new);
my $bb=scalar(@dupes);
my $cc=$total_old_net-$bb;
my $dd=scalar keys(%ints);
my $ee=scalar keys(%prots);
v_say("## $aa New interactions in $ARGV[1]");
v_say("## $cc Interactions unique to $ARGV[0]");
v_say("## $bb Common interactions");
v_say("## $dd Total interactions");
v_say("## $ee Total proteins");

sub usage{
    my $us="<NET FILE1> <NET FILE2>";
    my $desc="This script will merge two networks removing any duplicate interactions. It assumes that both networks have the same ID type.";

    my %opts=(
        "usage" => $us,
        "desc" => $desc,
        "h" => "Print this help and exit"
    );
    print_help_exit(\%opts,0);
}

        
