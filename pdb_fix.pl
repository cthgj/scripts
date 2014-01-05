#!/usr/bin/perl -w
use strict;
use Getopt::Long qw(:config bundling);
use feature qw(switch say);
use File::Path;
our $verbose;
our $debug;
our $force;
sub v_say;
require "MY_SUBS.pl";

my %rogids;
open(R, $ARGV[0]);
while (<R>) {
    chomp;
    my @a=(/rogid:(.+?)\|/g);
    map{$rogids{$_}++;}@a;

}
my $fh=check_file($ARGV[1],"r");
my %map;
while (<$fh>) {
    chomp;
    my @a=split(/\t/);
    my $l=$_;
    print "$a[4] : $_\n" if defined($rogids{$a[4]});
    $map{$a[4]}=$a[1] if defined($rogids{$a[4]});
}
print STDERR "Done\n";


open(B, $ARGV[2]);
my %map2;
while(<B>){
    chomp;
    my @a=split(/\t/);

}


