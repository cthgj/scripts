#!/usr/bin/perl -w
use strict;
use Getopt::Std;
use feature qw(switch say);
use File::Path;
our $verbose;
our $debug;
our $force;
sub v_say;
require "MY_SUBS.pl";
my %opts;
getopts('Bbh',\%opts);



my %hq;
if ($opts{B}) {
     %hq = %{get_hq_MIs('e')};
}
else {
    %hq = %{get_hq_MIs('b')};
}
while (<>) {
    my @mis=(/(MI:\d+)/);
    my $c=0;
    map{$c++ if defined($hq{$_})}@mis;
    print if $c>0;
}
