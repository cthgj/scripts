#!/usr/bin/perl -w

use strict;
use Getopt::Std;
getopts('a',\%opts) || do { print "Invalid option\n"; exit(1); };
my $network_file=$ARGV[0]|| die("Need a network file\n");
my $gaf_annotations_file=$ARGV[1]|| die("Need a GAF file\n");
my $geneology_file=$ARGV[2]||die("Need a geneology file\n");;



