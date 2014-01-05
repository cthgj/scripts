#!/usr/bin/perl -w
# translate_network_ids.pl --- 
# Author: cchapple 
# Created: 03 Apr 2012
use warnings;
use strict;
use Getopt::Std;
require "MY_SUBS.pl";

my (%ids, %opts);
getopts('hm:', \%opts);
usage() if $opts{h};

my $map=$opts{m}||die "Need a mapfile (-m)\n";


open(my $fh, "<", $map) or die "Could not open map file $map: $!\n";
while(<$fh>){
    chomp;
    my @a=split(/\t/);
    $ids{$a[0]}=$a[1];
}
while (<>) {
    chomp;
    my @a=split(/\t/);
    if ( defined($ids{$a[0]}) && defined($ids{$a[1]}) ) {
        my $p=join("\t", sort ($ids{$a[0]},$ids{$a[1]}));
        print "$p\n";
    }
    else {
	unless (defined($ids{$a[0]})){
	    print STDERR "No mapped ID for $a[0], SKIPPING\n";
	}
	unless (defined($ids{$a[1]})){	
	    print STDERR "No mapped ID for $a[1], SKIPPING\n";
	}

    }
    # next unless /729603/;
    #print "$ids{$a[0]}\t$ids{$a[1]}\n"
    
    
}


sub usage{
    my $us="-m MAP_FILE <NET FILE>";
    my $desc="This script will take a map and a network file and will replace ids in the network file with their counterparts from the map file. Any network IDs that are not present in the map file will be SKIPPED.";

    my %opts=(
        "usage" => $us,
        "desc" => $desc,
        "h" => "Print this help and exit",
        "m" => "Map file."
    );
    print_help_exit(\%opts,0);
}

        
