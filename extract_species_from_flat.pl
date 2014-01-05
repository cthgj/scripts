#!/usr/bin/perl
use strict;
use warnings;
use SWISS::Entry;
use Getopt::Std;
my %opts;
getopts('s:',\%opts);

my $species=$opts{s}||die();
# Change the line termination string so we read an entire entry at a time
local $/ = "\n//\n";


while(<>){
    next unless /OS\s+.*$species.*/i;
    print;
    
    # my $entry = SWISS::Entry->fromText($_);
    # if($entry->OSs->head->text =~ /$species/i){
    # 	print  $entry->toText;
    # }
    # else {
    # 	next;
    # }
}
