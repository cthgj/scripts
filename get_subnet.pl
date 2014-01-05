#!/usr/bin/perl -w
###########################################################################
# This script will take a list of proteins, parse a class file		  #
# for their classes and return the subnetwork consisting of the 	  #
# nodes in each protein's classes, as well as each protein's interactors. #
###########################################################################
use strict;
use Getopt::Std;
my (%net,%net_prots,%opts,%wanted, %class_wanted,%classes);
getopts('n:c:h',\%opts);

#usage() unless $ARGV[1];
my $network_file=$opts{n}|| die "Need a network file (-n)\n";
my $prots=$ARGV[0] || die "Need a list of desired proteins as the 1st argument\n";
my $class_file=$opts{c}||die "Need a class file (-c)\n";

########################
# Collect wanted nodes #
########################
if (-e $prots) {
    open(PROTS,"$prots")||die("Cannot open list of nodes ($prots) : $!\n");
    while(<PROTS>){
	chomp;
	$wanted{$_}++;
    }
}
else{$wanted{$prots}++;}

## Collect wanted classes
open(CLAS,"$class_file")||die("Cannot open class file $class_file : $!\n");
my $class;
while (<CLAS>) {
    chomp;
    if(/^\[CLASS:\s*(\d+)/){$class=$1}
    my $w=0;
    if(/^PN\s*(.+)/){
	my @prots=split(/,\s+/,$1);
	map{$w++ if defined($wanted{$_})}@prots;
	next unless $w>0;
	map{
	    $class_wanted{$_}++;
	    $classes{$_}{$class}++;
	}@prots;
    }
}

load_network();
#################################
# Build the attribute files and #
# print the network.	        #
#################################

foreach my $pair (keys(%net)) {
   print "$pair\n";
}

open(A, ">$prots.atr");
print A "True_Class (class=String)\n";
foreach my $pp (keys(%net_prots)) {
    my @c=keys(%{$classes{$pp}});
    local $"="::";
   print A "$pp=(@c)\n";
}

sub load_network {
    open(my $A,"<", "$network_file")|| die("Cannot open $network_file : $!\n");
    while(<$A>){
	next if /^\d+$/;
	next if /^!/;
	next if /^\s*$/;
	chomp;
	my ($bait,$target)=split(/\s+/,$_);
	if (defined($wanted{$bait}) || defined($wanted{$target})) {
	    my $pair=join("\t",sort(split(/\s+/,$_)));
	    $net{$pair}++;
	    #print "$_\n";
	    $net_prots{$bait}++;
	    $net_prots{$target}++;
	}
	elsif (defined($class_wanted{$bait}) && 
	       defined($class_wanted{$target})) {
	    my $pair=join("\t",sort(split(/\s+/,$_)));
	    $net{$pair}++;
	    $net_prots{$bait}++;
	    $net_prots{$target}++;
	    #print "$_\n";
 
	}
    }
    close($A);
}
