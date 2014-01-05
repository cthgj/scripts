#!/usr/bin/perl -w
use strict;

my (%accs,%nums,%interactors,%syns,%del);

my $gr=$ARGV[1];
my $map=$ARGV[0];

open(A,"$map")||die("cannot open $map : $!\n");
while(<A>){
    chomp; 
    my @b=split(/\s+/,$_); 
    map{$accs{$_}=$b[1]}@b;
    if ($b[1] eq "deleted"){$del{$b[0]}++ }
    else{
	$syns{$b[1]}{$b[0]}++;
#	$syns{$b[1]}{$b[2]}++;
#	push @{$syns{$b[1]}}, $b[0];
#	push @{$syns{$b[1]}}, $b[2] unless $b[0] eq $b[2];
	
    }
}

open(G,"$gr")||die("cannot open $gr : $!\n");
while(<G>){
    next if /^\d+$/; 
    chomp; 
    my @b=split(/\s+/,$_); 
    map{$nums{$_}++}@b; ## count interactions for each prot
    push @{$interactors{$b[0]}},$b[1];
    push @{$interactors{$b[1]}},$b[0];
}
my $c=0;
my %seen;
foreach my $name (keys(%accs)){
    my $total=0;
    next unless defined $nums{$name};
    defined $del{$name} ?  $seen{$name}++ : $seen{$accs{$name}}++;
#    print "see $accs{$name} $name: $seen{$accs{$name}}\n";
    
    $c++;
    if(defined $del{$name}) {print "$name\t$nums{$name}\tdeleted\n" unless $seen{$name}>1; next}
    else{
	next if $seen{$accs{$name}}>1;
	print " $accs{$name}\t";
	foreach my $syn (keys(%{$syns{$accs{$name}}})){
	    next unless defined($interactors{$syn});
	    print "\t$syn (" . scalar(@{$interactors{$syn}}) . ")";
	    $total+=scalar(@{$interactors{$syn}});
	}
    }
    
    print "\t$total\n";
}
