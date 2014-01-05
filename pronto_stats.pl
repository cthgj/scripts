#!/usr/bin/env perl

my %onto;
## Open GO.terms_alt_ids file
open(A,"$ARGV[0]"); 


while(<A>){
    next if /^!/; 
    chomp; 
    @a=split(/\t/); 
    $onto{$a[0]}=$a[$#a];
}
close(A);
my %counts;
open(B,"$ARGV[1]");
while (<B>) {
    print STDERR "$.\r" if $. % 100 == 0;
    next if /^#/;
    @a=split(/\t/);
    my $oo=join("", sort ($onto{$a[1]},$onto{$a[2]}));
    $counts{$oo}++;
}

foreach my $oo (keys(%counts)) {
    print "$oo\t$counts{$oo}\n";
}
print STDERR "Done";
