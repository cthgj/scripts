#!/usr/bin/env perl

## This script expects an input file like

# CC	171405
# CF	687960
# CP	1678365
# FF	692076
# FP	3373944
# PP	4117015

## And will print out an html table

my $spc;
my %ontos;
my %sort_hash;

my $header=shift;
foreach my $file (@ARGV) {
    my $fname=`basename $file`;
    if ($fname=~/([^.]+).*\.stats/) {
	$spc=$1;
    }
    else {
	print STDERR "Unkown species for : $file\n";
    }
    open(my $fh, '<', $file)||die("Could not open file $file: $!\n");
    while (<$fh>) {
	chomp;
	@a=split(/\s+/);
	$k{$spc}{$a[0]}=$a[1];
	$ontos{$a[0]}++;
    }
}
print "<div class=\"statsTable\"><table class=\"stats\"><thead><tr><th colspan=7>$header</th></tr><tr><th class=\"noborder\">&nbsp;</th>";
my @o=sort keys(%ontos);
foreach (@o) {
    print "<th>$_</th>";
}
print "</tr></thead><tbody>";

foreach my $species (keys(%k)) {
    $sort_hash{$species}=$k{$species}{'PP'}
}


foreach my $species (sort {$sort_hash{$foo} gt $sort_hashk{$bar}} keys(%sort_hash)) {
    print "<tr><td class=\"speciesTD\">$species</td>";
    foreach my $oo(@o) {
	print "<td>$k{$species}{$oo}</td>";
    }
    print "</tr>" ;   
}
print "</tbody></table></div>\n";
