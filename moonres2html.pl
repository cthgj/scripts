#!/usr/bin/perl -w
###################################################################################################
# This script will take a results directory (containing the subdirectories Interactome and Both), #
# scan for <SPECIES>.<METHOD>.moon files and print an HTML table of the results			  #
###################################################################################################
use strict;
use Getopt::Std;
require "MY_SUBS.pl";
my (%opts,%kk);
getopts('hvM:',\%opts) || do { print "Try '" .  this_script_name() . " -h' for more information\n"; exit(1); };


my $results_dir=$ARGV[0]|| die("Need to specify a results directory ($ARGV[0]) : $!\n");
my $MODES=$opts{M}||"I,B";




######################################
# Get the modes we are interested in #
######################################
my @mm=split(/,/, $MODES);
my %modes;
map{$modes{$_}++}@mm;
@mm=sort  { "\L$a" cmp "\L$b" }keys(%modes);
print STDERR "m : @mm\n";
##########################################################
# Collect moonGO out files for all species. 		 #
##########################################################
opendir(RES, $results_dir) || die "Could not open results directory $results_dir for reading: $!\n";
my @files=grep{(/^(.+?)\.(\w)\.moon$/) && -f "$results_dir/$_" && defined($modes{$2})} readdir(RES);
close(RES);

my %sps;

foreach my $dir (".","interactome","both") {
    foreach my $file (@files) {
	my $num=get_count($dir,$file);
	my $real_dir=$dir;
	$real_dir="annotations" if $dir eq'.';
	$file=~/(.+?)\.(\w)\.moon$/;
	my ($sp,$m)=($1,$2);
	$sps{$sp}++;
	$kk{$real_dir}{$sp}{$m}=$num;
    }
}
my @species=keys(%sps);

print "<table  class=\"holder\"><tr class=\"heading\"><th>Annotations</th><th>Interactome</th><th>Both</th></tr><tr>";
foreach my $dir ("annotations","interactome","both"){
    print "<td><table class=\"res\">";
    print "<tr><td>&nbsp;</td>" if $dir eq 'annotations';
    map{print "<th>$_</th>";}@mm;
    print "</tr>";
    my $real_dir=$dir;
    $real_dir="." if $dir eq 'annotations';
    foreach my $s (@species) {
	print "<tr><th>$s</th>" if $dir eq 'annotations';
	map{print "<td><a href=\"./results/$real_dir/$s.$_.moon.html\">$kk{$dir}{$s}{$_}</a></td>"}@mm;
	print "</tr>\n";
    }
    print "</table></td>";
}
print "</table>";
#my @species_res=grep { (/^(.+?)\.(\w)\.moon$/) && -f "$results_dir/$_" && defined($modes{$2})} readdir(RES);
close(RES);


#echo "<table><tr><th>Annotations</th><th>Interactome</th><th>Both</th></tr>"; for t in . interactome both .; do for n in human worm fly yeast mouse; do echo "<tr>"; for m in a b B c i I P; do  echo "<td>`grep -v "#" $t/$n.$m.moon | cut -f 1 | src`</td>"; done; echo "</tr>"; done; done

sub get_count{
   my $dir=shift;
   my $file=shift;
   my $n=`grep -v "#" $dir/$file | cut -f 1 | sort | uniq | wc -l`;    
   chomp($n);
   return($n);
}
