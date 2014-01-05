#!/usr/bin/perl -w
########################################################################
####  This program is free software; you can redistribute it and/or modify
####  it under the terms of the GNU General Public License as published by
####  the Free Software Foundation; either version 3 of the License, or
####  (at your option) any later version.
####
####  This program is distributed in the hope that it will be useful,
####  but WITHOUT ANY WARRANTY; without even the implied warranty of
####  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
####  GNU General Public License for more details.
####
####  You should have received a copy of the GNU General Public License
####  along with this program.  If not, see <http://www.gnu.org/licenses/>.
####
####  If you don't understand what Free Software is, please read (or reread)
####  this page: http://www.gnu.org/philosophy/free-sw.html
########################################################################
use strict;
use Getopt::Std;
my %opts;
getopts('hvfcmdk:', \%opts);
my $missing=$opts{m}||undef;
my $column=$opts{k}||undef;
my $common=$opts{c}||undef;
my $verbose=$opts{v}||undef;
my $fast=$opts{f}||undef;
my $dupes=$opts{d}||undef;
$missing=1 unless $common || $dupes;;
&usage() unless $ARGV[1];
&usage() if $opts{h};
my (%found,%k,%fields);
if ($column) {
    die("The -k option only works in fast (-f) mode\n") unless $fast;
    $column--; ## So I don't need to count from 0
}

open(F1,"$ARGV[0]")||die("Cannot open $ARGV[0]: $!\n");
while(<F1>){
    chomp;
    if ($fast){ 
	my @a=split(/\s+/,$_);
	$k{$a[0]}++;	
        $found{$a[0]}++;
    }
    else {
	$k{$_}++;	
        $found{$_}++;
    }
}
close(F1);
my $n=0;
open(F2,"$ARGV[1]")||die("Cannot open $ARGV[1]: $!\n");
my $size=0;
if($verbose){
    while(<F2>){
	$size++;
    }
}
close(F2);
open(F2,"$ARGV[1]")||die("Cannot open $ARGV[1]: $!\n");

while(<F2>){
    next if /^\s+$/;
    $n++;
    chomp;
    print STDERR "." if $verbose && $n % 10==0;
    print STDERR "[$n of $size lines]\n" if $verbose && $n % 800==0;
    if($fast){
	my @a=split(/\s+/,$_);
	$k{$a[0]}++ if defined($k{$a[0]});
	$fields{$a[0]}=\@a if $column;
    }
    else{
	my @keys=keys(%k);
	foreach my $key(keys(%found)){
	    if (/$key/){
		$k{$key}++ ;
		$found{$key}=undef unless $dupes;
	    }
	}
    }
}
close(F2);
print STDERR "[$n of $size lines]\n" if $verbose;
#$missing && do map{print "$_ : $k{$_}\n" }keys(%k);
if ($column) {
    $missing && do map{my @a=@{$fields{$_}}; print "$a[$column]\n" unless $k{$_}>1}keys(%k);
    $common &&  do map{my @a=@{$fields{$_}}; print "$a[$column]\n" if $k{$_}>1}keys(%k);
    $dupes &&   do map{my @a=@{$fields{$_}}; print "$a[$column]\n" if $k{$_}>2}keys(%k);
}
else {
    $missing && do map{print "$_\n" unless $k{$_}>1}keys(%k);
    $common &&  do map{print "$_\n" if $k{$_}>1}keys(%k);
    $dupes &&   do map{print "$_\n" if $k{$_}>2}keys(%k);
}
sub usage{
    print STDERR <<EndOfHelp;

  USAGE: compare_lists.pl FILE1 FILE2

      This script will compare FILE1 and FILE2, searching for the 
      contents of FILE1 in FILE2 (and NOT vice versa). FILE one must 
      be one search pattern per line, the search pattern need only be 
      contained within one of the lines of FILE2.

    OPTIONS: 
      -c : Print patterns COMMON to both files
      -f : Search only the first characters of each line of FILE2
      for the search patern given in FILE1
      -d : Print duplicate entries     
      -m : Print patterns MISSING in FILE2 (default)
      -h : Print this help and exit
EndOfHelp
      exit(0);
}
