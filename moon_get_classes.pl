#!/usr/bin/perl -w

use strict;
use Getopt::Std;
my (%annots,%opts,%prots,%classes,%want);
getopts("c:ha", \%opts);
my $wanted_classes=$opts{c}||undef;

if($wanted_classes){
    open(A,"$wanted_classes")||die("Cannot open wanted classes file $wanted_classes: $!\n");
    while(<A>){
	chomp;
	$want{$_}++;
    }
}

if(-e $ARGV[0]){
    open(PROTS,"$ARGV[0]")||die("Cannot open list of names (ARGV[0]) : $!\n");
    while(<PROTS>){
	chomp;
	$prots{$_}++;
    }
    close(PROTS);
}
else{$prots{$ARGV[0]}++;}


open(CLASS,"$ARGV[1]")||die("Cannot open class file (ARGV[1]) : $!\n");
my $class;
while(<CLASS>){
    chomp;
    if(/^\[CLASS:\s*(\d+)/){$class=$1}
    if(/^PN\s*(.+)/){
	my $kk=$1;
	my @prots=split(/,\s+/,$kk);
	if($wanted_classes){$class=0 unless defined($want{$class});}
	map{$classes{$_}{$class}++}@prots;
    }
    if(/^CA\s+(.+)/){
	my @annot=split(/\s+/,$1);
	$annots{$class}=\@annot;
    }
    elsif(/^CM\t(.+?)\s+(\d+.\d+)/){
	my @annot=(split(/\s+/,$1),$2);
	push @annot, "NO CONSENSUS";
	$annots{$class}=\@annot;
	
    }
}
## If we want each class and its annotations
if ($opts{a}) {
  foreach my $prot (keys(%prots)){
      print "$prot\n";
      foreach my $class (keys(%{$classes{$prot}})){
	  print "\t$class\t@{$annots{$class}}\n";
      }
  }
}
else {    
    print "Classes\n";
    foreach my $prot (keys(%prots)){
	my $line= $prot . " = (";
	foreach my $class (keys(%{$classes{$prot}})){
	    $line.= $class . "::";
	}
	$line=~s/::$//;
	print "$line)\n";
    }
}
