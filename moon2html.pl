#!/usr/bin/perl -w

use strict;
use Getopt::Std;
use Math::Round;
require "MY_SUBS.pl";
my %opts;

getopts('hs:',\%opts);
my $file=$ARGV[0];
$file=~s/\.html/\.txt/;
$file=~/(.+)\.(\w)\.moon/||die("File $file does not look like a .moon file!\n");
my $species=guess_species($1);
my $mode=$2;
my %prec;
if ($mode eq 'P'){
    my %gg;
    open(A,"$ARGV[0]");
    while(<A>){
	my @a=(/(GO:\d+)/g);
	map{$gg{$_}++}@a;
    }
    my @n=keys(%gg);
    close(A);
    open(B,">/tmp/$$.tmp");
    map{print B "$_\n";}@n;
    my $a=`term_precision.pl ~/research/GO/all_three_Jun27-2012.genealogy /tmp/$$.tmp`;
    my @b=split(/\n/,$a);
    map{my @a=split(/\t/); $prec{$a[0]}=nearest(.0001,$a[1])}@b;
    @n=keys(%prec);
    close(B);
}


my %modes=(
    "a" => "Distant from class, look for proteins with annotations dissimilar to those of their class.",
    "b" => "Bridge, look for cases where a protein interacts with 2 partners each in different and distant classes.",
    "B" => "Bridge unique, look for cases where a protein interacts with 2 partners each in different and distant classes that are not otherwise connected. That is, none of the proteins in class1 interacts with any of class2 and vice versa.",
    "i" => "Intersection, look for proteins found at the intersection of two classes annotated to distant GOs. ",
    "I|Inter" => "Intesection unique, look for proteins found at the intersection of two classes annotated to disimilar GOs ONLY.",
    "p" => "Percentage, identify the most common GO (>=50% of all interactors) and look for interactors annotated to a distant GO.",
    "P" => "Pairs, look for proteins whose interaction partners are annotated to dissimilar GOs",
    "d|direct" => "Direct, find ALL interactions between dissimilar GOs"
);

print<<EndOf;
<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN" "http://www.w3.org/TR/html4/loose.dtd">
<html>
<head>
     <title>MoonGO.pl Results</title>
    <link rel=stylesheet type="text/css" href="../general.css" title="Genome">
</head>
<body>
<a href="$file">Download the raw output.</a>
<table class="moon" border="0">
  <tr class="moon_head">
EndOf

my $c=0;
while(<>){
    chomp;
    next if /\#/;
    s/[\)\(]//g;
    ##For mode a 
    s/\s*\t:\t\s*/\t/;
    s/,//g;

    my @a=split(/\t/);

    my @b=split(/\s+/, $a[1]);
    my @c=split(/\s+/, $a[2]);
    ## Check what kind of file this is (which method)
    ## Modes B,b,i,I
    if($mode=~/^[BbiI]$/){
	if ($c==0){
	    print<<EndOf;
    <th>Cand Name</th>
    <th>Class1</th>
    <th>Class1 Card.</th>
    <th>Class1 GO</th>
    <th>GO1 prec.</th>
    <th>Class2</th>
    <th>Class2 Card.</th>
    <th>Class2 GO</th>
    <th>GO2 prec.</th>
    <th>e-value</th>
</tr>
    
EndOf
	}
	print "<tr><td><a href=\"http://www.uniprot.org/uniprot/$a[0]\">$a[0]</a></td>";
	print "<td><a href=\"../$species.BP.clas.html#$b[0]\">$b[0]</a></td><td>$b[1]</td><td><a href=\"http://www.ebi.ac.uk/QuickGO/GTerm?id=$b[2]\">$b[2]</a></td><td>$b[3]</td>";
	print "<td><a href=\"../$species.BP.clas.html#$c[0]\">$c[0]</a></td><td>$c[1]</td><td><a href=\"http://www.ebi.ac.uk/QuickGO/GTerm?id=$c[2]\">$c[2]</a></td><td>$c[3]</td>";
    }
    ## other modes
    elsif($mode eq 'P'){
	if ($c==0){
    print<<EndOf;
    <th>Cand Name</th>
    <th>Interactor1</th>
    <th>Interactor1 GO</th>
    <th>Int1 GO precision</th>
    <th>Interactor2</th>
    <th>Interactor2 GO</th>
    <th>Int2 GO precision</th>
    <th>e-value</th>
</tr>
EndOf
	}
	print "<tr><td><a href=\"http://www.uniprot.org/uniprot/$a[0]\">$a[0]</a></td>";
	print "<td><a href=\"http://www.uniprot.org/uniprot/$b[0]\">$b[0]</a></td><td><a href=\"http://www.ebi.ac.uk/QuickGO/GTerm?id=$b[1]\">$b[1]</a></td><td>$prec{$b[1]}</td>";
	print "<td><a href=\"http://www.uniprot.org/uniprot/$c[0]\">$c[0]</a></td><td><a href=\"http://www.ebi.ac.uk/QuickGO/GTerm?id=$b[1]\">$c[1]</a></td><td>$prec{$c[1]}</td>";
    }
    elsif($mode eq 'a'){
	if ($c==0){
    print<<EndOf;
    <th>Cand Name</th>
    <th>Cand GO</th>
    <th>Cand GO prec</th>
    <th>Class</th>
    <th>Class GO</th>
    <th>Class GO prec</th>
    <th>e-value</th>
</tr>
EndOf
	}
	print "<tr><td><a href=\"http://www.uniprot.org/uniprot/$a[0]\">$a[0]</a></td>";
	print "<td><a href=\"http://www.ebi.ac.uk/QuickGO/GTerm?id=$b[0]\">$b[0]</a></td><td>$b[1]</td>";
	print "<td><a href=\"../$species.BP.clas.html#$c[0]\">$c[0]</a></td>";
	print "<td><a href=\"http://www.ebi.ac.uk/QuickGO/GTerm?id=$c[1]\">$c[1]</a></td><td>$c[2]</td>";	

    }
    print "<td>$a[$#a]</td></tr>\n";

    $c=1;
}
print "</table></body></html>\n";

