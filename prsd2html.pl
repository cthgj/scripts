#!/usr/bin/perl -w

use strict;
use Getopt::Std;
require "MY_SUBS.pl";
my %opts;
getopts('hg:s:',\%opts);
my $species=$opts{s}||'human';
my %gos2terms;
my $gold_file=$opts{g}||$ENV{"HOME"} . "/research/moonlight/data/gold_simple.txt";
print<<EndOf;
<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN" "http://www.w3.org/TR/html4/loose.dtd">
<html>
<head>
    <link rel=stylesheet type="text/css" href="../general.css" title="Genome">
    <script src="../sorttable-small.js"></script>
</head>
<body>

<div class="scrollWrapper">
EndOf
my $b=0;
my %seen;
my %gold;
############################
# Read gold standard prots #
############################
my $fh=check_file($gold_file,"r");
while (<$fh>) {
    next if $.==1;
    /^(.+?)\s/;
    $gold{$1}++;
}


while(<>){
    chomp;
    next if /\#.*TESTS/;
    s/,(\w)/, $1/g;
    ##Parse descriptor line
    if (s/^#+\s*//) {
	print "<table class=\"moon sortable\" style =\" border-collapse:collapse;  margin-bottom:5px\"><tr >";
	my @line=split(/\t/);
	map{
	    print "<th>$_</th>";
	}@line;
	if ($line[0]=~/pair/) {
	    print "<th>Prob.</th>";
	}

	print "</tr>\n";
    }
    ## If this is not the prot name
    elsif (s/^\s+//){
	print "<tr>";
	s/\"//g;
	s/(GO:\d+)/<a href=http:\/\/www.ebi.ac.uk\/QuickGO\/GTerm?id=$1>$1<\/a>/g;
	s/(Class.: )(\d+)/<a href=\.\.\/$species.BP.clas.html\#$2>$1$2<\/a>/g;
	my @line=split(/\t/);
	map{
	    print "<td>$_</td>";
	}@line;
	print "</tr>\n"
    }
    ## If this is the prot name or if these are GOs
    else{
	my @line=split(/\t/);
	####################
        # If these are GOs #
        ####################
	if (/^GO:/) {
	    ## If these are pairs
	    if (/(GO:\d+_.+?)\s/) {
		my ($o,$f);
		$o=$f=$1;
		$o=~s/://g;
		## Link term descriptions to quickGO
		my @gos=split(/_/,$line[0]);
		$line[8]="<a href=http://www.ebi.ac.uk/QuickGO/GTerm?id=$gos[0]>$line[8]</a>";
		$line[9]="<a href=http://www.ebi.ac.uk/QuickGO/GTerm?id=$gos[1]>$line[9]</a>";
		## Link GO pairs to their candidate lists
		$line[0]=~s/$f/<a href=\"$o.html\">$f<\/a>/;

	    }
	    ## Link single GOs to their quickGO page
	    elsif (/(GO:\d+)/) {
		my ($o,$f);
		$o=$f=$1;
		$o=~s/://g;
		$line[7]="<a href=http://www.ebi.ac.uk/QuickGO/GTerm?id=$line[0]>$line[7]</a>";
		$line[0]=~s/$f/<a href=\"$o.html\">$f<\/a>/;

	    }
	    print "<tr>";
	    map{
		print "<td>$_</td>";
	    }@line;
	    print "</tr>\n";
	}
	else{
	    ## Check for gold
	    if (defined($gold{$_})) {
		print "</tr><tr class=moon_gold><th  colspan=\"9\" ><a href=\"http://www.uniprot.org/uniprot/$_\">$_ (Gold)</a></th></tr>\n";	
	    }
	    else {
		print "</tr><tr class=moon_head><th  colspan=\"9\" ><a href=\"http://www.uniprot.org/uniprot/$_\">$_</a></th></tr>\n";	
	    } 
	}
	
    }
}
print "</tr></table></body></html>\n";

