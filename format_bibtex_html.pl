#!/usr/bin/perl -w
use strict;

my $c=0; # counter
my $a=-1; # counter
my @file=<STDIN>;
foreach my $line (@file){
    $a++;
    $line =~ s/<body>/<body><ul class="biblist">/;
    $line =~ s/<p><a.+?<\/a>/<li>/;

    if($line =~ /<em>/){
	$line =~ s/<em>/<br><em>/g;
	$c++;
    }
    if($line =~ /<li>/){$c=0};
    $line =~ s/\s([\d\)\(]+):/ <b>$1<\/b>:/;
    ## Add a new line after the list of names
    ## but only the list of names ($c)
    $line =~ s/(,\s*and\s*.+?\.)/$1<br>/ && do {$c++} if $c==0; 
    
    ##Sometimes there is an 'and' at the end of the line
    ## and that screws up the later substitution (list of names)
    if ($line =~ /and\s*\n/ && $file[$a+1]=~ /(^[^\.]+\.)/){
	$file[$a+1]=~s/(^[^\.]+\.)/$1<br>/ if $c==0;
    }
    ##again, sometimes just one name on new line.
    if($line=~/^\s*[^\s]+\.\s*$/){
	$line=~s/^\s*([^\s]+?\.)\s*\n/$1<br>/ if $c==0;
    }
    $line=~s/<br>\s*<br>/<br>/ig;
    print $line;
}
