#!/usr/bin/perl5.8.0 -w
## Script to copy files from the iPod using gnupod_search.pl.
## Run gnupod_search.pl witn --view=u and pipe through ipodcopy

##                           EXAMPLE                                          ## 
## gnupod_search.pl -m /mnt/iPod -t"title" --view=u | ipodcopy /music/dirname ##
##                                                                            ##
use strict;

my $dest = $ARGV[0] || die("need at least one argument");

while(<STDIN>)
{
    next unless /^\//;
    /^(.*\d\/)(.*mp3)/; 
    my $n = $2; 
    my $p = $1; 
    $n =~ s/F(\d\d)/f$1/; 
    my $nn = $n; 
    $nn =~ s/\s/_/g; 
    system("cp -v \"$p$n\" $dest");
}
