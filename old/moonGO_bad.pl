#!/usr/bin/perl

my $cmnd=`ps aux | grep $$ | grep -v grep`;
chomp($cmnd);
#$cmnd=~s/^.+(moonG.+?)$/$1/;
print "\n$cmnd : " . scalar(@ff) . " gold standard proteins were f\n";
