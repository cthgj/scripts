#!/usr/bin/perl

sub fac {
  $_[0]>1?$_[0]*fac($_[0]-1):1;
}

#print fac($ARGV[0]) . "\n";


my $a = (fac(25) * fac(75) *  fac(70) * fac(30) )/( fac(100) * fac(15) *  fac(10) * fac(55) * fac(20) );
# $a = (fac(16) * fac(84) *  fac(79) * fac(21) )/( fac(100) * fac(15) *  fac(1) * fac(64) * fac(21) );

print "a : $a\n";
