#!/usr/bin/perl
# this does not work; perl hasn't seen the 'hello'
# subroutine yet

hello();

sub hello
{
  print "Hello, world.\n";
}

