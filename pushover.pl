#!/usr/bin/perl -w
use strict;
use LWP::UserAgent;

my $title=`echo \$HOSTNAME`|| "none";
LWP::UserAgent->new()->post
    (
     "https://api.pushover.net/1/messages", 
     [
      "token" => "f0ZI0bXdKoUNjN4VS6caePiUcxbKKW",
      "user" => "7M7bVuK1ievp6J5ztYOxJxGrGFLkIu",
      "message" => "$ARGV[0]",
      "title" => "$title",
     ]);
