#!/usr/bin/perl -w
use strict;
use Number::Bytes::Human qw(format_bytes);
use v5.10;
use Math::Round;
my $max_mem=0;
my %MAX;
my $name=$ARGV[0];
## If we have given a process name to monitor
while (1){
    my @pids;
    if ($name=~/^\d+$/) {
	$pids[0]=$name
    }
    else {
	@pids=split(/\n/,`pgrep $name`);
    }
    my $running=$#pids;
    exit(0) unless $pids[0];
    foreach my $pid (@pids) {
	my $is_running=`ps aux | grep -v grep | grep -wc $pid `;
	chomp($is_running);
	unless ($is_running){
	    $running--;
	    next;
	}
	my $com=`ps -o args -p $pid --no-headers`;
	chomp($com);
	my $size=`ps  -o size  -p $pid --no-headers`;
	chomp $size;
	my $bytes=(1024 * 1024);
	my $gigs=$size/$bytes;
	$gigs=nearest(0.01,$gigs);
	$MAX{$pid}//=0;
	my $time=`date +'%d/%m %H:%M'`;
	chomp($time);
	if ($size > $MAX{$pid}){
	    $MAX{$pid}=$size;
	    say "Max so far for command \"$com\" at $time : $size, $gigs GB";
	}
    }
}



