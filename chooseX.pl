#!/usr/bin/perl -w

use strict;
use feature "switch";


my $layout=`disper -l | grep -oP 'Philips|Dell'`;
chomp($layout);
print "$layout\n";

my $command="/usr/bin/Xorg -br -verbose -audit 0 -novtswitch ";
given($layout){
    when (/Dell/){
	$command.="-config xorg.conf.dell";
    }
    when (/Philips/){
	$command.="-config xorg.conf.philips";
    }
}

system("$command &");
