#!/usr/bin/perl -w

use strict;

my $speed=`cat /proc/i8k | cut -d " " -f 6`;
my $r_speed=`cat /proc/i8k | cut -d " " -f 8`;
my $text;
chomp($speed);
chomp($r_speed);
if($speed==2){
    $text="Fast";
}
elsif($speed==1){
    $text="Slow";
}
elsif($speed==0){
    $text="Off";
}
else{
    $text="$speed, Huh?";
}

print "$text [ $r_speed ]";
