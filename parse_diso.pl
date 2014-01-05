#!/usr/bin/perl -w
use strict;
use Getopt::Std;
use Switch;
use Math::Round;
sub v_say;
sub debug;
require "MY_SUBS.pl";
my (%opts,%prots);
getopts('tw:',\%opts) || do { print "Invalid option, try '$0 -h' for more information\n"; exit(1); };
my $LENGTH=$opts{w}||20;
my $total=$opts{t}||undef;
#usage() unless $ARGV[0];

foreach my $file (@ARGV) {
    open(F,"gawk '(NR>5){printf \$3}END{print \"\"}' $file |")||die "Could not open $file : $!\n";
    my $string=<F>;
    close(F);
    chomp($string);
    my $str_len=length($string);
    my @bb=($string=~/\*/g);
    my $perc=nearest(.001,($#bb+1)*100/$str_len);   
#    my $perc=($#bb+1)*100/$str_len;   
    ## If we just want the overall % of disorder
    if ($total) {
	print $#bb+1 . " " . "$str_len $perc\n";
	next;
    }
 
#    print "SS -$string- ::: $str_len : $p : @bb\n";
    my %mean=();
    for (my $win_len=$LENGTH; $win_len<=100;$win_len++) {
	my $dis=0;
	my $i=0;
#	for ($i=0; $i<=length($string)-$win_len; $i++) {
	i:while ($i<=length($string)-$win_len) {
	    my $p=get_window($i,$win_len,$string);
	    $dis+=$p;
	    # if ($p>=0){
	    # 	$dis+=$p;
	    # }
	    # else {
	    # 	print "aaaaaaaaaaaaaaaaaaa\n";
	    # 	last i;
	    # }
	    $i++;
	}
	if ((length($string)-$win_len+1)>0) {
	    $mean{$win_len}=$dis/$i;
	    print "$win_len\t$mean{$win_len}\n";
#	    print "\$mean{$win_len}=$dis/$i=$mean{$win_len}\n";
	} 
	else {
	    $mean{$win_len}="NaN";
	}
    }
}

#12

sub get_window{
    my $start=shift;
    my $len=shift;
    my $str=shift;
    my $win=substr($str,$start,$len);
    my @bb=($win=~/\*/g);
    my $p=nearest(.001,($#bb+1)*100/$len);
    return($p);
}
