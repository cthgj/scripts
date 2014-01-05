#!/usr/bin/perl -w

###
# Simple little script which takes pep's parseblast.out file and turns it intosomething usefull.
# 
###

use Getopt::Std;
use strict();

my %opts;
getopts('p',\%opts) || die();

my (%score, %exp, %name, %sel) =();
my ($score, $exp, $name, $sel) = 0;
my $noplants = $opts{p} || 1;

system("gawk \'{print \$7\" \" \$9\" \"\$12\" \"\$18}\' $ARGV[0] > tmp");


open(TMP,"tmp") || die "cannot open tmp : $!\n";
while(<TMP>)
{
    /(.*?)\s(.*?)\s(.*?)\s(.*?)\n/og;
    print STDERR "$& is \$&\n";
#    if ($3 eq $name){next;}
    if($3 =~ /arab || barley || cotton || grape || ice\.plant || lettuce || maize || pine || rice || potato || rice || soybean || sunflower || wheat/){print STDOUT "###$&\n"; next;}
    else{
	print "%%%%%%%%%%%%%%%%%%%%%%%%%%%\n";
#	$score = $1;
#	$exp = $2;
	$name = $3;
#	$sel = $4;    
	
	$score{$3} = $1;
	$exp{$3} = $2;
	$name{$3} = $3;
	$sel{$3} = $4;     
    }
}

close(TMP);
unlink("tmp");

my @keys = keys(%name);


map {print STDOUT "$sel{$_} $name{$_}   score: $score{$_} e-value : $exp{$_}\n"} @keys

