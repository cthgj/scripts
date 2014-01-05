#!/usr/bin/perl -w

use Math::Round;
use strict;
use Getopt::Std;
my %opts;
getopts('a:d:jkmh',\%opts);
&usage() if $opts{h};
sub usage{
    print STDERR<<EndOfHelp;
    -h : print this help and exit
    -k : Calculate for Karolos
    -j : Calculate for Juli
    -m : always pay minimum 
    -a : always paid specified amount (def:100)
EndOfHelp
exit();
}

my $juli=$opts{j}||undef;
my $karo=$opts{k}||undef;
my $deuda=$opts{d}||0;
my ($t1,$t2)=(1,1);
my ($deuda_t1, $deuda_t2)=0;

print "Date\t\t  T1\tT2\tTotal\t        Paid    Min\tLeft\n";
my @timeData = localtime(time);
my $month=$timeData[4]+1;
die("Change the date stuff, new year!\n") if $month==0;
 if($juli){
 	$t2=1458.98;
 	$t1=561.35;
 }
elsif($karo){
    $t1=237.65;
    $t2=719.57;
}
else{die("Pa juli o pa karo?\n");}
my $total=$t1+$t2;
my $mm=0;
my $k=0;
my $last_total=0;
my $n=0;
while(($t1+$t2)>0){
    my $year=$timeData[5]+1900+$k;
    $mm=($timeData[4]+1+$n) % 12;
    if ($mm==0){
	$mm=12 ;
	$k++;
    };
    $n++;
    $t1=nearest(.01,$t1);
    $t2=nearest(.1,$t2);
    my ($min_this_month,$t1,$t2)=&calculate($t1,$t2);
    my $padding=6-length($t2);
    my $pad="";
    for(my $n=$padding;$n>0;$n--){
	$pad .=" ";
    }
    
    $padding=10-length($t2+$t2);
    $pad="";
    for(my $n=$padding;$n>0;$n--){
	$pad .=" ";
    }
    my $this_month;
    print "$mm/$year\t\t$t1\t$t2";
    print "\t" . ($t1+$t2) . "$pad\t";

    $last_total=$total;
    
    if($opts{a}){
	$this_month=$opts{a};
    }
    elsif($opts{m}){
	$this_month=$min_this_month;
    }
    else{
	$min_this_month<100 ? ($this_month=100) : ($this_month=$min_this_month) ;
    }
    $this_month=100 if $this_month<100;
    my $o=($t1+$t2)-$this_month;
    #print "  :($t1+$t2)-$this_month=$o\t";
    my $left=nearest(.01,($t1+$t2)-$this_month);
    print "$this_month\t$min_this_month\t$left\n";

}

sub calculate{
    $t1=shift;
    $t2=shift;
    if($juli){
	my $pp=$t1 * 0.03;
	if($pp>30){
	    $deuda_t1=$pp;
	}
	else{$deuda_t1=30;}
	if(($t1 * 0.1)<50){ $deuda_t1=50 unless $opts{m}}
	if(($t2 * 0.1)<50){ $deuda_t2=50 unless $opts{m}}
	#print "dd : $deuda_t1\n";
	$deuda_t2=$t2 * 0.10;
    }
    elsif($karo){
	$deuda_t1=$t1 * 0.10;
	$deuda_t2=$t2 * 0.10;
	if(($t1 * 0.1)<50){ 
	    $deuda_t1=50 unless $opts{m} ;
	    
	    
	}
	if(($t2 * 0.1)<50){ $deuda_t2=50 unless $opts{m}}

    }
    if($deuda_t1<=0){
	$deuda_t2+=50;
    }
    else{
	$deuda_t1=$deuda_t2=50 if ($deuda_t1+$deuda_t2)<100;
    }
    $deuda=nearest(.1,$deuda_t1 + $deuda_t2);
    $t1-=$deuda_t1;
    $t2-=$deuda_t2;
    unless($opts{m}){
	$deuda_t1=$t1 if $t1<=100;
	$deuda_t2=$t2 if $t1<=100;
    }
    $t1=0 if $t1<0;
    $t2=0 if $t2<0;
   # print STDERR "B $t1:$t2 ::: $deuda_t1:$deuda_t2\n";
    return($deuda,$t1,$t2);
}

# $i = &round($f) - round to nearest integer
# sub round {
#     int($_[0] + ($_[0] >=0 ? 0.5 : -0.5)); 
# }


## juli  1458.98  561.35;
## karolos 217,65, 604,59
