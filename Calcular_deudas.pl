#!/usr/bin/perl -w

use Math::Round;
use strict;
use Getopt::Std;
my %opts;
getopts('a:d:jkbmh',\%opts);
&usage() if $opts{h};

my $juli=$opts{j}||undef;
my $karo=$opts{k}||undef;
my $deuda=$opts{d}||0;
my ($t1,$t2)=(1,1);
my ($deuda_t1, $deuda_t2)=0;
my @timeData = localtime(time);
my $month=$timeData[4]+1;
die("Change the date stuff, new year!\n") if $month==0;
if($juli){
 	$t2=1458.98;
 	$t1=561.35;
 }
elsif($karo){
    $t1=201.11;
    $t2=767.40;
}
my $both=$opts{b}||undef;
if($both){
	$t2=1458.98+237.65;
 	$t1=561.35+719.57;
}
my $debt=$t1+$t2;
my $total=$t1+$t2;
my $mm=0;
my $k=0;
my $last_total=0;
my $n=0;
print "Date\t\t  T1\tT2\tTotal\tPaid    Min\tLeft\n";

while($debt>0){
    my $year=$timeData[5]+1900+$k;
    $mm=($timeData[4]+1+$n) % 12;
    if ($mm==0){
	$mm=12 ;
	$k++;
    };
    $n++;
    $t1=nearest(.01,$t1);
    $t2=nearest(.1,$t2);
    
    my ($deuda_t1,$deuda_t2)=&calculate($t1,$t2);
    my $min=$deuda_t1+$deuda_t2;    
    $last_total=$total;
    $total=$t1+$t2;
    print "$mm/$year\t\t$t1\t$t2\t";
    my $paid;
    # if($deuda_t2+$deuda_t1<=100){
    # 	  $deuda_t2=$deuda_t1=50;
    #   };
    # if( $t1<=0){
    # 	$deuda_t2+=50;
    # }
    if($opts{a}){
	$paid=$opts{a};
    }
    elsif($opts{m}){
	$paid=$deuda_t1+$deuda_t2;
    }
    else{
	if(($deuda_t1+$deuda_t2)<100){
	    $paid=100;
	}
	else{
	    $paid=$deuda_t1+$deuda_t2 ;
	}
    }
    if($debt<100){$paid=$debt; print "\t(dep : $debt)\t"}
#    my $left=nearest(.01,$total-$this_month);
    $debt=$debt-$paid;
    print "$total\t$paid($deuda_t1+$deuda_t2)\t$min\t$debt\n";
  
    if($t1<=0){
    $t2-=$paid;
    }
    else{
	$t1-=($paid/2);
	$t2-=($paid/2);
    }
    $debt=$t1+$t2;
    $t1=0 if $t1<0;
    $t2=0 if $t2<0;

}
###########################################################3

sub calculate{
    $t1=shift;
    $t2=shift;
    if($juli){
	my $pp=$t2 * 0.03;
	if($pp>30){
	    $deuda_t2=$pp;
	}
	else{$deuda_t2=30;}
	#print "dd : $deuda_t1\n";
	$deuda_t1=$t1 * 0.10;
    }
    elsif($karo){
	$deuda_t1=$t1 * 0.10;
	$deuda_t2=$t2 * 0.10;
	
    }
    elsif($both){
	$deuda_t1=$t1 * 0.10;
	$deuda_t2=$t2 * 0.10;	
    }
    return($deuda_t1,$deuda_t2);
}



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
