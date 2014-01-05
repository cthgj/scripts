#!/usr/bin/perl -w
use strict;
use Getopt::Std;

my @line;
my($sp,$k,$length,$prog,$type,$st,$en,$score,$str,$ph,$name,$diff);
my $prot;# = $ARGV[0];
#$prot =~ s/(.*)_(.*)\.gff/$2/;
my $ff=0;
my %p;
my @lines;
open(A,"$ARGV[0]");
while(<A>)
{
    my($a,$prog,$type,$st,$en,$score,$str,$ph,$prot) = split;
    $a =~ /^(.*?)_.+?\((\d+),/; 
    ($sp,$length) = ($1, $2);
    if($str eq "+")
    {
	$st = $length + $st;
	$en = $length + $en;
    }
    else
    {
	$st = $length + $en;
	$en = $length + $st;
    }
    
    if(defined($p{$st}{$sp})){ 
	if($st < $p{st}{$sp}){
	    $p{$sp} = $prot;
	    $p{st}{$sp} = $st;    
	    
	} 
    }
    else{$p{$sp} = $prot; $p{st}{$sp} = $st; }
    $score = "0.5";
    my $pp =  "$sp\t$sp\t$type\t$st\t$en\t$score\t$str\t$ph\t$prot\n";
    push @lines, $pp;
}
close(A);
#open(A,"$ARGV[0]") || die("shit\n");
#open(A,"sort -k9 $ARGV[0] |") || die("shit\n");
#while(<A>)
map
{
    my($a,$prog,$type,$st,$en,$score,$str,$ph,$prot) = split;
    # $a =~ /^(.*?)_.+?\((\d+),/; 
#     print "aaaaaaaaaaaa : $a : $_\n"; die();
#     ($sp,$length) = ($1, $2);
#     if($str eq "+")
#     {
# 	$st = $length + $st;
# 	$en = $length + $en;
#     }
#     else
#     {
# 	$st = $length + $en;
# 	$en = $length + $st;
#     }
#     $score = "0.5";

   #  $k = $p{st}{$sp} -1;
#     $st = 1;
#     $en = $en - $k  ;
#     $ff++;
	if(($prot eq "$p{$sp}") && ($ff == 0))
	{
	    $k = $st -1;
	    $st = 1;
	    $en = $en - $k  ;
	    $ff++;
	}
	else#(($prot eq "CG31715") && ($ff != 0))
	{
	    $st = $st - $k;
	    $en = $en - $k;
	}

    

    # if(($prot eq "$p{$sp}") && ($ff == 0))
#     {	
# 	$k = $st -1;
# 	$st = 1;
# 	$en = $en - $k  ;
# 	$ff++;
#     }
#     else#(($prot eq "CG31715") && ($ff != 0))
#     {
# 	$st = $st - $k;
# 	$en = $en - $k;
#     }
    $ff= 0 unless $prot eq "$p{$sp}";
    
# 	if(($prot eq "CG31715") && ($ff == 0))
# 	{
# 	    $k = $st -1;
# 	    $st = 1;
# 	    $en = $en - $k  ;
# 	    $ff++;
# 	}
# 	else#(($prot eq "CG31715") && ($ff != 0))
# 	{
# 	    $st = $st - $k;
# 	    $en = $en - $k;
# 	}
# 	$ff= 0 unless $prot eq "CG31715";
    
print STDOUT "$sp\t$sp\t$type\t$st\t$en\t$score\t$str\t$ph\t$prot\n";
}@lines;
