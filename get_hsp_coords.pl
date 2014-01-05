#!/usr/bin/perl -w
##
## This script will return the the region (start-end) of each
## subject sequence which gave an HSP. So, if there are 3 HSPs for subject x
## eg 1-30 56-79 and 87-100, the script will return queryname : sub (1, 100, strand) 
use strict;
my $subject = "a";
my $frame;
my $sstart;#9999999999999999;
my ($skip, $query, $qstart, $qend, $send) = (0,0,0,0,0); 
my $hh;
my $sub_num = 0;
my %fr;
while(<>)
{
    if (/y=\s*(.+?)\s/)
    {
	$sub_num = 0;
	$query=$1;
	$subject = "a";
    }
    elsif(/>(.+?)\s/)
    {
	$sub_num++;
	$hh=0;
	if($subject ne $1 && $subject ne "a")
	{
	    print "$query $subject $sstart $send $frame\n" if $sub_num < 13;
	}
	$subject = $1;
	($qstart, $qend, $sstart, $send) = (0,0,0,0,0); 	
    }
    elsif(/Expect.*?=\s*(.*)\s/)
    {
	$skip=0;
	my $evalue = $1;
	if ($evalue =~ /^\s?e/)
	{
	    $evalue = "1" . $evalue;
	}
	$skip=1 if $evalue > 0.001;
    }
    elsif(/Identities\s=\s.*?\((\d+)%\)/)
    {
	$skip=1 if $1 < 50;
    }
    elsif(/Frame\s=\s([+-])/ && $hh == 0)
    {
	
	$fr{$subject} = $1 if $hh == 0;
	$frame=$1;
	$skip = 1 if $frame ne $fr{$subject};
	$hh=1;
	$sstart = 9999999999999999;

    }
    elsif(/y:/)
    {   
	next if $skip==1;
	/:\s?(\d+)[^\d]*(\d+)/ || die();
	if($frame eq "+")
	{
	    $qstart = $1 unless $1 > $qstart;
	    $qend = $2 unless $2<$qend;
	}
	else
	{
	    
	}
	
    }
    elsif(/jct:/)
    {
	next if $skip==1;
	/:\s?(\d+)[^\d]*(\d+)/ || die();
	if($frame eq "+")
	{
	    $sstart = $1 unless $1 > $sstart;
	    $send = $2 unless $2<$send;	
	}
	else
	{
	    $sstart = $2 unless $2 > $sstart;
	    $send = $1 unless $1<$send;	
	}
	
    }

}
print "$query $subject $sstart $send $frame\n" if $sub_num < 13;





