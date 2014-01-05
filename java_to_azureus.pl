#!/usr/bin/perl -w
use strict;
## Silly script to change PID/java in netstat output to "vuze" for my conky

my @fields;
my %h;
my $vuze_pid;
while(<>){
    next unless /ESTABLISHED/;
    @fields=split();

    $fields[4]=~/(.+?)\.local/ && do{
	my $k=$1;
	$fields[8]=$k;
    };
    if(/(\d+).java\s*$/){
	if($vuze_pid){
	    if($1 eq $vuze_pid){
		$fields[8]=~s/java/vuze/;
	    }
	}
	else{
	    $vuze_pid=$1;
	    open(A,"ps aux | grep $vuze_pid | egrep 'azureus|vuze' |") && do{
		$fields[8]=~s/java/vuze/;
		close(A);
	    }
	}
    }
#    for (my $n=0; $n<scalar(@fields); $n++){print "$n : $fields[0]\n";}
    $fields[8]=~s/\d+\///;
    if($fields[8] eq "-"){
	if($fields[4]=~/192.168.0.72/){
	    $fields[8]="badabing";
	}
	else{$fields[8]="unknown";}
	
    }
#    print "$fields[8]\n";
    $h{$fields[8]}++;
}
map{print "\$alignr $h{$_} $_\n"}keys(%h);
#$alignr${execi 10 /home/cchapple/scripts/java_to_azureus.pl ~/netstat  | cut -d: -f1 | cut -d\/ -f2 | sort | uniq -c | sort -nr}
