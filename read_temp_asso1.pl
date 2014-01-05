#!/usr/bin/perl -w
use strict;
use Getopt::Std;
use Math::Round;

my (%all_terms,%ancestors,%stats,%gos,%go_count,%opts,%counted);
getopts('pv',\%opts);
my $verbose=$opts{v}||undef;
my $print_to_file=$opts{p}||undef; ## do not run rscript, just print it
my $subonto="P";
my $gaf_file=$ARGV[0]||"gene_association.goa_human";
my $geneology_file=$ARGV[1]||"biological_process.genealogy";
my $prot_name="Bob";



open(A,"$geneology_file")||die("cannot open $geneology_file: $!\n");
while(<A>){
    chomp;
    my @terms=split(/\t/);
#    map {$all_terms{$_}=1}@terms;
    my $child=shift(@terms);
    @{$ancestors{$child}}=@terms;
}
close(A);
my $silivar="!";
open(A,"$gaf_file")||die("cannot open $gaf_file: $!\n");
while(<A>){
    next if /^$silivar/; ## silly thing cause emacs screws up the syntax highlihting when /!/;
    chomp;
    my @tmparray=split(/\t/);
    next unless $tmparray[8] eq $subonto;
    next unless $tmparray[3]=~/^$/;  ## skip the NOT annotations
    my $nn=$tmparray[2]; ## $nn == protein name
    $gos{$nn}{$tmparray[4]}++; ## $tmparray[4] == GO
}
close(A);
my @prots=keys(%gos);
my $total=scalar(@prots); ## total is the number of annotated proteins
my $c=0;

## Inherit parent terms
foreach my $prot (@prots) {
    $c++;
    printf(STDERR "[Proteins: $c of $total]\r") if $c % 1000 == 0;
    foreach my $goterm (keys(%{$gos{$prot}})){
	$go_count{$goterm}++;
	map{$gos{$prot}{$_}++; $go_count{$_}++}@{$ancestors{$goterm}};
    }
## @GOs is now all the terms of $prot
    my @GOs=keys(%{$gos{$prot}});
    for (my $n=0; $n<=$#GOs;$n++){
	my $jojo=$n+1;
	for (my $k=$n+1; $k<scalar(@GOs);$k++){
	    #print "aa : $n:$k " . scalar(@GOs) ." dd\n";
	    my $go1=$GOs[$n];
	    my $go2=$GOs[$k];
	    $go1=~/GO:(\d+)/;
	    my $l=$1;
	    $go2=~/GO:(\d+)/;
	    my $k=$1;
	    my $GOpair;
	    $k>$l ? ($GOpair=$go1 . "xxx" . $go2) : ($GOpair=$go2 . "xxx" . $go1);
	    $stats{$GOpair}++;
	}##end for my $k
    }##end for my $n
}##end foreach prot

## OK, now deal with all go pairs
## R stuff
## $stats{$GOpair} : number of times these 2 GOs coincide
## $go_count{$go1} : number of prots with go1 
## $total-$go_count{$go1} : number of prots without go1 
## $go_count{$go2} : number of prots with go2

## need number of go2!go1
## $go_count{$go2}-$stats{$GOpair}


## Try forking
my @bob=keys(%stats);
my $num=round($#bob/8);
my $kl=1;
my @children;

for (1..4){
    my $pid = fork();
    if ($pid) {
	push(@children, $pid);
	print STDOUT "PPID : $pid\n";
	$num=$num*$kl;
	$kl++;
    }
    elsif ($pid == 0) {
	for (my $n=0;$n<=$num;$n++){
	    printf(STDERR "[GOs: $n of $k]\r") if $k % 100 == 0;
	    my $GOpair=$bob[$n];
	    my ($go1,$go2)=split(/xxx/,$GOpair);
		print STDOUT "hhyper=phyper($stats{$GOpair},$go_count{$go1}," . ($total-$go_count{$go1}) .",$go_count{$go2},lower.tail = FALSE)\n";
	    print STDOUT "lhyper=phyper($stats{$GOpair},$go_count{$go1}," . ($total-$go_count{$go1}) .",$go_count{$go2},lower.tail = TRUE)\n";
	    print STDOUT "hypers=c(\"$go1\",\"$go2\",hhyper,lhyper)\n";
	    print STDOUT "write(hypers,file=\"/ptitbbackup/cchapple_temp/hyper\",ncolumns=4,append=TRUE,sep=\"\t\")\n";
	}
	exit(0);
    }
}
foreach (@children) {
    waitpid($_, 0);
}

print "num : $num : $#bob\n"; exit();






printf(STDERR "[Proteins: $c of $total]\r");
print STDERR "\n-----------------\n";
my $n=0;
my $k=scalar(keys(%stats));
foreach my $GOpair (keys(%stats)){
    $n++;
    printf(STDERR "[GOs: $n of $k]\r") if $k % 100 == 0;
    my ($go1,$go2)=split(/xxx/,$GOpair);                                     
#go2               $total-$go_count{$go1}
    print STDOUT "hhyper=phyper($stats{$GOpair},$go_count{$go1}," . ($total-$go_count{$go1}) .",$go_count{$go2},lower.tail = FALSE)\n";
    print STDOUT "lhyper=phyper($stats{$GOpair},$go_count{$go1}," . ($total-$go_count{$go1}) .",$go_count{$go2},lower.tail = TRUE)\n";
    print STDOUT "hypers=c(\"$go1\",\"$go2\",hhyper,lhyper)\n";
    print STDOUT "write(hypers,file=\"/ptitbbackup/cchapple_temp/hyper\",ncolumns=4,append=TRUE,sep=\"\t\")\n";
    
}
printf(STDERR "[GOs: $n of $k]\r");
print STDERR "\n";
