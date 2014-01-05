#!/usr/bin/perl -w

use strict;

my $map=$ARGV[0];
my $graph=$ARGV[1];

my (%newname,%interactions,%all);

unless($ARGV[0]){
    print STDERR "This script will take a network.gr file and a map file (old_name\\tnew_name\\taccession) and replace any names in the network that have changed with their new versions. \n\nUSAGE: $0 map graph.gr\n\n";
    exit(0);
}

open(MAP,"$map")|| die("Could not open $map :$!\n");
while(<MAP>){
    chomp;
    my @a=split(/\s+/,$_);
    ## $accs{old_acc}=new_acc
    $newname{$a[0]}=$a[1];
}
close(MAP);
open(GRAPH,"$graph")|| die("Could not open $graph :$!\n");
while(<GRAPH>){
    next if /^\d+$/;
    my $line=$_;
    my @a=split(/\s+/,$line);
    ## ids{acc}=id, eg $id{A0A111}=PAFA_MYCTU
    # map{defined($newname{$_}) && do {$line=~s/$_/$newname{$_}/g;}@a;

    # @a=split(/\s+/,$line);
#    unless ($line=~/deleted/){
	$all{$a[1]}=$all{$a[0]}=1;
	if (exists($interactions{$a[1]})){
	    push @{$interactions{$a[1]}},$a[0];
	}
	else{
	    push @{$interactions{$a[0]}},$a[1];
	}
    
#    print "$line";
}
close(GRAPH);
my @out;
my (%merged,%deleted,%int);
my $d=0;
bait:foreach my $bait (keys(%interactions)){
    $newname{$bait}=$bait unless defined($newname{$bait});
    if ($newname{$bait}=~/deleted/){
	$deleted{$bait}+=scalar(@{$interactions{$bait}});
	$d+=scalar(@{$interactions{$bait}});
	next bait;
    }
    target:foreach my $target (@{$interactions{$bait}}){
	$newname{$target}=$target unless defined($newname{$target});
	die("huh?\n") if $newname{$bait}=~/deleted/;
	if ($newname{$target}=~/deleted/){
	    $deleted{$target}++;
	    $d++;
	    next target ;
	}
        if($newname{$bait} ne $bait){
	    $merged{$bait}=$newname{$bait};
	}
	if($newname{$target} ne $target){
	    $merged{$target}=$newname{$target};
	}
	$int{"$newname{$bait}\t$newname{$target}"}++;

    }
}

print STDOUT scalar(keys(%int)) . "\n";
map{print "$_\n"}keys(%int);
### print stats   32385 => 31923
my @ll;
print STDERR "\n\n";
map{
    push @ll, $_, if defined($all{$merged{$_}})
}keys(%merged);

print STDERR scalar(@ll) . " duplicate entries were merged (@ll)\n\tand\n$d interactions involving a deleted prot were deleted\n";

#map{print STDERR "$_\n" if $int{$_}>1}keys(%int);
