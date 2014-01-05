#!/usr/bin/perl -w
my %interactions;


foreach my $graph(@ARGV){
    open(A,$graph);
    while(<A>){
	next if /^\d+$/;
	my @a=split(/\s+/,$_);    
	die("$. : $_") unless $a[1];
	push @{$interactions{$graph}{$a[0]}}, $a[1];
	push @{$interactions{$graph}{$a[1]}}, $a[0];

    }
    close(A);
}
print STDERR "Read ". scalar(@ARGV) . " graphs\n";

foreach my $graph(@ARGV){
    my $count=0;
    foreach my $prot (keys(%{$interactions{$graph}})){
	$count+=scalar(@{$interactions{$graph}{$prot}});
	
    }
   print STDERR "$graph : $count interactions, " . scalar(keys(%{$interactions{$graph}})) . " proteins\n"; 
}
