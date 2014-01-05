#!/usr/bin/perl -w

my (@cm,@ca); 
open(A,"$ARGV[0]"); 
while(<A>){
    if(/^CA\s+(.+?)$/){
	@a=split(/\s+/,$1); 
	push @ca,@a;
    } 
    if(/^CM\s+(.+?)\s+\d/){
	@a=split(/\s+/,$1); 
	push @cm,@a;
    } 
}
close(A);
open(A,"$ARGV[0]"); 

    while(<A>){
	if(/^CA\s+(.+?)$/){
	    chomp;
	    print "CA\t";
	    $k=$1; 
	    @a=split(/\s+/,$k); 
	    $num=$#a+1; 
	    my %seen;
	    for(my $i=0; $i<=$num; $i++){
		my $annot=$cm[int(rand($#cm))];
		while(defined($seen{$annot})){
		    $annot=$cm[int(rand($#cm))];
		}
		print "$annot ";
	    }
	    print "\n";

	}
	elsif(/^CM\s+(.+?)(\s+\d.+)$/){
	    my $rest=$2;
	    @a=split(/\s+/,$k); 
	    print "CM\t";
	    $num=$#a+1; 
	    my %seen;
	    for(my $i=0; $i<=$num; $i++){
		my $annot=$cm[int(rand($#cm))];
		while(defined($seen{$annot})){
		    $annot=$cm[int(rand($#cm))];
		}
		print "$annot ";
	    }
	    print "$rest\n";

	}
	else{print}
    }

