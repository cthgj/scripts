#!/usr/bin/perl -w

use strict;

my (%scores,%hits,%best);
while(<>){
    my @a=split(/\s+/); 
    $a[11]=~/.+_(.+?)$/; ##query
    my $q_sp=$1;  ## query species
    $a[17]=~/.+_(.+?)$/; ##subject
    my $s_sp=$1; ## subject species
    if(($q_sp ne $s_sp) && ($a[5]>=50)){ 
	## If we have already found one good hit for this prot/species
	if(exists($hits{$a[11]}{$s_sp})){
	    ## And the score for the current sbject i >= score 
	    ## of existing homolog, add current seq to best hits 
	    ## $a[8]=e-val, $a[7]=score,$a[11]=query,$a[17]=subject
	    if(($scores{$a[11]}{$s_sp}{SCORE}<=$a[7]) && 
	       ($scores{$a[11]}{$s_sp}{EVAL}>=$a[8])){
		my @k=keys(%{$hits{$a[11]}});
		$hits{$a[11]}{$s_sp}{$a[17]}++;
		$scores{$a[11]}{$s_sp}{SCORE}=$a[7];
		$scores{$a[11]}{$s_sp}{EVAL}=$a[8];
	    }
	}
	## If this is the first hit found for this 
	## query prot and this subject species
	else{
	    $hits{$a[11]}{$s_sp}{$a[17]}++;
	    $scores{$a[11]}{$s_sp}{SCORE}=$a[7];
	    $scores{$a[11]}{$s_sp}{EVAL}=$a[8];
	}
    }
}
my @species= ('HUMAN','MOUSE','YEAST','CAEL','DROME');
my %seen;
my $c=0;
foreach my $prot (keys(%hits)){
    my $printed=0;
    next if $seen{$prot};
    $prot=~/.+_(.+?)$/; 
    my $q_sp=$1;  ## query species
    foreach my $species (@species){
	foreach my $hit (keys(%{$hits{$prot}{$species}})){
	    if(defined($hits{$hit}{$q_sp}{$prot})){
		if($printed==0){
		    $c==0 ? 
			print "$prot" :
			print "\n$prot";
		    $printed=1;
		    $c=1;
		}
		print "\t$hit";
		$seen{$hit}++;
	    }
	}
    }
}
print "\n";
