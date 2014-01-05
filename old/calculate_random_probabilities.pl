#!/usr/bin/perl 
## Get the numbers necessary to calculate the probabilities of association of each GO pair

use strict;
use Getopt::Std;
my %opts;
getopts('vs:',\%opts);
# use GO::Basic;
# use GO::Parser;
# use GO::Model::Association;
my $gaf_file=$ARGV[0]||"gene_association.goa_human";
my $obo_file=$ARGV[1]||"biological_process.obo";
my $gen_file=$ARGV[2];
my $subonto=$opts{s}||'P';


my %papas=&load_geneology($gen_file);
## parse GAF
open(GAF,"$gaf_file")||die("cannot open $gaf_file: $!\n");
my (%gos, %count,%gocount);
   while(<GAF>){
	next if /^!/;
	print STDERR "Gaf : $.\r";
	chomp;
	my @tmparray=split(/\t/);
	next unless $tmparray[8] eq $subonto;
	next unless $tmparray[3]=~/^$/;  ## skip the NOT annotation
	$gos{$tmparray[4]}{$tmparray[1]}++;## $tmparray[1]== name, $tmparray[4] == GO
	$count{$tmparray[1]}{$tmparray[4]}++ ;

}
close(GAF);
my (%seen,%go_count,%good_genes);
##count prots with >1 dir annot
foreach my $prot (keys(%count)){
    my @gos=keys(%{$count{$prot}});
    if (scalar(@gos)>1){ ## we are only interested in prots with >1 dir annot
	$good_genes{$prot}=1 ;
	foreach my $ggo (@gos){
	    $go_count{$ggo}++ unless $seen{$_}{$prot};
	    $seen{$ggo}{$prot}++;
	    $gos{$ggo}{$prot}++;
	    map{
		$go_count{$_}++ unless $seen{$_}{$prot};
		$gos{$_}{$prot}++;
		$seen{$_}{$prot}++;
	    }@{$papas{$ggo}}
	}
    }
    else{next}
}
my $pcount=scalar(keys(%good_genes));
print STDERR "\nPROT : $pcount\n";
my @GOs=keys(%go_count);
#Now, calculate all numbers
    for (my $n=0; $n<=$#GOs;$n++){
	print STDERR "$n of $#GOs\r";
	next if $GOs[$n] eq 'GO:0008150' ;
	next if $GOs[$n] eq 'GO:0005575' ;
	for (my $k=$n+1; $k<scalar(@GOs);$k++){
	    next if $GOs[$k] eq 'GO:0008150' ;
	    next if $GOs[$k] eq 'GO:0005575' ;
	    my $both=0;

	    ######
	    ## Skip this, we will choose random overlaps
	    ######
	    # foreach my $prot (keys(%{$gos{$GOs[$n]}})){
	    # 	next unless $good_genes{$prot};
	    # 	$both++ if defined $gos{$GOs[$k]}{$prot}
	    # }
#$o1 = $n*$Ngg/$nocc{$K[0]}/$nocc{$K[1]};
 	   
	    $go_count{$GOs[$n]}>$go_count{$GOs[$k]} ? 
	    ($both=int(rand($go_count{$GOs[$n]}))) : 
	    ($both=int(rand($go_count{$GOs[$k]}))) ; 
	    

	    my @b=sort {$b lt $a} ($GOs[$n],$GOs[$k]);
	    my $pair = join("_",@b);
	    my $o=$both*$pcount/$go_count{$GOs[$n]}/$go_count{$GOs[$k]};
	    $o <1 && do {
		print "$pair\t$pcount\t$go_count{$GOs[$n]}\t$go_count{$GOs[$k]}\t$both\t$o\n";
	    }

	}

}

#Now, calculate all numbers
#foreach my $go (keys(%gos)){
    ## Calculate number of prots annotated
    ## to each GO (implicit and direct)
#    my $term = $obo->get_term($go);   # fetch a term by ID
# #    my $ancestor_terms = $obo->get_parent_terms($term->acc);
#     $gocount{$go}=scalar(keys(%{$gos{$go}}));
    
#     print STDOUT "GO: $go\n";
#     foreach my $anc_term (@$ancestor_terms) {
# 	print STDOUT "AA : $anc_term\n";
# 	print $anc_term->acc, " ";
#     }
#     print "\n";
#}


print STDERR "\n\n-------------------------------------------\n\n";


sub load_geneology{
    my $file=shift;
    my %ancestors;
    open(A,"$file")||die("Cannot open $file : $!\n");
    while(<A>){
	chomp;
	my @terms=split(/\t/);
	my %seen;
	@terms = grep { ! $seen{$_} ++ } @terms;
	my $child=shift(@terms);
	$ancestors{$child}=[@terms];
    }
    print SDTERR "GEN : $.\r";

    return(%ancestors);
}
