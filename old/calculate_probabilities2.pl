#!/usr/bin/perl 
## Get the numbers necessary to calculate the probabilities of association of each GO pair

use strict;
use Getopt::Std;
my %opts;
getopts('Aaovc:s:g:',\%opts);
# use GO::Basic;
# use GO::Parser;
# use GO::Model::Association;
my $gaf_file=$ARGV[0]||"gene_association.goa_human";
my $obo_file=$ARGV[1]||"biological_process.obo";
my $gen_file=$ARGV[2];
my $subonto=$opts{s}||'P';
my $over=$opts{o}||undef; ## print OVERrepresented pairs (control)
my $all=$opts{a}||undef;
my $All_ontos=$opts{A}||undef; ## calculate cross ontology probabilities
my $calculated=$opts{c}||undef; ## useful when I already have some probabilities calculated 
                                ## and I only want the rest. This option will pass a file
                                ## in which the 1st characters up to the 1st space are the
                                ## gopairs I have. 
my %cross_onto;  
my $alt_go_file=$opts{g}||'./data/GO.terms_alt_ids';
my %ontos;
my %have=&read_the_ones_I_have($calculated) if $calculated;
my %papas=&load_geneology($gen_file);
my %pcounts;
my %synonyms;
## get ontologies foreach go
if($alt_go_file){
    open(G,"$alt_go_file")||die("Could not open $alt_go_file : $!\n");
    while(<G>){
	chomp;
	next if /^\!/;
	/\t([PCF])\t/;
	my $oo=$1;
	my @a=split(/\t/);
	## $a[0] is the current term
	$ontos{$a[0]}=$oo;
	$synonyms{$a[0]}=$a[0];
	##$a[1] is the alternate?obsolete term names (if any)
	my @b=split(/\s+/,$a[1]);
	map{
	    next unless /^GO:\d+$/;
	    $ontos{$_}=$oo;
	    $synonyms{$_}=$a[0];
	}@b;
    }
}



## parse GAF
open(GAF,"$gaf_file")||die("cannot open $gaf_file: $!\n");
my (%gos, %count,%gocount);
   while(<GAF>){
	next if /^!/;
	print STDERR "Gaf : $.\r";
	chomp;
	my @tmparray=split(/\t/);
	unless($All_ontos){next unless $tmparray[8] eq $subonto;}
	next unless $tmparray[3]=~/^$/;  ## skip the NOT annotation
        ## If I am working with all ontologies, the
	## background in each case must be the num of 
	## prots that have at least one direct annotation
	## from EACH ontology
	if($All_ontos){
	    ## $tmparray[1]== name, $tmparray[4] == GO $tmparray[8] == ontology
	    # print "\$gos{$tmparray[4]}{$tmparray[8]}{$tmparray[1]}++;\n";
	    $gos{$tmparray[4]}{$tmparray[8]}{$tmparray[1]}++;
	    $count{$tmparray[1]}{$tmparray[8]}{$tmparray[4]}++ ;
	}
	else{
	    $gos{$tmparray[4]}{$tmparray[1]}++;## $tmparray[1]== name, $tmparray[4] == GO
	    $count{$tmparray[1]}{$tmparray[4]}++ ;
	}
	
	
}
close(GAF);
my (%seen,%go_count,%good_genes);
## populates pcount with the counts
## for each pair of ontos
&get_pcounts();
print "\nPF : $pcounts{PF}\n";
print "CF : $pcounts{CF}\n";
print "PC : $pcounts{PC}\n";
print "PP : $pcounts{PP}\n";
print "CC : $pcounts{CC}\n";
print "FF : $pcounts{FF}\n";
die();
##count prots with >1 dir annot
foreach my $prot (keys(%count)){
    ## If I am working with all ontologies, the
    ## background in each case must be the num of 
    ## prots that have at least one direct annotation
    ## from EACH ontology
    my (@gosP,@gosF,@gosC);
    if($All_ontos){
	foreach my $o ('P','F','C'){
	    my @gos=keys(%{$count{$prot}{$o}});
	    if($#gos>=0 ){
		$good_genes{$prot}{$o}=1;
		foreach my $ggo (@gos){
		    $ggo=$synonyms{$ggo};
		    unless ($seen{$ggo}{$prot}){
			## If I am working with all ontologies, I need, foreach 
			## gopair of combined ontology OO (eg, MP or CP), the # of prots 
			## annotated to each go that have at least one annotation in the 
			## ontology of the other. eg: go1(MF)_go2(BP) need number of prots
			## with >=1 annot BP that are annotated to go1
			
			$go_count{$ggo}{'P'}++ if $count{$prot}{'P'}{$ggo}>0;
			$go_count{$ggo}{'F'}++ if $count{$prot}{'F'}{$ggo}>0;
			$go_count{$ggo}{'C'}++ if $count{$prot}{'C'}{$ggo}>0;
			$seen{$ggo}{$prot}++;

			map{
			    $go_count{$_}{'P'}++ if $count{$prot}{'P'}{$_}>0;
			    $go_count{$_}{'F'}++ if $count{$prot}{'F'}{$_}>0;
			    $go_count{$_}{'C'}++ if $count{$prot}{'C'}{$_}>0;
			    $gos{$_}{$o}{$prot}++;
			    $seen{$_}{$prot}++;
			}@{$papas{$ggo}};
		    }
		    

		}
	    }
	    else{next}
	}
    } ## end if allontos
    else{
	my @gos=keys(%{$count{$prot}});
	if (scalar(@gos)>1){ ## we are only interested in prots with >1 dir annot
	    $good_genes{$prot}=1 ;
	    foreach my $ggo (@gos){
		$ggo=$synonyms{$ggo};
		$go_count{$ggo}++ unless $seen{$ggo}{$prot};
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
}## end foreach $prot


my $pcount=scalar(keys(%good_genes));
print STDERR "\nPROT : $pcount\n";
my @GOs=keys(%go_count);
print STDERR "GOs: " . scalar(@GOs) . "\n"; 
#Now, calculate all numbers
for (my $n=0; $n<=$#GOs;$n++){
    my $go1=$synonyms{$GOs[$n]};
    print STDERR "$n of $#GOs\r";
    next if $go1 eq 'GO:0008150' ; ## skip P root
    next if $go1 eq 'GO:0005575' ; ## skip C root
    next if $go1 eq 'GO:0003674' ; ## skip F root
    for (my $k=$n+1; $k<scalar(@GOs);$k++){
	my $go2=$synonyms{$GOs[$k]};
	next if $go2 eq 'GO:0008150' ;
	next if $go2 eq 'GO:0005575' ;
	next if $go2 eq 'GO:0003674' ; 
	my @b=sort {$b lt $a} ($go1,$go2);
	my $pair = join("_",@b);
	my $go1=$b[0];
	my $go2=$b[1];
	next if defined($have{$pair});
	#next unless $pair eq 'GO:0009987_GO:0044449';
	
	my $both=0;
	if($All_ontos){
	    my $onto1=$ontos{$go1}||die("No ontology for term : $go1\n");
	    my $onto2=$ontos{$go2}||die("No ontology for term : $go2\n");
	    foreach my $prot (keys(%{$gos{$go1}{$onto1}})){
		next unless $good_genes{$prot}{$onto1};
		next unless $good_genes{$prot}{$onto2};
		$both++ if defined $gos{$go2}{$onto2}{$prot};
	    }

	    my $oo=$onto1 . $onto2;
	    my $go1_count=$go_count{$go1}{$onto2}||0;
	    my $go2_count=$go_count{$go2}{$onto1}||0;
	    #print "aa $oo : $pcounts{$oo}\n"; die();
	    next unless $pair eq 'GO:0010181_GO:0040010';
	    print "$pair : my \$o=$both*$pcounts{$oo}/$go1_count/$go2_count\n";
	    die();
	    my $o=$both*$pcounts{$oo}/$go_count{$go1}{$onto2}/$go_count{$go2}{$onto1};
#	    die("aa $go1,$go2 :: $pair\t$pcounts{$oo}\t$go1_count\t$go2_count\t$both\t$o\n") unless $go1_count;
	    if($over){##if we want OVERrepresented pairs (control)
		$o >=1 && do {
		    print "$pair\t$pcounts{$oo}\t$go_count{$go1}\t$go_count{$go2}\t$both\t$o\n";
		}
	    }
	    ## if we want all, over and underrepresented
	    elsif($all){
#		print "$pair\t$pcounts{$oo}\t$go_count{$go1}\t$go_count{$go2}\t$both\t$o\n";
		print "$pair\t$pcounts{$oo}\t$go1_count\t$go2_count\t$both\t$o\n";
	    }
	    else{
		$o <1 && do {
		    print "$pair\t$pcounts{$oo}\t$go_count{$go1}\t$go_count{$go2}\t$both\t$o\n";
		}
	    }
	}
	else{
	    foreach my $prot (keys(%{$gos{$go1}})){
		next unless $good_genes{$prot};
		$both++ if defined $gos{$go2}{$prot}
	    }
#$o1 = $n*$Ngg/$nocc{$K[0]}/$nocc{$K[1]};
	    my $o=$both*$pcount/$go_count{$go1}/$go_count{$go2};
	    if($over){##if we want OVERrepresented pairs (control)
		$o >=1 && do {
		    print "$pair\t$pcount\t$go_count{$go1}\t$go_count{$go2}\t$both\t$o\n";
		}
	    }
	    ## if we want all, over and underrepresented
	    elsif($all){
		print "$pair\t$pcount\t$go_count{$go1}\t$go_count{$go2}\t$both\t$o\n";
	    }
	    else{
		$o <1 && do {
		    print "$pair\t$pcount\t$go_count{$go1}\t$go_count{$go2}\t$both\t$o\n";
		}
	    }
	}
    } ## end for k $#GOs
} ## end for n $#GOs
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

sub read_the_ones_I_have{
    my %hash;
    my $file=shift;
    print STDERR "Loading have...";
    open(C,"$file")||die("Could not open $file : $!\n");
    while(<C>){
	/^(.+?)\s/;
	$hash{$1}++;
    }
    close(C);
    print STDERR "Done\n";
    return(%hash);
}

sub get_pcounts(){
    
    foreach my $prot (keys(%count)){
	# print "pp : $prot : $count{$prot}{P} : $count{$prot}{F} : $count{$prot}{C}\n";
	if(defined($count{$prot}{'P'}) && defined($count{$prot}{'F'})){
	    $pcounts{'FP'}++;
	    $pcounts{'PF'}=$pcounts{'FP'};
	}
	if(defined($count{$prot}{'P'}) && defined($count{$prot}{'C'})){
	    $pcounts{'CP'}++;
	    $pcounts{'PC'}=$pcounts{'CP'};
	}
	if(defined($count{$prot}{'C'}) && defined($count{$prot}{'F'})){
	    $pcounts{'FC'}++;
	    $pcounts{'CF'}=$pcounts{'FC'};
	}
	if(scalar(keys(%{$count{$prot}{'P'}}))>1){$pcounts{'PP'}++}
	if(scalar(keys(%{$count{$prot}{'F'}}))>1){$pcounts{'FF'}++}
	if(scalar(keys(%{$count{$prot}{'C'}}))>1){$pcounts{'CC'}++}
    }
}
