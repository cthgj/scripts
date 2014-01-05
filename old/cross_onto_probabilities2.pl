#!/usr/bin/perl 
## Get the numbers necessary to calculate the probabilities of association of each GO pair

use strict;
use Getopt::Std;
my %opts;
getopts('Aaovc:s:g:',\%opts);


my $gaf_file=$ARGV[0]||"gene_association.goa_human";
my $obo_file=$ARGV[1]||"biological_process.obo";
my $gen_file=$ARGV[2];
my $subonto=$opts{s}||'P';
my $over=$opts{o}||undef; ## print OVERrepresented pairs (control)
my $all=$opts{a}||undef;
my $verbose=$opts{v}||undef;
my $All_ontos=$opts{A}||undef; ## calculate cross ontology probabilities
my $calculated=$opts{c}||undef; ## useful when I already have some probabilities calculated 
                                ## and I only want the rest. This option will pass a file
                                ## in which the 1st characters up to the 1st space are the
                                ## gopairs I have. 
my $alt_go_file=$opts{g}||'./data/GO.terms_alt_ids';

## Must do this first, as it loads the GO synonyms as well
my %go_ontos=&get_go_ontologies($alt_go_file);

my %papas=&load_genealogy($gen_file);
my %prots=&parse_gaf($gaf_file);
my %go_counts=&count_gos();


################## Main Program ##########################
my @GOs=keys(%{$go_ontos{ONTO}});
print STDERR "GOs: " . scalar(@GOs) . "\n"; 
#Now, calculate all numbers
my %have;
for (my $n=0; $n<=$#GOs;$n++){
    die() if $n==10;
    my $go1=$go_ontos{SYNONYM}{$GOs[$n]};
    print STDERR "$n of $#GOs\r";
    #unless($go1 eq 'GO:0009987' || $go1 eq 'GO:2001001'){next}
    next if $go1 eq 'GO:0008150' ; ## skip P root
    next if $go1 eq 'GO:0005575' ; ## skip C root
    next if $go1 eq 'GO:0003674' ; ## skip F root 
    for (my $k=$n+1; $k<scalar(@GOs);$k++){
	my $go2=$go_ontos{SYNONYM}{$GOs[$k]};
	#unless($go2 eq 'GO:0009987' || $go2 eq 'GO:2001001'){next}
	next if $go2 eq 'GO:0008150' ;
	next if $go2 eq 'GO:0005575' ;
	next if $go2 eq 'GO:0003674' ; 
	my @b=sort {$b lt $a} ($go1,$go2);
	my $pair = join("_",@b);
	$go1=$b[0];
	$go2=$b[1];
	my $oo = join("", sort {$b lt $a} ($go_ontos{ONTO}{$go1},$go_ontos{ONTO}{$go2}));
	# my $onto1=$go_ontos{ONTO}{$go1}||die("No ontology for term : $go1\n");
	# my $onto2=$go_ontos{ONTO}{$go2}||die("No ontology for term : $go2\n");
	# my @c=sort {$b lt $a} ($onto1,$onto2);
	# my $oo=$c[0] . $c[1];
	my $both=0;
	foreach my $prot (keys(%{$prots{GOs}{$go1}})){
	    next unless defined($prots{ONTOs}{$prot}{$go_ontos{ONTO}{$go1}});
	    next unless defined($prots{ONTOs}{$prot}{$go_ontos{ONTO}{$go2}});
	    $both++ if defined($prots{GOs}{$go2}{$prot});
	}
	my $a=$go_counts{GOs}{$go1}{$oo};
	my $b=$go_counts{GOs}{$go2}{$oo};
	$a=1 if $a==0;
	$b=1 if $b==0;
	my $o=$both*$go_counts{ONTOs}{$oo}/$a/$b;
	$go_counts{GOs}{$go1}{$oo}=0 unless $go_counts{GOs}{$go1}{$oo};
	$go_counts{GOs}{$go2}{$oo}=0 unless $go_counts{GOs}{$go2}{$oo};

	#####################################################################################
        # This is the string that I will pass to R. R will then do		            # 
	# phyper(x,m,n,k), where						            #
	#       x=the number of white balls drawn       = #both			            #
	#       m=the number of white balls in the urn  = #go1			            #
	#       n=the number of black balls in the urn  = #Total-#go1		            #         
	#       k=the number of balls drawn from the urn= #go1+#go2-#both                   #
        #####################################################################################
	#                                                                                   #
	#####################################################################################
        # $pair=GO pair									    #
	# $both=overlap						                            #
	# $go_counts{ONTOs}{$oo}= num of prots with >=1 annot in each ontology              #
	# $go_counts{GOs}{$go1}{$oo} = num of prots with >=1 annot in each ontology 	    #
	#                              that are NOT annotated to go1			    #
	# $go_counts{GOs}{$go2}{$oo} = num of prots with >=1 annot in each ontology 	    #
	#                              that ARE annotated to go2			    #
        #####################################################################################
                				     

	my $string="$pair\t$go_counts{ONTOs}{$oo}\t$go_counts{GOs}{$go1}{$oo}\t$go_counts{GOs}{$go2}{$oo}\t$both\t$o\n";

	if($over){##if we want OVERrepresented pairs (control)
	    $o >=1 && do {
		print $string;
	    }
	}
	## if we want all, over and underrepresented
	elsif($all){
	    print $string;	  
	}
	##if we want UNDERrepresented pairs
	else{
	    $o <1 && do {
		print $string;
	    }
	}
	
    }
}

######################## Subroutines ######################

############################################################
## get_go_ontologies                                      ##
############################################################
sub get_go_ontologies{
    print STDERR "Getting ontologies..." if  $verbose;
    my $alt_go_file=shift;
    my %ontos;
    ## get ontologies foreach go
    if($alt_go_file){
	open(G,"$alt_go_file")||die("Could not open $alt_go_file : $!\n");
	while(<G>){
	    chomp;
	    next if /^\!/;
	    /\t([PCF])\t/||die("Cannot parse alt GO file line : $_\n");
	    my $oo=$1;
	    my @a=split(/\t/);
	    ## $a[0] is the current term
	    $ontos{SYNONYM}{$a[0]}=$a[0];
	    $ontos{ONTO}{$a[0]}=$oo;
	    ##$a[1] is the alternate?obsolete term names (if any)
	    my @b=split(/\s+/,$a[1]);
	    map{
		next unless /^GO:\d+$/;
		$ontos{SYNONYM}{$_}=$a[0];
	    }@b;
	}
    }
    print STDERR "Done\n" if $verbose;
    return(%ontos);
}

############################################################
## load_genealogy                                         ##
############################################################
sub load_genealogy{
    print STDERR "Loading genealogy..." if  $verbose;
    my $file=shift;
    my %ancestors;
    open(A,"$file")||die("Cannot open $file : $!\n");
    while(<A>){
	chomp;
	my @t=split(/\t/);
	my @terms;
	## make sure we are using the most recent term synonym
	map{push @terms, $go_ontos{SYNONYM}{$_}}@t;
	my %seen;
	@terms = grep { ! $seen{$_} ++ } @terms;
	my $child=shift(@terms);
	$ancestors{$child}=[@terms];
    }
    print STDERR "Done\n" if $verbose;
    return(%ancestors);
}

############################################################
## parse GAF                                              ##
############################################################
sub parse_gaf{
    my $gaf_file=shift;
    my %seen;
    open(GAF,"$gaf_file")||die("cannot open $gaf_file: $!\n");
    my %prots;
    while(<GAF>){
	next if /^!/;
	print STDERR "Gaf : $.\r" if $verbose;
	chomp;
	my @tmparray=split(/\t/);
	## $tmparray[1]== name, $tmparray[4] == GO $tmparray[8] == ontology
	my $prot=$tmparray[1];
	my $onto=$tmparray[8];
	my $go=$go_ontos{SYNONYM}{$tmparray[4]};
	unless($All_ontos){next unless $onto eq $subonto;}
	next unless $tmparray[3]=~/^$/;  ## skip the NOT annotations
	## This collects each prot's GOs
	$prots{PROTS}{$prot}{GOs}{$go}++ unless $seen{$prot}{$go};
	## This collects each prot's ontos
	$prots{ONTOs}{$prot}{$onto}++ unless $seen{$prot}{$go};
	## This collects each GO's prots
	$prots{GOs}{$go}{$prot}++ unless $seen{$prot}{$go};
	$seen{$prot}{$go}++;
	## Now add ancestor GOs
	map{
	    $prots{PROTS}{$prot}{GOs}{$_}++;
	    $prots{GOs}{$_}{$prot}++;
	}@{$papas{$go}};
    }
    close(GAF);
    print STDERR "\n" if $verbose;
    return(%prots);
}

############################################################
## count_gos                                              ##
############################################################
sub count_gos{
    print STDERR "Counting GOs...\n" if $verbose;
    my %go_counts;
    foreach my $prot (keys(%{$prots{PROTS}})){
	my %this_prot;
	my @os= ('P','F','C');
	for (my $n=0; $n<=$#os;$n++){
	    for (my $k=$n; $k<=$#os;$k++){
		my @c=sort {$b lt $a} ($os[$n],$os[$k]);
		my $oo=$c[0] . $c[1];
		if ($os[$n] eq $os[$k]){
		    ## Count prots with >1 DIRECT annotations 
		    ## in the current ontology
		    $go_counts{ONTOs}{$oo}++ if $prots{ONTOs}{$prot}{$os[$n]}>1;
		    ## This protein has at least two
		    ## DIRECT annotations in the current onto
		    $this_prot{$oo}++ if $prots{ONTOs}{$prot}{$os[$n]}>1;
		}
		## Count prots with >=1 DIRECT annotation 
		## in each of the current ontologies
		else{
		    if(defined($prots{ONTOs}{$prot}{$os[$n]}) && 
		       defined($prots{ONTOs}{$prot}{$os[$k]})){
			$go_counts{ONTOs}{$oo}++;
			## This protein has an annotation
			## in each of the current ontos
			$this_prot{$oo}++;
		    }
		   
		}
	    }
	}
	    

	## Count occurences of each GO term
	## in EACH ontology/ontology combination
	## iff this prot has at least 2 DIRECT annotations
	## in that ontology/ontology combination
       	foreach my $go (keys(%{$prots{PROTS}{$prot}{GOs}})){
	    ## they should all be their synonym at this point, check just in case
	    die("BAD synonym $go : $go_ontos{SYNONYM}{$go}\n") if $go ne $go_ontos{SYNONYM}{$go};
	    for (my $n=0; $n<=$#os;$n++){
		for (my $k=$n; $k<=$#os;$k++){
		    my @c=sort {$b lt $a} ($os[$n],$os[$k]);
		    my $oo=$c[0] . $c[1];  
		    #print "$prot \$go_counts{GOs}{$go}{$oo}++\n";
		    $go_counts{GOs}{$go}{$oo}++ if $this_prot{$oo}>0;
		    #print "$go ($go_ontos{ONTO}{$go}): $prot : $oo : $this_prot{$oo}\n";
		}
	    }
		    #print "-----------------------\n";

	    # foreach my $onto (keys(%this_prot)){
	    # 	print "$go ($go_ontos{ONTO}{$go}): $prot : $onto : $this_prot{$onto}\n";
	    # 	$go_counts{GOs}{$go}{$onto}++ if $this_prot{$onto}>1;
	    # }
	    # print "-----------------------\n";
	}
	
    } ## end foreach $prot
    map{print STDERR "$_ : $go_counts{ONTOs}{$_}\n"}keys(%{$go_counts{ONTOs}}) if $verbose;
   
    return(%go_counts);
}


