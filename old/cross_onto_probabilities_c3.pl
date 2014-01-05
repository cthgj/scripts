#!/usr/bin/perl 
## Get the numbers necessary to calculate the probabilities of association of each GO pair
#use strict;
use Getopt::Std;
use Inline C;
use IO::File;
my %opts;
getopts('haOvCio:n:c:s:g:S:d:',\%opts);


my $gaf_file=$ARGV[0]||"gene_association.goa_human";
my $gen_file=$ARGV[1];
my $interactome=$opts{i}||undef; ## Calculate interactome-based probs
my $net_file=$opts{n}||'undef';
my $synfile=$opts{s}||'undef';
my $subonto=$opts{o}||'PFC';
my $over=$opts{O}||undef; ## print OVERrepresented pairs (control)
my $all=$opts{a}||undef;
my $verbose=$opts{v}||undef;
&usage() if $opts{h}; 
#my $calculated=$opts{c}||undef; ## useful when I already have some probabilities calculated 
## and I only want the rest. This option will pass a file
## in which the 1st characters up to the 1st space are the
## gopairs I have. 
my $bdir=$opts{d}||'.'; ## base-dir for output files
my $alt_go_file=$opts{g}||"./data/GO.terms_alt_ids";

## Check which subontology we want
my %want;
$subonto=join("", sort {$b lt $a} split(//,$subonto));
map{$want{$_}++}(split(//,$subonto));
if($verbose){
    print STDERR "Calculating $subonto ontologies\n";
} 
## Many go pairs have exactly the same stats (eg, 12159_0_0_0_0), 
## to avoid doing the same calculation multiple times in R, 
## I condense them so that instead of having a line like:
##      go1_go2\t$tot\$go1\t$go2\t$both etc
## I have a line like:
##     12159_0_0_0_0\t12159\t0\t0\t0 etc 
## and create a file mapping each go to the appropriate number combination
my $condense=$opts{C}||undef; 

my $species=$opts{S}||undef;
die("Need a species name (-S) for a file suffix if condensing (-C) : $species\n") unless $species;




my %interactors;
if($interactome){
    die("Need a synonyms file (-s) if working with a network\n") unless $synfile;
    %interactors=&load_network($net_file); 
}
## Must do this first, as it loads the GO synonyms as well
my %go_ontos=&get_go_ontologies($alt_go_file);
my %papas=&load_genealogy($gen_file);
my %prots=&parse_gaf($gaf_file);
my %go_counts=&count_gos();

# print STDERR "\$go_counts{GOs}{GO:0043231}{CC} : $go_counts{GOs}{'GO:0043231'}{CC}\n";
# print STDERR "\$go_counts{GOs}{GO:0016020}{CC} : $go_counts{GOs}{'GO:0016020'}{CC}\n";

################## Main Program ##########################
## Chose only those gos we want ($subonto)
#my @GOs=keys(%{$go_ontos{ONTO}});


my @ooo=keys(%want);
my $onto1=$ooo[0];
my $onto2=$ooo[1];
if($interactome){
    my %temp_hash;
    foreach my $prot (keys(%{$prots{PROTS}})){
	$prot=&get_name($prot,"NET",2);
	next unless defined($interactors{$prot});
	foreach my $go (keys(%{$prots{PROTS}{$prot}{GOs}})){
	    $temp_hash{$go}++;
	}
    }
    my @GOs=grep{defined($want{$go_ontos{ONTO}{$_}})}keys(%temp_hash);
    print STDERR "GOs: " . scalar(@GOs) . " ($#GOs)\n" if $verbose; 
    my @pps=keys(%interactors);
    # for(my $i=0; $i<$#GOs; $i++){
    # 	for(my $k=$i+1; $k<=$#GOs; $k++){
    # 	    print "$i/$k $GOs[$i]:$GOs[$k]:\t";
    # 	    map{print "$_,"}keys(%{$prots{GOs}{$GOs[$i]}});
    # 	    print "\n\t\t\t\t";
    # 	    map{print "$_,"}keys(%{$prots{GOs}{$GOs[$k]}});
    # 	    print "\n";
    # 	}
    # }
    #c_loop_inter($#GOs,$verbose,$subonto,\%{$go_ontos{ONTO}},\%{$go_counts{GOs}},\%{$prots{GOs}},\%{$go_counts{ONTOs}},\@GOs);
    c_loop_inter2($#GOs,$verbose,$subonto,\%{$go_ontos{ONTO}},\%{$go_counts{GOs}},\%{$prots{GOs}},\%{$go_counts{ONTOs}},\%interactors,\@GOs);
    #c_loop_inter2($#GOs,$verbose,\%interactors);

    
}
else{
    #my @GOs=grep{defined($want{$go_ontos{ONTO}{$_}})}keys(%{$go_ontos{ONTO}});
    my @GOs=grep{defined($want{$go_ontos{ONTO}{$_}})}keys(%{$prots{GOs}});


  
    print STDERR "GOs: " . scalar(@GOs) . "\n" if $verbose; 
    
    foreach my $subgr (keys (%prots)){
	print "subgrGO : $subgr\n";
	foreach my $ss (keys (%{$prots{$subgr}})){
#	    print "ss : $subgr $ss\n";
	    print "\$prots{$subgr}  :: $ss\n";
	    foreach my $val (keys (%{$prots{$subgr}{$ss}})){
		print "\t\$prots{$subgr}{$ss}{$val}  :: $val\n";
		foreach my $subval (keys (%{$prots{$subgr}{$ss}{$val}})){
		    print "\t\t\$prots{$subgr}{$ss}{$val}{$subval} :: $prots{$subgr}{$ss}{$val}{$subval}\n";
		    foreach my $subsubval (keys (%{$prots{$subgr}{$ss}{$val}{$subval}})){
		    print "aaaaa\t\t\$prots{$subgr}{$ss}{$val}{$subval} :: $subsubval\n";
		    }
		}
	
	    }
	}
        }
die();

    print STDERR $prots{GOs}{"GO:AAA"} . " nn \n";
    map{print STDERR "\$prots{GOs}{GO:AAA}{$_} : $prots{GOs}{'GO:AAA'}{$_}\n"}keys(%{$prots{GOs}{"GO:AAA"}});
    print STDERR "==========1================\n";
    map{print STDERR "\$prots{GOs}{GO:BAA}{$_} : $prots{GOs}{'GO:BAA'}{$_}\n"}keys(%{$prots{GOs}{"GO:BAA"}});
    print STDERR "==========================\n";
    print STDERR "\$prots{GOs}{GO:AAA}: $prots{GOs}{'GO:AAA'}\n";
    print STDERR "\$prots{GOs}{GO:BAA}: {$prots{GOs}{'GO:BAA'}\n";
    print STDERR "Ancestors ABB : @{$papas{'GO:ABB'}}\n";
    print STDERR "==========================\n";
    my @lk=keys(%papas);
    for(my $i=0; $i<=$#lk; $i++){print "$i : $lk[$i]\n"}
    map{print "a $_ : @{$papas{$_}}\n"}keys(%papas);
   

 c_loop($#GOs,$verbose,$subonto,\%{$go_ontos{ONTO}},\%{$go_counts{GOs}},\%{$prots{GOs}},\%{$go_counts{ONTOs}},\@GOs);
}


######################## Subroutines ######################

############################################################
## get_go_ontologies                                      ##
############################################################
sub get_go_ontologies{
    print STDERR "Getting ontologies..." if  $verbose;
    my $alt_go_file=shift;
    my %ontos;
    # $ontos{SYNONYM}{'GO:0008150'}='GO:0008150';
    # $ontos{SYNONYM}{'GO:0005575'}='GO:0005575';
    # $ontos{SYNONYM}{'GO:0003674'}='GO:0003674';
    # $ontos{ONTO}{'GO:0008150'}='P';
    # $ontos{ONTO}{'GO:0005575'}='C';
    # $ontos{ONTO}{'GO:0003674'}='F';
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
	    next if $a[0] eq 'GO:0008150' ; ## skip P root
	    next if $a[0] eq 'GO:0005575' ; ## skip C root
	    next if $a[0] eq 'GO:0003674' ; ## skip F root 
	    
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
      next if /^\#/;
      next unless /\w/;
      chomp;
      
      my @t=split(/\t/);
      my @terms;
      ## make sure we are using the most recent term synonym
      $go_ontos{SYNONYM}{$_}=$_ unless defined($go_ontos{SYNONYM}{$_});
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
    print STDERR "Gaf_file : $gaf_file\n";
    my %prots;
    my @ontos=split(//,$subonto);
    while(<GAF>){
	next if /^!/;
	print STDERR "Gaf : $.\r" if $verbose;
	chomp;
	my @tmparray=split(/\t/);
	## $tmparray[1]== name, $tmparray[4] == GO $tmparray[8] == ontology
	my $prot=&get_name($tmparray[1],"NET",0);
	if($interactome){
	    next unless defined($interactors{$prot});
	}
	my $onto=$tmparray[8];
	my $go=$go_ontos{SYNONYM}{$tmparray[4]};
	next unless defined($want{$onto});
	next unless $tmparray[3]=~/^$/;  ## skip the NOT annotations
	## This collects each prot's GOs
	$prots{PROTS}{$prot}{GOs}{$go}++ unless $seen{$prot}{$go};
	print STDOUT "xx \$prots{PROTS}{$prot}{GOs}{$go}\n";
	## This collects each prot's ontos
	$prots{ONTOs}{$prot}{$onto}++ unless $seen{$prot}{$go};
	## This collects each GO's prots
	$prots{GOs}{$go}{$prot}++ unless $seen{$prot}{$go};
	$seen{$prot}{$go}++;
	## Now add ancestor GOs
	my @pp=@{$papas{$go}};
	for(my $i=0; $i<=$#pp; $i++){
	    $prots{PROTS}{$prot}{GOs}{$pp[$i]}++;
	    $prots{GOs}{$pp[$i]}++;
	}
	## For some reason, this adds an extra
	## empty entry to @{$papas{$go}}
	## replaced it with the above loop
	# map{
	#     $prots{PROTS}{$prot}{GOs}{$_}++;
	#     $prots{GOs}{$_}{$prot}++;
	# }@{$papas{$go}};
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
    my $NONE=0;
    my @PROTS;
    if($interactome){@PROTS=keys(%interactors);}
    else{@PROTS=keys(%{$prots{PROTS}})}
    foreach my $prot (@PROTS){
	$prot=&get_name($prot,"NET",2);
	unless ($prots{PROTS}{$prot}{GOs}){
	    $NONE++ ;
	    #print "$prot\txx\n";
	}
	my %this_prot;
	my @os= ('P','F','C');
	for (my $n=0; $n<=$#os;$n++){
	    for (my $k=$n; $k<=$#os;$k++){
		my $oo = join("", sort {$b lt $a} ($os[$n],$os[$k]));
		if ($os[$n] eq $os[$k]){
		    if($interactome){
			## Count prots with >=1 DIRECT annotations 
			## and >=1 INTERACTORS in the current ontologies 
			if ($prots{ONTOs}{$prot}{$os[$n]}>=1){
			    inter:foreach my $inter (keys(%{$interactors{$prot}})){
				if(defined($prots{ONTOs}{$inter}{$os[$n]})){
				    $go_counts{ONTOs}{$oo}++ ;
				    $this_prot{$oo}++;
				    last inter;
				}
			    }
			}
		    }
		    else{
			if($prots{ONTOs}{$prot}{$os[$n]}>1){
			    ## Count prots with >1 DIRECT annotations 
			    ## in the current ontologies
			    $go_counts{ONTOs}{$oo}++ ;
			    ## This protein has at least two
			    ## DIRECT annotations in the current onto
			    $this_prot{$oo}++ ;
			}
		    }
		}
		else{
		    if($interactome){
			## If $prot is annotated to one of the ontologies
			## and has >=1 INTERACTOR annotated to the other
			if ($prots{ONTOs}{$prot}{$os[$n]}>=1){
			  inter:foreach my $inter (keys(%{$interactors{$prot}})){
			      if(defined($prots{ONTOs}{$inter}{$os[$k]})){
				  $go_counts{ONTOs}{$oo}++ ;
				  $this_prot{$oo}++;
				  last inter;
			      }
			  }
			}
			elsif($prots{ONTOs}{$prot}{$os[$k]}>=1){
			  inter:foreach my $inter (keys(%{$interactors{$prot}})){
			      if(defined($prots{ONTOs}{$inter}{$os[$n]})){
				  $go_counts{ONTOs}{$oo}++ ;
				  $this_prot{$oo}++;
				  last inter;
			      }
			  }
			}	
		    }
		    else{
			## Count prots with >=1 DIRECT annotation 
			## in each of the current ontologies
			if(defined($prots{ONTOs}{$prot}{$os[$n]}) && 
			   defined($prots{ONTOs}{$prot}{$os[$k]})){
			    $go_counts{ONTOs}{$oo}++;
			    ## This protein has an annotation
			    ## in each of the current ontos
			    $this_prot{$oo}++;
			}
		    }
		}
	    } ## end for my $k
	} ## end for my $n
	

	## Count occurences of each GO term
	## in EACH ontology/ontology combination
	## iff this prot has at least 2 DIRECT annotations
	## in that ontology/ontology combination
       	foreach my $go (keys(%{$prots{PROTS}{$prot}{GOs}})){
	    ## they should all be their synonym at this point, check just in case
	    die("BAD synonym $go : $go_ontos{SYNONYM}{$go}\n") if $go ne $go_ontos{SYNONYM}{$go};
	    for (my $n=0; $n<=$#os;$n++){
		for (my $k=$n; $k<=$#os;$k++){
		    my $oo = join("", sort {$b lt $a} ($os[$n],$os[$k]));

		    ## If annotations, $this_prot{$oo}>0 iff 
		    ## this prot has at least 2 DIRECT annotations 
		    ## in this $oo. Elsif interactome, iff $prot is 
		    ## annotated to one of the ontologies and has 
		    ## >=1 INTERACTOR annotated to the other
		    $go_counts{GOs}{$go}{$oo}++ if $this_prot{$oo}>0;
		}
	    }
	}
	
    } ## end foreach $prot
    map{print STDERR "$_ : $go_counts{ONTOs}{$_}\n"}keys(%{$go_counts{ONTOs}}) if $verbose;
    print STDERR "NONE : $NONE\n";
    return(%go_counts);
}

############################################################
## Load Network                                           ##
############################################################
sub load_network {
    my %proteins;
    my $network_file=shift;
    open(A, "$net_file")|| die("Cannot open $net_file : $!\n");
    print STDERR "Loading Network...";
    while(<A>){
	next if /^\d+$/;
	next if /^!/;
	next if /^\s*$/;
	chomp;
	my ($bait,$target)=split(/\t/,$_);
	$bait=&get_name($bait,"NET",2);
	$target=&get_name($target,"NET",3);
	## %proteins will hold all the info on all the 
	## network's proteins
	#my $key=join("-",sort {$b lt $a} ($bait,$target));
	$proteins{$bait}{$target}++;
	$proteins{$target}{$bait}++;
	#$proteins{$key}++;
    }
    print STDERR "Done\n";
    close(A);
    return(%proteins)
}
############################################################
## get_name                                               ##
############################################################
sub get_name{
    return($_[0]) unless $synfile; 
    my $name=$_[0];
    my $old_name=$name;
    my $want=$_[1]||"null";
    my $call=$_[2]||undef;

    if ($synonyms{LOADED}==0){
	if(-e $synfile){
	    open(S,$synfile);
	    while(<S>){
		chomp;
		my @a=split(/\t/);
		map{s/\s//g}@a;
		$synonyms{NET}{$a[1]}=$a[0];
		$synonyms{GAF}{$a[1]}=$a[1];
		$synonyms{NET}{$a[0]}=$a[0];
		$synonyms{GAF}{$a[0]}=$a[1];
	    }
	    close(S);
 	$synonyms{LOADED}=1;

	}
    }
    $synonyms{NET}{$name}=$name unless(defined($synonyms{NET}{$name}));
    $synonyms{GAF}{$name}=$name unless(defined($synonyms{GAF}{$name}));
    $want eq 'NET' ? 
	return($synonyms{NET}{$name}) :
	return ($synonyms{GAF}{$name});   
}


############################################################
## Usage                                                  ##
############################################################

sub usage{
    print STDERR "Not written yet\n"; exit();


}



 __END__
__C__
#include <string.h>
#include <stdlib.h>

SV* c_loop_inter2(int max,int verbose, SV* SVontology, SV* go_ontos_ONTO, SV* go_counts_GOs, SV* prots_GOs, SV* go_counts_ONTOs, SV* inter_ref,SV* GOsSV){
   char pair[100], go1[12], go2[12], onto[2], oo[3], wanted_onto[3], id[22];
   int i,k,l,m,both,o,go1count,go2count;

   HV *go_ontos = (HV*)SvRV(go_ontos_ONTO);
   HV *go_counts = (HV*)SvRV(go_counts_GOs);
   int onto_keys = hv_iterinit(go_ontos); // start hash iteration
   int go_counts_keys = hv_iterinit(go_counts);
   HV *prots_GOs_hash=(HV*)SvRV(prots_GOs); // get the %prots{GOS} hash
   HV* inters_hash=(HV*)SvRV(inter_ref);
   HV *go1_count_hash, *go2_count_hash;
   SV* newhashref, *GO2, *GO1,**go1_protss;
  // Get the desired ontology
  strcpy(wanted_onto,SvPV(SVontology,PL_na));

  if (! SvROK(GOsSV))
    croak("GOsSV is not a reference");
  SV *GOsSV_ref1 = (AV*)SvRV(GOsSV); // Get a ref to @GOs

  /****************************************************/
  /* This is the main loop that will iterate through  */
  /* the list of GO terms and build ALL possible      */
  /* GO term pairs				      */
  /****************************************************/
  for(i=0;i<max; i++){
    if(verbose==1)
      fprintf(stderr,"%d of %d\r",i,max);
        //get go1
    SV** GOsSV_ref2=av_fetch(GOsSV_ref1,i,0);
    SV* GO1=*GOsSV_ref2 ;
    strcpy(go1,SvPV(GO1,PL_na));  
    
    // get ontology for go1
    SV* go1_onto=SvPV(*hv_fetch( go_ontos , go1 , strlen(go1), 0),PL_na);
    
    /*******************************************************************/
    /* Check if we want this ontology				       */
    /* 	                                                               */
    /* size_t strspn( char *s1, const char *s2) :: returns the length  */
    /* of the longest substring of s1 that begins at the start of s1   */ 
    /* and consists only of  the characters found in s2.               */
    /*******************************************************************/
    int want=strspn(go1_onto, wanted_onto);
    if(want==0)
      continue;
    
    // Get the %go_count{$go1} hash. It will exist
    // if at least one protein is annotated to go1

    if(hv_exists_ent(go_counts,GO1,0)>0){
      SV* newhashref=*hv_fetch( go_counts, go1, strlen(go1), 0);
      go1_count_hash=(HV*)SvRV(newhashref);
    }
    // Get the number of proteins annotated to go1
    go1_protss=hv_fetch( prots_GOs_hash, go1 , strlen(go1), 0);
    for(k = i+1; k <=max; k++) {
      //get go2
      SV** GOsSV2_ref2=av_fetch(GOsSV_ref1,k,0);
      GO2=*GOsSV2_ref2 ;
      strcpy(go2,SvPV(GO2,PL_na));
      //Get GO pair
      strcpy(pair,go1);
      strcat(pair,"_");
      strcat(pair,go2);
      // get ontology for go2
      SV* go2_onto=SvPV(*hv_fetch( go_ontos , go2 , strlen(go2), 0),PL_na);
      /***************************************/
      /* now sort to get oo, both ontologies */
      /* concatenated in alphabetical order  */
      /***************************************/
      if(strcmp(go1_onto,go2_onto)<0){
      	strcpy(oo,go1_onto);
      	strcat(oo,go2_onto);
      }
      else{
      	strcpy(oo,go2_onto);
      	strcat(oo,go1_onto);
      }
      /* Check if we want this ontology */
      want=strspn(oo, wanted_onto);
      //printf("want : %d, %s,%s\n",want,oo,wanted_onto);
      if(want==0)
	continue;
      // Get the %go_count{$go1} hash. It will exist
      // if at least one protein is annotated to go1
      if(hv_exists_ent(go_counts,GO1,0)>0){
	newhashref=*hv_fetch( go_counts, go1, strlen(go1), 0);
	go1_count_hash=(HV*)SvRV(newhashref);
      }
      /***************************************/
      /* Get the number of proteins anotated */
      /* to each of the go terms	     */
      /***************************************/
      /****************************************************/
      /* The relevant hash entry in %go_counts{GOs}{$go1} */
      /* will only exist if at least one prot (of those   */
      /* that have >1 DIRECT annotation from each of the  */
      /* two current ontologies (oo)) is annotated to go1 */
      /****************************************************/
      
       go1count=0;
       if(hv_exists_ent(go_counts,GO1,0)>0){
	 if(hv_exists(go1_count_hash,oo,strlen(oo))>0){
	   SV* go1cnt=*hv_fetch( go1_count_hash, oo , strlen(oo), 0);
	   go1count=SvIV(go1cnt);
      	}
      }
      go2count=0;
      if(hv_exists_ent(go_counts,GO2,0)>0){
      	newhashref=*hv_fetch(go_counts, go2, strlen(go2), 0);
      	go2_count_hash=(HV*)SvRV(newhashref);
      	if(hv_exists(go2_count_hash,oo,strlen(oo))>0){
      	  SV* go2cnt=*hv_fetch( go2_count_hash, oo , strlen(oo), 0);
      	  go2count=SvIV(go2cnt);
      	}
      }

      both=0; 
      //foreach bait annot to go1
      // (prots_GOs=\%{$prots{GOs}})    $prots{GOs}{$go}{$prot}++
      SV* go1_prots_ref=*hv_fetch(prots_GOs_hash , go1, strlen(go1), 0);
      SV* go2_prots_ref=*hv_fetch(prots_GOs_hash , go2, strlen(go2), 0);
      // go1_prots_hash=$prots{GOs}{$go1}
      HV* go1_prots_hash=(HV*)SvRV(go1_prots_ref);
      HV* go2_prots_hash=(HV*)SvRV(go2_prots_ref);
      int go1_prots=hv_iterinit(go1_prots_hash);
      char bait[20], target[20];
      for(l=0;l<go1_prots;l++){
	HE* bait_entry = hv_iternext(go1_prots_hash);
	// get bait prot name
	strcpy(bait,HePV(bait_entry,PL_na));
	// Now, get bait's interactors
	SV* bait_inters_ref=*hv_fetch(inters_hash, bait, strlen(bait), 0);
	HV* bait_inters_hash=(HV*)SvRV(bait_inters_ref);
	int targets=hv_iterinit(bait_inters_hash);
	for (m = 0; m < targets; m++) {
	  HE* target_entry = hv_iternext(bait_inters_hash);
	  strcpy(target,HePV(target_entry,PL_na));
	  // Check if bait/target pair is annotated to go1 && go2
	  if(hv_exists(go1_prots_hash,bait,strlen(bait))>0 &&
	     hv_exists(go2_prots_hash,target,strlen(target))>0){
	    both++;
	  }
	}
      }
      /******************************************************/
      /* Get the total prots for this ontology combination */
      /******************************************************/
      HV* onto_count=(HV*)SvRV(go_counts_ONTOs);
      int tot_onto=SvIV(*hv_fetch( onto_count, oo , strlen(oo), 0));
      /********************************/
      /* Get over/underrepresentation */
      /********************************/
      int b_both=both;
      int b_tot_onto=tot_onto;
      int b_go2count=go2count;
      int b_go1count=go1count;
      if(both==0)
      	b_both=1;
      if(b_tot_onto==0)
      	b_tot_onto=1;
      if(b_go1count==0)
      	b_go1count=1;
      if(b_go2count==0)
      	b_go2count=1;
      o=b_both*b_tot_onto/b_go1count/b_go2count;
 
      //printf("%s\t%d\t%d\t%d\t%d\t%d\n",pair,tot_onto,go1count,go2count,both,o);
      printf("%s\t%d\t%d\t%d\t%d\t%d\n",pair,tot_onto,go1count,go2count,both,o);
    } // end for k

  } // end for i


   /* int nums=hv_iterinit(inters_hash); */
   /* int ll,kk; */
   /* char bait[20], target[20]; */
   /* for (ll = 0; ll < nums; ll++) { */
   /*   HE* bait_entry = hv_iternext(inters_hash); */
   /*   //SV* key=HePV(bait_entry); */
   /*   strcpy(bait,HePV(bait_entry,PL_na)); */
   /*   //SV* target_hash=HeVAL(bait_entry); */
   /*   HV* target_hash=(HV*)SvRV(HeVAL(bait_entry)); */
   /*   int targets=hv_iterinit(target_hash); */
   /*      for (kk = 0; kk < targets; kk++) { */
   /* 	  HE* target_entry = hv_iternext(target_hash); */
   /* 	  strcpy(target,HePV(target_entry,PL_na)); */
	  

   /* 	  printf("kk : %s %s\n",bait,target); */
   /* 	} */
   /* } */
   fflush(stdout);
 }



SV* c_loop(int max,int verbose, SV* SVontology, SV* hash_ref1, SV* hash_ref2, SV* prots_hash_ref, SV* onto_countref, SV* GOsSV){
  HV* go_ontos, *go_counts;
  HE* onto_entry, *hash_entry2;
  char pair[100], go1[12], go2[12], onto[2], oo[3], wanted_onto[3], id[22];
  int onto_keys, i, k, l, num_keys2,go1count,go2count,both,tot_onto,b_go1count,b_go2count,b_both,b_tot_onto,want;
  float o; //overrepresentation
  SV *go1_onto, *go2_onto, *newhashref, *go1cnt, *go2cnt, **go1_protss , **go2_protss, *go1_prots, *go2_prots, *GOsSV_ref1, **GOsSV_ref2, *GO1, **GOsSV2_ref2, *GO2;
  HV *go1_count_hash, *go2_count_hash,*prots,*onto_count;
  go_ontos = (HV*)SvRV(hash_ref1);
  go_counts = (HV*)SvRV(hash_ref2);
  onto_keys = hv_iterinit(go_ontos); // start hash iteration
  num_keys2 = hv_iterinit(go_counts);
  prots=(HV*)SvRV(prots_hash_ref); // get the %prots{GOS} hash

    printf("go1_go2\t\ttot\t#go1\t#go2\t#both\to\n");


  // Get the desired ontology
  strcpy(wanted_onto,SvPV(SVontology,PL_na));

  if (! SvROK(GOsSV))
    croak("GOsSV is not a reference");
  GOsSV_ref1 = (AV*)SvRV(GOsSV); // Get a ref to @GOs
  
  /****************************************************/
  /* This is the main loop that will iterate through  */
  /* the list of GO terms and build ALL possible      */
  /* GO term pairs				      */
  /****************************************************/

  for(i=0;i<max; i++){
    if(verbose==1)
      fprintf(stderr,"%d of %d\r",i,max);
    
    //get go1
    GOsSV_ref2=av_fetch(GOsSV_ref1,i,0);
    GO1=*GOsSV_ref2 ;
    strcpy(go1,SvPV(GO1,PL_na));
    
    /* if(strcmp(go1,"GO:0043231")!=0 && */
    /*    strcmp(go1,"GO:0016020")!=0){ */
    /*   continue; */
    /* } */
    // get ontology for go1
    go1_onto=SvPV(*hv_fetch( go_ontos , go1 , strlen(go1), 0),PL_na);
    
    /*******************************************************************/
    /* Check if we want this ontology				       */
    /* 	                                                               */
    /* size_t strspn( char *s1, const char *s2) :: returns the length  */
    /* of the longest substring of s1 that begins at the start of s1   */
    /* and consists only of  the characters found in s2.               */
    /*******************************************************************/
    want=strspn(go1_onto, wanted_onto);
    if(want==0)
      continue;
    // Get the %go_count{$go1} hash. It will exist
    // if at least one protein is annotated to go1
    if(hv_exists_ent(go_counts,GO1,0)>0){
      newhashref=*hv_fetch( go_counts, go1, strlen(go1), 0);
      go1_count_hash=(HV*)SvRV(newhashref);
    }
    // Get the number of proteins annotated to go1
    //printf("iiii %d,%d \n",max,i);
    go1_protss=hv_fetch( prots, go1 , strlen(go1), 0);
    for(k = i+1; k <=max; k++) {   
      //get go2
      GOsSV2_ref2=av_fetch(GOsSV_ref1,k,0);
      GO2=*GOsSV2_ref2 ;
      strcpy(go2,SvPV(GO2,PL_na));
      //Get GO pair
      strcpy(pair,go1);
      strcat(pair,"_");
      strcat(pair,go2);
      //printf("A1A GOs %s %s\n",go1,go2);

      // get ontology for go2
      //printf("aaaaaaa max,k,i %d,%d,%d :: %s  \n",max,k,i,go2);
      go2_onto=SvPV(*hv_fetch( go_ontos , go2 , strlen(go2), 0),PL_na);
       
       
      /***************************************/
      /* now sort to get oo, both ontologies */
      /* concatenated in alphabetical order  */
      /***************************************/
      if(strcmp(go1_onto,go2_onto)<0){
      	strcpy(oo,go1_onto);
      	strcat(oo,go2_onto);
      }
      else{
      	strcpy(oo,go2_onto);
      	strcat(oo,go1_onto);
      }
      /* Check if we want this ontology */
      want=strspn(oo, wanted_onto);
      //printf("want : %d, %s,%s\n",want,oo,wanted_onto);
      if(want==0)
	continue;
      //printf("xxx %s\n",oo);
  
      /***************************************/
      /* Get the number of proteins anotated */
      /* to each of the go terms	     */
      /***************************************/
      /****************************************************/
      /* The relevant hash entry in %go_counts{GOs}{$go1} */
      /* will only exist if at least one prot (of those   */
      /* that have >1 DIRECT annotation from each of the  */
      /* two current ontologies (oo)) is annotated to go1 */
      /****************************************************/
      
      go1count=0;
      if(hv_exists_ent(go_counts,GO1,0)>0){
      	if(hv_exists(go1_count_hash,oo,strlen(oo))>0){
      	  go1cnt=*hv_fetch( go1_count_hash, oo , strlen(oo), 0);
      	  go1count=SvIV(go1cnt);
      	}
      }
      go2count=0;
      if(hv_exists_ent(go_counts,GO2,0)>0){
      	newhashref=*hv_fetch(go_counts, go2, strlen(go2), 0);
      	go2_count_hash=(HV*)SvRV(newhashref);
      	if(hv_exists(go2_count_hash,oo,strlen(oo))>0){
      	  go2cnt=*hv_fetch( go2_count_hash, oo , strlen(oo), 0);
      	  go2count=SvIV(go2cnt);
      	}
      }
      // Get the number of proteins annotated to go2
      go2_protss=hv_fetch( prots, go2 , strlen(go2), 0);

      /***************************************************************/
      /* Now, if both hashes exist, that is if at least one protein  */
      /* is annotated to each of the gos, then dereference the SV**  */
      /* to an SV* and count the prots annotated to both gos         */
      /***************************************************************/
      both=0;
      if(go1_protss != NULL  &&
      	 go2_protss != NULL){
	HV* hash1 = (HV*)SvRV(*go1_protss); //$prots{GOs}{$go1}
	HV* hash2 = (HV*)SvRV(*go2_protss); //$prots{GOs}{$go2}
	int num_keys1 = hv_iterinit(hash1); // start hash iteration
	printf("YES GOs %s %s %d\n",go1,go2, num_keys1);
	  /* HE* hash_entry1 = hv_iternext(hash1); */
	  /* printf("x1x hash1 %s \n",SvPV(hv_iterkeysv(hash_entry1), PL_na)); */
	for (l = 0; l < num_keys1; l++) {
	  HE* hash_entry1 = hv_iternext(hash1);
	  //printf("xx hash1 %s l:%d, num1:%d\n",SvPV(hv_iterkeysv(hash_entry1), PL_na), l,num_keys1);
	  num_keys2 = hv_iterinit(hash2); // start hash2 iteration
	  printf("out %d \n",l);
	  int ii;
	  for (ii = 0; ii < num_keys2; ii++) {
	    //printf("in, %d %d %d\n",l,ii,num_keys2);
	    HE* hash_entry2 = hv_iternext(hash2);
	    //printf("xx hash2 %s %s\n",SvPV(hv_iterkeysv(hash_entry1), PL_na),SvPV(hv_iterkeysv(hash_entry2), PL_na));
	    // SvPV(hv_iterkeysv(hash_entry1) == value of current key (prots{GOs}{go1} hash)
	    if(strcmp(SvPV(hv_iterkeysv(hash_entry1), PL_na),SvPV(hv_iterkeysv(hash_entry2), PL_na))==0)
	      both++;
	    //printf("out2, %d %d\n",l,ii);
	  }
	}
	hv_undef (hash1);
	hv_undef (hash2);     
      }
      /******************************************************/
      /* Get the total prots for this ontology combination */
      /******************************************************/
      onto_count=(HV*)SvRV(onto_countref);
      tot_onto=SvIV(*hv_fetch( onto_count, oo , strlen(oo), 0));
      /********************************/
      /* Get over/underrepresentation */
      /********************************/
      b_both=both;
      b_tot_onto=tot_onto;
      b_go2count=go2count;
      b_go1count=go1count;
      if(both==0)
      	b_both=1;
      if(b_tot_onto==0)
      	b_tot_onto=1;
      if(b_go1count==0)
      	b_go1count=1;
      if(b_go2count==0)
      	b_go2count=1;
      o=b_both*b_tot_onto/b_go1count/b_go2count;
      
      /* sprintf(id,"%d_%d_%d_%d",tot_onto,go1count,go2count,both); */
      /* printf("AA %s\n",id); */
      /* exit(0); */

      /****************/
      /* Print output */
      /****************/
      
      printf("%s\t%d\t%d\t%d\t%d\t%d\n",pair,tot_onto,go1count,go2count,both,o);
    } // end for k
  } // end for i
  fflush(stdout);
}


/* /\************************************************************************************\/ */
/* /\*************************** Interactome ********************************************\/ */
/* /\************************************************************************************\/ */
/* SV* c_loop_inter(int max,int verbose, SV* SVontology, SV* hash_ref1, SV* hash_ref2, SV* prots_hash_ref, SV* onto_countref, SV* inter_hash_ref, SV* GOsSV){ */
/*   HV* go_ontos, *go_counts; */
/*   HE* onto_entry, *hash_entry2; */
/*   char pair[100], go1[12], go2[12], onto[2], oo[3], wanted_onto[3], id[22]; */
/*   int onto_keys, i, ii, k, l, num_keys1, num_keys2,go1count,go2count,both,tot_onto,b_go1count,b_go2count,b_both,b_tot_onto,want; */
/*   float o; //overrepresentation */
/*   SV *go1_onto, *go2_onto, *newhashref, *go1cnt, *go2cnt, **go1_protss , **go2_protss, *go1_prots, *go2_prots, *GOsSV_ref1, **GOsSV_ref2, *GO1, **GOsSV2_ref2, *GO2; */
/*   HV *go1_count_hash, *go2_count_hash,*prots,*onto_count, *inters; */
/*   go_ontos = (HV*)SvRV(hash_ref1); */
/*   go_counts = (HV*)SvRV(hash_ref2); */
/*   onto_keys = hv_iterinit(go_ontos); // start hash iteration */
/*   num_keys2 = hv_iterinit(go_counts); */
/*   prots=(HV*)SvRV(prots_hash_ref); // get the %prots{GOS} hash */
/*   inters=(HV*)SvRV(inter_hash_ref); // get the %inters hash (%inters{p1}{p2}) */
/*   // Get the desired ontology */
/*   strcpy(wanted_onto,SvPV(SVontology,PL_na)); */

/*   if (! SvROK(GOsSV)) */
/*     croak("GOsSV is not a reference"); */
/*   GOsSV_ref1 = (AV*)SvRV(GOsSV); // Get a ref to @GOs */
  
/*   /\****************************************************\/ */
/*   /\* This is the main loop that will iterate through  *\/ */
/*   /\* the list of GO terms and build ALL possible      *\/ */
/*   /\* GO term pairs				      *\/ */
/*   /\****************************************************\/ */

/*   for(i=0;i<max; i++){ */
/*     if(verbose==1) */
/*       fprintf(stderr,"%d of %d\r",i,max); */
    
/*     //get go1 */
/*     GOsSV_ref2=av_fetch(GOsSV_ref1,i,0); */
/*     GO1=*GOsSV_ref2 ; */
/*     strcpy(go1,SvPV(GO1,PL_na));   */
    
/*     // get ontology for go1 */
/*     go1_onto=SvPV(*hv_fetch( go_ontos , go1 , strlen(go1), 0),PL_na); */
    
/*     /\*******************************************************************\/ */
/*     /\* Check if we want this ontology				       *\/ */
/*     /\* 	                                                               *\/ */
/*     /\* size_t strspn( char *s1, const char *s2) :: returns the length  *\/ */
/*     /\* of the longest substring of s1 that begins at the start of s1   *\/  */
/*     /\* and consists only of  the characters found in s2.               *\/ */
/*     /\*******************************************************************\/ */
/*     want=strspn(go1_onto, wanted_onto); */
/*     if(want==0) */
/*       continue; */
/*     // Get the %go_count{$go1} hash. It will exist */
/*     // if at least one protein is annotated to go1 */
/*     if(hv_exists_ent(go_counts,GO1,0)>0){ */
/*       newhashref=*hv_fetch( go_counts, go1, strlen(go1), 0); */
/*       go1_count_hash=(HV*)SvRV(newhashref); */
/*     } */
/*     // Get the number of proteins annotated to go1 */
/*     go1_protss=hv_fetch( prots, go1 , strlen(go1), 0); */
/*     for(k = i+1; k <=max; k++) { */
      
/*       //get go2 */
/*       GOsSV2_ref2=av_fetch(GOsSV_ref1,k,0); */
/*       GO2=*GOsSV2_ref2 ; */
/*       strcpy(go2,SvPV(GO2,PL_na)); */
/*       //Get GO pair */
/*       strcpy(pair,go1); */
/*       strcat(pair,"_"); */
/*       strcat(pair,go2); */
/*       // get ontology for go2 */
/*       go2_onto=SvPV(*hv_fetch( go_ontos , go2 , strlen(go1), 0),PL_na); */
       
       
/*       /\***************************************\/ */
/*       /\* now sort to get oo, both ontologies *\/ */
/*       /\* concatenated in alphabetical order  *\/ */
/*       /\***************************************\/ */
/*       if(strcmp(go1_onto,go2_onto)<0){ */
/*       	strcpy(oo,go1_onto); */
/*       	strcat(oo,go2_onto); */
/*       } */
/*       else{ */
/*       	strcpy(oo,go2_onto); */
/*       	strcat(oo,go1_onto); */
/*       } */
/*       /\* Check if we want this ontology *\/ */
/*       want=strspn(oo, wanted_onto); */
/*       //printf("want : %d, %s,%s\n",want,oo,wanted_onto); */
/*       if(want==0) */
/* 	  continue; */
/*       //printf("xxx %s\n",oo);   */
/*       /\***************************************\/ */
/*       /\* Get the number of proteins anotated *\/ */
/*       /\* to each of the go terms	     *\/ */
/*       /\***************************************\/ */
/*       /\****************************************************\/ */
/*       /\* The relevant hash entry in %go_counts{GOs}{$go1} *\/ */
/*       /\* will only exist if at least one prot (of those   *\/ */
/*       /\* that have >1 DIRECT annotation from each of the  *\/ */
/*       /\* two current ontologies (oo)) is annotated to go1 *\/ */
/*       /\****************************************************\/ */
      
/*       go1count=0; */
/*       if(hv_exists_ent(go_counts,GO1,0)>0){ */
/*       	if(hv_exists(go1_count_hash,oo,strlen(oo))>0){ */
/*       	  go1cnt=*hv_fetch( go1_count_hash, oo , strlen(oo), 0); */
/*       	  go1count=SvIV(go1cnt); */
/*       	} */
/*       } */
/*       go2count=0; */
/*       if(hv_exists_ent(go_counts,GO2,0)>0){ */
/*       	newhashref=*hv_fetch(go_counts, go2, strlen(go2), 0); */
/*       	go2_count_hash=(HV*)SvRV(newhashref); */
/*       	if(hv_exists(go2_count_hash,oo,strlen(oo))>0){ */
/*       	  go2cnt=*hv_fetch( go2_count_hash, oo , strlen(oo), 0); */
/*       	  go2count=SvIV(go2cnt); */
/*       	} */
/*       } */
/*       // Get the number of proteins annotated to go2 */
/*       go2_protss=hv_fetch( prots, go2 , strlen(go2), 0); */
/*       /\***************************************************************\/ */
/*       /\* Now, if both hashes exist, that is if at least one protein  *\/ */
/*       /\* is annotated to each of the gos, then dereference the SV**  *\/ */
/*       /\* to an SV* and count the interactions between these gos      *\/ */
/*       /\***************************************************************\/ */
/*       both=0; */
/*       if(go1_protss != NULL  || */
/* 	 go2_protss != NULL){ 	  */
/* 	HV* hash1 = (HV*)SvRV(*go1_protss); //$prots{GOs}{$go1} */
/* 	HV* hash2 = (HV*)SvRV(*go2_protss); //$prots{GOs}{$go2} */
/* 	num_keys1 = hv_iterinit(hash1); // start hash iteration */
/* 	char *prot1[20], *prot2[20]; */
/* 	printf("numkeys1: %d\n",num_keys1); */
/* 	//foreach prot annot to go1 */
/* 	for (l = 0; l < num_keys1; l++) { */
/* 	  printf("AA  (%s-%s)\n",go1,go2); */
	  
/* 	  HE* hash_entry1 = hv_iternext(hash1); */
/* 	  strcpy(prot1,SvPV(hv_iterkeysv(hash_entry1),PL_na)); */
/* 	  num_keys2 = hv_iterinit(hash2); // start hash2 iteration */
/* 	  //foreach prot annot to go2 */
/* 	  for (ii = 0; ii < num_keys2; ii++) { */
/* 	    HE* hash_entry2 = hv_iternext(hash2); */
/* 	    //SvIV(*hv_fetch( onto_count, oo , strlen(oo), 0)); */
/* 	    strcpy(prot2,SvPV(hv_iterkeysv(hash_entry2),PL_na)); */
/* 	    printf("AA %s %s (%s-%s)\n",prot1,prot2,go1,go2); */
	    
/* 	    /\* SV *p1_inters=hv_fetch( inters, prot1, strlen(prot1),0); *\/ */
/* 	    /\* // get %{$interactors{p1}} *\/ */
/* 	    /\* HV* interactors=(HV*)SvRV(p1_inters); *\/ */
/* 	    /\* int num_keys3=hv_iterinit(interactors); *\/ */
/* 	    // SvPV(hv_iterkeysv(hash_entry1) == protname; value of current key (prots{GOs}{go1} hash) */
/* 	  } */
/* 	} */
/* 	hv_undef (hash1); */
/* 	hv_undef (hash2); */
/*       } */
/*       /\******************************************************\/ */
/*       /\* Get the total prots for this ontology combination *\/ */
/*       /\******************************************************\/ */
/*       onto_count=(HV*)SvRV(onto_countref); */
/*       tot_onto=SvIV(*hv_fetch( onto_count, oo , strlen(oo), 0)); */
/*       /\********************************\/ */
/*       /\* Get over/underrepresentation *\/ */
/*       /\********************************\/ */
/*       b_both=both; */
/*       b_tot_onto=tot_onto; */
/*       b_go2count=go2count; */
/*       b_go1count=go1count; */
/*       if(both==0) */
/*       	b_both=1; */
/*       if(b_tot_onto==0) */
/*       	b_tot_onto=1; */
/*       if(b_go1count==0) */
/*       	b_go1count=1; */
/*       if(b_go2count==0) */
/*       	b_go2count=1; */
/*       o=b_both*b_tot_onto/b_go1count/b_go2count; */
      
/*       /\* sprintf(id,"%d_%d_%d_%d",tot_onto,go1count,go2count,both); *\/ */
/*       /\* printf("AA %s\n",id); *\/ */
/*       /\* exit(0); *\/    */

/*       /\****************\/ */
/*       /\* Print output *\/ */
/*       /\****************\/ */
/*       //printf("%s\t%d\t%d\t%d\t%d\t%d\n",pair,tot_onto,go1count,go2count,both,o); */
/*     } // end for k */
/*   } // end for i */
/*   fflush(stdout); */
/* } */
 
 


























 /* SV* c_interactome(int max, SV* hash_ref1, SV* prots){ */
/*   int i,k; */
/*   char pair[100], go1[12], go2[12], prot1[12],prot2[12]; */
/*   SV  *prots_ref; */
  
/*   HV *interactors=(HV*)SvRV(hash_ref1); */
/*   // Get a ref to @pps */
/*   SV* proteins=(AV*)SvRV(prots); */
/*   for(i=0;i<max; i++){ */
/*     SV **protein_pair = av_fetch(proteins,i,0); */
/*     char *p_pair=SvPV(*protein_pair,PL_na); */
    
/*     /\***************************************************************************\/ */
/*     /\* char *strtok(char *s1, const char *s2)				       *\/ */
/*     /\* repeated calls to this function break string s1 into "tokens"--that     *\/ */
/*     /\* is the string is broken into substrings, each terminating with a '\0',  *\/ */
/*     /\* where the '\0' replaces any characters contained in string s2. The      *\/ */
/*     /\* first call uses the string to be tokenized as s1; subsequent calls use  *\/ */
/*     /\* NULL as the first argument. A pointer to the beginning of the current   *\/ */
/*     /\* token is returned; NULL is returned if there are no more tokens.	       *\/ */
/*     /\***************************************************************************\/ */
/*     char *prot1=strtok(p_pair,"-"); */
/*     char *prot2=strtok(NULL,"-"); */

/*     /\*******************************************************\/ */
/*     /\* Now, go through the gos of each protein and get the *\/ */
/*     /\* number that annotate both of them		   *\/ */
/*     /\*******************************************************\/ */
/*     //$prots{PROTS}{$prot}{GOs}{$go} */
/*   } */
/*   printf("YES\n"); */

/*   fflush(stdout); */
/* } */
