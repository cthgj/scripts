#!/usr/bin/perl -w

###############################################################################
#                           Data Structures				      #
###############################################################################
# 									      #
# $go_ontos{SYNONYM}{$GO:term}= go term synonym	    		              #
# $go_ontos{ONTO}{$GO:term}= go term ontology				      #
# @{$papas{$GO:term}} = @term_ancestor_list				      #
# $prots{GOs}{$PROT}{$GO:term} = will exist if $PROT is annotated to $GO:term #
# $prots{$PROT}{ONTOs}{$onto} = will exist if $PROT is annotated to $onto     #
# @{$prots{$prot}{GOList}}= list of terms a protein is annotated to	      #
# @{$gos{$GO:term}}= list of $PROTS annotated to $GO:term    		      #
# $go_counts{ONTOs}{$o}= num of interactions annotated to $o                  #
# $go_counts{GOs}{$go}{$oo}= num of proteins annotated to $go (onto1) that    #
#                            also have an annotation in onto2                 # 
# $wanted_proteins{$oo}{$PROT}=will exist if this protein is wanted, i.e., if #
#                             it has >=1 DIRECT annotation in EACH of the     #
#                             Desired ontologies.                             #
###############################################################################

use Getopt::Std;
# use Inline (C => 'DATA',
# 	    CC => 'gcc',
# 	    CCFLAGS => '-O3',
#     );
use Inline C;
use IO::File;
use strict;
my %opts;
getopts('eva:d:G:g:o:',\%opts);
my $gaf_file=$opts{a}||die("Need a gene association file (-a)\n");
my $gen_file=$opts{g}||die("Need a genealogy file (-g)\n");
my $subonto=$opts{o}||'PFC';
my $verbose=$opts{v}||0;
my $debug=$opts{d}||0;
my $alt_go_file=$opts{G}||"./data/GO.terms_alt_ids";
my $exp_codes=$opts{e}||undef; ## Use only "good" evidence codes
usage() if $opts{h}; 

#################################
# Get the "good" evidence codes #
#################################
my %good_codes=
(
 "IDA"=>1,
 "IEP"=>1,
 "IGI"=>1,
 "IMP"=>1,
 "IPI"=>1,
 "TAS"=>1
);


###################################
# Check which subontology we want #
###################################
my %want;
$subonto=join("", sort {$b lt $a} split(//,$subonto));
map{$want{$_}++}(split(//,$subonto));
if ($verbose) {
  print STDERR "Calculating $subonto ontologies\n";
} 
#################################
# Declare protein synonyms hash #
#################################
my %synonyms=('LOADED' => 0);

##############################################
# Read Ontology information for each go term #
##############################################
my %go_ontos=get_go_ontologies($alt_go_file);

###############################################
# Read geneology (ancestors) for each go term #
###############################################
my %papas=load_genealogy($gen_file);

###################################
# Parse gene ontology annotations #
###################################
my ($r_prots,$r_gos)=parse_gaf($gaf_file);
our %prots=%{$r_prots};
our %gos=%{$r_gos};

####################
# Collect go stats #
####################
my ($rc1,$rc2)=get_go_stats();
my %go_counts=%{$rc1};
my %wanted_proteins=%{$rc2};

#######################################################
# At this point all information has been collected    #
# except the overlap between each go ($both). So, get #
# this and print				      #
#######################################################
my @GOs=keys(%{$go_counts{GOs}});
c_loop($#GOs,$verbose,\@GOs,\%gos,\%{$go_ontos{ONTO}}, \%{$go_counts{ONTOs}}, \%{$go_counts{GOs}}, \%{$prots{GOs}}, \%wanted_proteins);
print STDERR "Finished correctly\n";
exit;
#####################################################
# This is the algorithm coded in perl for debugging #
#####################################################
# for (my $i=0; $i<scalar(@GOs);$i++) {
#   print STDERR "$i of $#GOs\r" if $verbose;
#   for (my $k=$i+1; $k<scalar(@GOs);$k++) {
#     my @a=sort {$b lt $a} ($GOs[$i],$GOs[$k]);
#     my $pair = join("_", @a);
#     my $oo = join("", sort {$b lt $a} ($go_ontos{ONTO}{$a[0]},$go_ontos{ONTO}{$a[1]}));
#     ## Now go through the prots annotated to each go
#     my $both=0;
#     my %seen;
#     foreach my $prot (@{$gos{$GOs[$i]}} ) {
#       next if $seen{$prot};
#       next unless $wanted_proteins{$oo}{$prot};
#       $seen{$prot}++;
#       if (defined($prots{GOs}{$prot}{$GOs[$k]})) {
# 	$both++;
#       }
#     }
#     $go_counts{GOs}{$a[0]}{$oo}=0 unless defined($go_counts{GOs}{$a[0]}{$oo});
#     $go_counts{GOs}{$a[1]}{$oo}=0 unless defined($go_counts{GOs}{$a[1]}{$oo});
#     $go_counts{ONTOs}{$oo}="0" unless defined($go_counts{ONTOs}{$oo});
#     print "$pair\t$go_counts{ONTOs}{$oo}\t$go_counts{GOs}{$a[0]}{$oo}\t$go_counts{GOs}{$a[1]}{$oo}\t$both\n";
#   }
# }











#####################################################################
#                          SUBROUTINES                              #
#####################################################################

############################################################
## get_go_ontologies                                      ##
############################################################
sub get_go_ontologies{
  print STDERR "Getting ontologies..." if  $verbose;
  my $alt_go_file=shift;
  my %ontos;
  ## get ontologies foreach go
  if ($alt_go_file) {
    open(my $G,"$alt_go_file")||die("Could not open $alt_go_file : $!\n");
    while (<$G>) {
      chomp;
      next if /^\!/;
      /\t([PCF])\t/||die("Cannot parse alt GO file line : $_\n");
      my $oo=$1;
      my @a=split(/\t/);
      $ontos{SYNONYM}{$a[0]}=$a[0];
      $ontos{ONTO}{$a[0]}=$oo;
      ##$a[1] is the alternate?obsolete term names (if any)
      my @b=split(/\s+/,$a[1]);
      map{
	next unless /^GO:\d+$/;
	$ontos{SYNONYM}{$_}=$a[0];
      }@b;
    }
    close($G);
  }
   
  print STDERR "Done\n" if $verbose;
  return(%ontos);
}
############################################################
## load_genealogy                                         ##
############################################################
sub load_genealogy{
  my $file=shift;
  print STDERR "Loading genealogy ($file)..." if  $verbose;
  my %ancestors;
  open(my $A,"$file")||die("Cannot open $file : $!\n");
  while (<$A>) {
    next if /^\#/;
    next unless /\w/;
    chomp;
      
    my @t=split(/\t/);
    my @terms;
    ## make sure we are using the most recent term synonym
    map{
      $go_ontos{SYNONYM}{$_}=$_ unless defined($go_ontos{SYNONYM}{$_});
      push @terms, $go_ontos{SYNONYM}{$_};
    }@t;
    my %seen;
    @terms = grep { ! $seen{$_} ++ } @terms;
    my $child=shift(@terms);
    $ancestors{$child}=[@terms];
  }
  close($A);
  print STDERR "Done\n" if $verbose;
  return(%ancestors);
}
############################################################
## parse GAF                                              ##
############################################################
sub parse_gaf{
  my $gaf_file=shift;
  open(my $GAF,"$gaf_file")||die("cannot open $gaf_file: $!\n");
  my (%prots,%gos,%seen);
  my @ontos=split(//,$subonto);
  while (<$GAF>) {
    next if /^!/;
    print STDERR "Gaf: $.\r" if $verbose;
    chomp;
    my @tmparray=split(/\t/);
    ## $tmparray[1]== name, $tmparray[4] == GO $tmparray[8] == ontology
    my $prot=$tmparray[1];
    my $onto=$tmparray[8];
    next unless defined($want{$onto});
    my $go=$go_ontos{SYNONYM}{$tmparray[4]}||$tmparray[4];    
    next unless $tmparray[3]=~/^$/; ## skip the NOT annotations
    
    #########################################
    # If we only want "good" evidence codes #
    #########################################
    if ($exp_codes) {
	next unless defined($good_codes{$tmparray[6]});
    }

    ## Collect each prot's GOs
    $prots{GOs}{$prot}{$go}++ unless $seen{$prot}{$go};
    push @{$prots{$prot}{GOList}}, $go unless $seen{$prot}{$go};
    
    ## Collect each prot's ontos
    $prots{$prot}{ONTOs}{$onto}++ unless $seen{$prot}{$go};
    ## This collects each GO's prots
    push @{$gos{$go}},$prot unless $seen{$prot}{$go};
    
    ## Now add ancestor GOs
    if (exists($papas{$go})) {
      my @pp=@{$papas{$go}};
      for (my $i=0; $i<=$#pp; $i++) {
	next if $seen{$prot}{$pp[$i]};
	## Some part_of relationships involve different ontologies
	## e.g. GO:0042910(MF) part_of GO:0042908(BP). Do not count
	## these as parents.
	#print STDERR "\$go_ontos{ONTO}{$pp[$i]} eq \$go_ontos{ONTO}{$go} :: $go_ontos{ONTO}{$pp[$i]} eq $go_ontos{ONTO}{$go}\n";
	next unless $go_ontos{ONTO}{$pp[$i]} eq $go_ontos{ONTO}{$go};

	$prots{GOs}{$prot}{$pp[$i]}++;
	#$gos{$pp[$i]}{$prot}++;
	push @{$gos{$pp[$i]}},$prot;
	push @{$prots{$prot}{GOList}}, $pp[$i];
	$seen{$prot}{$pp[$i]}++;
      }
    }
    ## The GAF file can be redundant
    $seen{$prot}{$go}++;
    
  }
  print STDERR "\n" if $verbose;
  return(\%prots,\%gos);
}
############################################################
## get_go_stats                                           ##
############################################################
sub get_go_stats{
  my (%wanted_ps,%count);
  my $NONE=0;
  foreach my $prot (keys(%{$prots{GOs}})) {
    unless ($prots{GOs}{$prot}) {
      $NONE++ ;
      print STDERR "Prot $prot has no annots\n";
    }

    ## Check if this protein has >=1 annotation in each ontology
    my %seen_ontos;
    foreach my $o1 (keys(%want)) {
      foreach my $o2 (keys(%want)) {
	my $oo = join("", sort {$b lt $a} ($o1,$o2));
	##Skip if we have already seen this onto combination (for this prot)
	next if exists($seen_ontos{$oo});
	if (defined($prots{$prot}{ONTOs}{$o1}) &&
	    defined($prots{$prot}{ONTOs}{$o2})) {
	  ############################################
	  # If  $o1 and $o2 are the same, we want    #
	  # this protein only if it has >=2 DIRECT   #
	  # annotations in this ontology	     #
	  ############################################
	  if ($o1 eq $o2) {
		    
	    if ($prots{$prot}{ONTOs}{$o1}>1) {
	      $count{ONTOs}{$oo}++;
	      $wanted_ps{$oo}{$prot}=1;
	      ## If this protein has >=1 DIRECT annotation
	      ## in each ontology, $count{GOs}{$go}{$oo}++;
	      foreach my $go (@{$prots{$prot}{GOList}}) {
		$count{GOs}{$go}{$oo}++;
	      }
	    }
	  }
	  ############################################
	  # If they are different, we want	     #
	  # this protein if it has >=1 DIRECT        #
	  # annotation in each ontology	             #
	  ############################################
	  else {
	    $count{ONTOs}{$oo}++;
	    $wanted_ps{$oo}{$prot}=1;
	    ## If this protein has >=1 DIRECT annotation
	    ## in each ontology, $count{GOs}{$go}{$oo}++;
	    foreach my $go (@{$prots{$prot}{GOList}}) {
	      $count{GOs}{$go}{$oo}++;
	    }
	  }
	}
	$seen_ontos{$oo}=1;
      }
    }
  }
  print STDERR "\n"if $verbose;
  map{print STDERR "$_ : $count{ONTOs}{$_}\n"}keys(%{$count{ONTOs}}) if $verbose;
  print STDERR "NONE : $NONE\n";

  return(\%count,\%wanted_ps);
}

############################################################
## Usage                                                  ##
############################################################

sub usage{
  print STDERR "Not written yet\n"; exit();
    

}
############################################################
## Debug                                                  ##
############################################################

sub debug{
  my $level=shift;
  if ($debug>=$level) {
    print STDOUT "@_";
  }
}

__END__
__C__
int c_loop(int max, int verbose, SV *SVGOlist, SV *SVgos, SV *SVgo_ontos, SV *SVonto_counts, SV *SVgo_counts, SV *SVprot_gos, SV *SVwanted_proteins){
  int counter1,counter2=0, inversed=0, both, cc1,cc2;
  char *go1[12], *go2[12], pair[30], oo[3];

  HV *go_ontos = (HV*)SvRV(SVgo_ontos); // Get a ref to %{$go_ontos{ONTO}}
  AV *GOsSV_ref1 = (AV*)SvRV(SVGOlist); // Get a ref to @GOs
  HV *go_prots=(HV*)SvRV(SVgos); //get a ref to %gos (list of prots per go)
  HV *onto_count=(HV*)SvRV(SVonto_counts); //get a ref to %{$go_counts{ONTOs}}
  HV *go_counts=(HV*)SvRV(SVgo_counts); //get a ref to \%{$go_counts{GOs}}
  HV *prot_gos=(HV*)SvRV(SVprot_gos); //get a ref to \%{$prots{GOs}}
  HV *wanted_prots=(HV*)SvRV(SVwanted_proteins); // get a ref to the %wanted_proteins hash 

  /****************************************************/
  /* This is the main loop that will iterate through  */
  /* the list of GO terms and build ALL possible      */
  /* GO term pairs				      */
  /****************************************************/
  for(counter1=0;counter1<=max; counter1++){ 
   if(verbose==1)
     fprintf(stderr,"%d of %d\r",counter1,max);

    /***********/
    /* Get go1 */
    /***********/
    strcpy(go1,SvPV( *av_fetch(GOsSV_ref1,counter1,0) ,PL_na));

    /************************/
    /* get ontology for go1 */
    /************************/
    SV *go1_onto=SvPV(*hv_fetch( go_ontos , go1 , strlen(go1), 0),PL_na);

    /*********************************************/
    /* Get the list of proteins annotated to go1 */
    /*********************************************/
    SV **go1_prot_list_ref=hv_fetch(go_prots,go1, strlen(go1), 0);
    AV *go1_prot_list;
    if(go1_prot_list_ref != NULL){
      go1_prot_list=(AV*)SvRV(*go1_prot_list_ref);
    }
    /******************************************************************/
    /* Get the $go_counts{GOs}{$go1} hash, keys are onto combinations */
    /******************************************************************/
    HV *go1_count_hash;
    SV **tmp1=hv_fetch( go_counts, go1, strlen(go1), 0);
    if(tmp1!=NULL){
      go1_count_hash=(HV*)SvRV(*tmp1);
    }
    else{fprintf(stderr, "GO %s has no go1_count_hash\n",go1); exit(0);}
    /*************************************/
    /* Start iterating GOlist to get go2 */
    /*************************************/
    for(counter2 = counter1+1; counter2 <=max; counter2++) {        
      
      /******************/
      /* Get go2	*/
      /******************/
      strcpy(go2,SvPV( *av_fetch(GOsSV_ref1,counter2,0) ,PL_na));
      /***************/
      /* Get GO pair */
      /***************/
      if(strcmp(go1,go2)<0){
	strcpy(pair,go1);
	strcat(pair,"_");
	strcat(pair,go2);
	inversed=0; 
      }
      else{
	strcpy(pair,go2);
	strcat(pair,"_");
	strcat(pair,go1);
	inversed=1; /* keeps track of which go is which for printing output */
      }     

      /***********************/
      /* get go2 onto	       */
      /***********************/
      SV *go2_onto=SvPV(*hv_fetch( go_ontos, go2 , strlen(go2), 0),PL_na);
    /***************************************/
      /* Sort to get oo, both ontologies     */
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
      /**************************************************************/
      /* Get the list of wanted prots for this ontology combination */
      /**************************************************************/
      SV **onto_wanted_prots_ref=hv_fetch(wanted_prots,oo, strlen(oo), 0);
      HV *onto_wanted_prots;
      if(onto_wanted_prots_ref != NULL){
	onto_wanted_prots=(HV*)SvRV(*onto_wanted_prots_ref);
      }

      /*********************************************/
      /* Get the list of proteins annotated to go2 */
      /*********************************************/
      SV **go2_prot_list_ref=hv_fetch(go_prots,go2, strlen(go2), 0);
      AV *go2_prot_list;
      if(go2_prot_list_ref != NULL){
	go2_prot_list=(AV*)SvRV(*go2_prot_list_ref);
      }
      /******************************************************************/
      /* Get the $go_counts{GOs}{$go2} hash, keys are onto combinations */
      /******************************************************************/
      HV *go2_count_hash;
      SV **tmp2=hv_fetch( go_counts, go2, strlen(go2), 0);
      if(tmp2!=NULL){
	go2_count_hash=(HV*)SvRV(*tmp2);
      }
      else{fprintf(stderr,"GO %s has no go2_count_hash\n",go2); exit(0);}

      /*********************************************************************/
      /* both counts the number of interactions connecting these two ontos */
      /* Reset to 0 at each iteration.					   */
      /*********************************************************************/
      both=0;
      /**************************************************************/
      /* Count how many of go1's proteins are also annotated to go2 */
      /**************************************************************/
      /*********************************/
      /* Iterate the list of go1 prots */
      /*********************************/
      for(cc1=0;cc1<=av_len(go1_prot_list); cc1++){ 
	      
	/***************************************/
        /* Get the name of the go1 prot	       */
        /***************************************/
	SV **go1_prot_ref=av_fetch(go1_prot_list,cc1,0);
	SV *go1_prot;
	char prot1[30];
	if(go1_prot_ref != NULL){
	    go1_prot=*go1_prot_ref;
	    strcpy(prot1,SvPV(go1_prot,PL_na));
	}
	/*******************************************************/
        /* Check if this protein is wanted in this onto	       */
        /*******************************************************/
	if(onto_wanted_prots_ref != NULL &&
	   hv_exists_ent(onto_wanted_prots,go1_prot,0)==0){
	  continue;
	}
 	/*********************************************************/
        /* Check if this protein is also annotated to go2	 */
        /*********************************************************/
 	SV **prot_gos_ref=hv_fetch(prot_gos,prot1, strlen(prot1), 0);
	HV *prot_gos_hash;
	if(prot_gos_ref != NULL){
	  prot_gos_hash=(HV*)SvRV(*prot_gos_ref); //%{$prots{GOs}{$PROT}}
	    if(hv_exists(prot_gos_hash,go2,strlen(go2))>0){
	      both++;
	    }
	}
	else{
	 
	  continue;
	}
      }  // end list of go1 prots
      /***********************************************************************/
      /* Get the total number of proteins annotated to this onto combination */
      /***********************************************************************/
      int tot_onto=0;
      if(hv_fetch( onto_count, oo , strlen(oo), 0) != NULL){
	tot_onto=SvIV(*hv_fetch( onto_count, oo , strlen(oo), 0));
      }
      /*************************************************************/
      /* For each of the two gos, get the number of prots with >=1 */
      /* DIRECT annotation in that go's onto				   */
      /*************************************************************/
      int go1count=0,go2count=0;
      SV **tmp3=hv_fetch( go1_count_hash, oo , strlen(oo), 0);
      if(tmp3!=NULL)
	go1count=SvIV(*tmp3);
      tmp3=hv_fetch( go2_count_hash, oo , strlen(oo), 0);
      if(tmp3!=NULL)
	go2count=SvIV(*tmp3);
      
      
      inversed==0 ?
      	printf("%s\t%d\t%d\t%d\t%d\n",pair,tot_onto,go1count,go2count,both) :
      	printf("%s\t%d\t%d\t%d\t%d\n",pair,tot_onto,go2count,go1count,both);
    }  //end for counter2
  } //end for counter1
  fflush(stdout);
}

