#!/usr/bin/perl -w
###############################################################################
#                           Data Structures				      #
###############################################################################
# 									      #
# $go_ontos{SYNONYM}{$GO:term}= go term synonym	    		              #
# $go_ontos{ONTO}{$GO:term}= go term ontology				      #
# @{$papas{$GO:term}} = @term_ancestor_list				      #
# $prots{$PROT}{GOs}{$GO:term} = will exist if $PROT is annotated to $GO:term #
# $prots{$PROT}{ONTOs}{$onto} = will exist if $PROT is annotated to $onto     #
# $wanted_prots{$oo}{$PROT} = will exist if $PROT has >=1 interaction         #
#                             connecting the two ontos of $oo                 #
# @{$prots{$prot}{GOList}}= list of terms a protein is annotated to	      #
# @{$gos{$GO:term}}= list of $PROTS annotated to $GO:term    		      #
# $interactors{$a}{$b}= exists if $a inter/ts with $b OR vice versa (sorted). #
# $inter_lists{$a}=[@list]= list of prots that interact with $a               #
# $go_counts{ONTOs}{$o}= num of interactions annotated to $o                  #
# $go_counts{GOs}{$go}{$oo}= num of interactions between any prot annotatated #
#                            to $go (onto1) and any annotated to onto2        #
###############################################################################

use Getopt::Std;
use Inline (C => 'DATA',
	    CC => 'gcc',
	    CCFLAGS => '-O3',
    );
#use Inline C;
use IO::File;
use strict;
my %opts;
getopts('pva:D:G:g:m:n:o:s:',\%opts);
my $gaf_file=$opts{a}||die("Need a gene association file (-a)\n");
my $gen_file=$opts{g}||die("Need a genealogy file (-g)\n");
my $net_file=$opts{n}||die("Need a network file (-n)\n");
my $synfile=$opts{m}||die("Need a map file (-m)\n");
my $subonto=$opts{o}||'PFC';
my $verbose=$opts{v}||0;
my $debug=$opts{D}||0;
my $alt_go_file=$opts{G}||"./data/GO.terms_alt_ids";
my $species=$opts{s}||undef;
usage() if $opts{h}; 

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


###############################
# Load and parse Network file #
###############################
my ($r1,$r2,$r3,$r4)=load_network($net_file,1);
our %go_counts=%{$r1};
our %interactions=%{$r2};
our %inter_lists=%{$r3};
our %wanted_prots=%{$r4};
#######################################################
# At this point all information has been collected    #
# except the overlap between each go ($both). So, get #
# this and print				      #
#######################################################
my @GOs=keys(%{$go_counts{GOs}});

#####################################################
# This is the algorithm coded in perl for debugging #
#####################################################
if($opts{p}){
    for (my $i=0; $i<scalar(@GOs);$i++) { 
	my $go1=$GOs[$i];
	print STDERR "$i of $#GOs\r" if $verbose;
	for (my $k=$i+1; $k<scalar(@GOs);$k++) {
	    my $go2=$GOs[$k];
	    
	    my @a=sort {$b lt $a} ($go1,$go2);
	    my $pair = join("_", @a);
	    my $oo = join("", sort {$b lt $a} ($go_ontos{ONTO}{$a[0]},$go_ontos{ONTO}{$a[1]}));
	    ## Now go through the prots annotated to each go
	    my ($both,$both2, $same)=(0,0,0);
	    my %seen2;
	    my @go1_proteins=@{$gos{$go1}};
	    foreach my $go1prot (@go1_proteins) {
		next unless defined($wanted_prots{$oo}{$go1prot});
		foreach my $go2prot (@{$inter_lists{$go1prot}}){
		    next unless defined($wanted_prots{$oo}{$go2prot});
		    print "$go1prot : $go2prot\n";
		    $same++ if $go1prot eq $go2prot;
		    my @pp=sort{$b lt $a} ($go1prot,$go2prot);
		    ## if go2prot is annotated to go2
		    if(defined ($prots{GOs}{$go2prot}{$go2})){
			next if $seen2{$pp[0] . $pp[1]};
			$both++;
			$seen2{$pp[0] . $pp[1]}++;
		    }
		 
		}
	    }
	    
	    $go_counts{GOs}{$a[0]}{$oo}=0 unless defined($go_counts{GOs}{$a[0]}{$oo});
	    $go_counts{GOs}{$a[1]}{$oo}=0 unless defined($go_counts{GOs}{$a[1]}{$oo});
	    print "$pair\t$go_counts{ONTOs}{$oo}\t$go_counts{GOs}{$a[0]}{$oo}\t$go_counts{GOs}{$a[1]}{$oo}\t$both\t$same\n";
	}
    }
    
    exit;
    
}
else{
    c_loop($#GOs,$verbose,\@GOs,\%gos,\%wanted_prots,\%interactions,\%inter_lists,\%{$go_ontos{ONTO}}, \%{$go_counts{ONTOs}}, \%{$go_counts{GOs}},\%{$prots{GOs}});
    print STDERR "Finished correctly\n";
}

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
  ## Quickly read throug the network file to
  ## collect a list of wanted proteins
  my $href=load_network($net_file,0);
  my %prots_in_network=%{$href};
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
    my $prot=get_name($tmparray[1],"NET",0);
    ##Skip prots that are not in the network
    next unless defined($prots_in_network{$prot});
    my $onto=$tmparray[8];
    my $go=$go_ontos{SYNONYM}{$tmparray[4]}||$tmparray[4];
    next unless defined($want{$onto});
    next unless $tmparray[3]=~/^$/; ## skip the NOT annotations
    ## Collect each prot's GOs
    $prots{GOs}{$prot}{$go}++ unless $seen{$prot}{$go};
    push @{$prots{$prot}{GOList}}, $go unless $seen{$prot}{$go};
    ## Collect each prot's ontos
    $prots{$prot}{ONTOs}{$onto}++ unless $seen{$prot}{$go};
    ## This collects each GO's prots
    push @{$gos{$go}},$prot unless $seen{$prot}{$go};
    debug(1,"G:$go ($prot): @{$gos{$go}}\n");
    ## Now add ancestor GOs
    if (exists($papas{$go})) {
      my @pp=@{$papas{$go}};
      for (my $i=0; $i<=$#pp; $i++) {
	  next if $seen{$prot}{$pp[$i]};
	  ## Some part_of relationships involve different ontologies
	  ## e.g. GO:0042910(MF) part_of GO:0042908(BP). Do not count
	  ## these as parents.


	  $prots{GOs}{$prot}{$pp[$i]}++;
	  #$gos{$pp[$i]}{$prot}++;
	  push @{$gos{$pp[$i]}},$prot ;
	  push @{$prots{$prot}{GOList}}, $pp[$i] ;
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
## Load Network                                           ##
############################################################
sub load_network {
  my $network_file=shift;
  ## Mode 0 will return a simple hash $h{prot_name} with
  ## all prots in the network. this is then passed to 
  ## parse_gaf(), allowing it to skip any unneeded prots
  ##
  ## Mode 1 will read the network in detail and return
  ## three hashes: %go_counts,%interactions and %wanted_prots
  
  my $mode=shift; 
  my (%interactors,,%inter_lists,%go_stats,%w, %h);
  my @wanted_ontos=keys(%want);
  my $NONE=0;
  open(my $A, "$net_file")|| die("Cannot open $net_file : $!\n");
  while (<$A>) {
      
    next if /^\d+$/;
    next if /^!/;
    next if /^\s$/;
    chomp;
    my ($bait,$target)=split(/\t/,$_);
    $bait=get_name($bait,"NET",2);
    $target=get_name($target,"NET",3);
    if ($mode==0) {
      $h{$bait}=$h{$target}=1;
    } 
    else {
      print STDERR "Net: $.\r" if $verbose;
      ## Sort the bait/target pair, this way
      ## even if it also present as a target/bait
      ## pair, it will be counted correctly
      my @aa=sort{$b lt $a} get_name($bait,"NET",2),get_name($target,"NET",3);
      $bait=$aa[0];
      $target=$aa[1];
      $interactors{$bait}{$target}++;
      push @{$inter_lists{$target}},$bait;
      unless($bait eq $target){push @{$inter_lists{$bait}},$target;}
      
 
      if ($bait eq $target) {
	  next if $interactors{$bait}{$target}>2;
	  $NONE++ unless defined($prots{GOs}{$bait});
      } 
      else {
	  next if $interactors{$bait}{$target}>1;
	  $NONE++ unless defined($prots{GOs}{$target});
	  $NONE++ unless defined($prots{GOs}{$bait});
      }
      ## Check which ontologies the bait/target
      ## pair is annotated to
      my %seen;
      # foreach my $o1 (keys(%want)) {
      # 	foreach my $o2 (keys(%want)) {
      my @ontos=keys(%want);
      for (my $i=0; $i<=$#ontos; $i++){
	  my $o1=$ontos[$i];
	  for (my $k=$i; $k<=$#ontos; $k++){
	      my $o2=$ontos[$k];
	      my $oo = join("", sort {$b lt $a} ($o1,$o2));
	      if ( (defined($prots{$bait}{ONTOs}{$o1}) && 
		    (defined($prots{$target}{ONTOs}{$o2})) ) ||
		   (defined($prots{$bait}{ONTOs}{$o2}) &&
		    (defined($prots{$target}{ONTOs}{$o1})) ) ) {
		  $go_stats{ONTOs}{$oo}++ unless $seen{$oo};
		  
		  ## This interaction bridges the two ontos.
		  ## Therefore, we want both the bait and the 
		  ## target for this onto
		  $w{$oo}{$bait}=1;
		  $w{$oo}{$target}=1;
		  $seen{$oo}++;
	      }
	  }
      } 
      ## Get a list of unique GOs annotating
      ## this interaction
      my %k;
      my @redundant_gos=(keys(%{$prots{GOs}{$bait}}),keys(%{$prots{GOs}{$target}}));
      map{$k{$_}++}@redundant_gos;
      my @ok=keys(%k);
      debug(1,"$bait,$target, gos:@ok");
      foreach my $oo (keys(%seen)) {     
	  map{$go_stats{GOs}{$_}{$oo}++;}keys(%k);	  
      }
    }
  }
  close($A);
  if($mode>0 && $verbose){
      print STDERR "\nNONE: $NONE\n";
      map{print STDERR "$_ : $go_stats{ONTOs}{$_}\n"}keys(%{$go_stats{ONTOs}}) if $verbose;
  }
    $mode==0 ? return(\%h) : return(\%go_stats,\%interactors,\%inter_lists,\%w);

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
    
  #    unless(defined($synonyms{NET}{$name})){
  #defined($synonyms{GAF}{$name})){
  # if($name eq 'all'){
  if ($synonyms{LOADED}==0) {
    if (-e $synfile) {
      open(my $S,$synfile);
      while (<$S>) {
	chomp;
	my @a=split(/\t/);
	map{s/\s//g}@a;
	$synonyms{NET}{$a[1]}=$a[0];
	$synonyms{GAF}{$a[1]}=$a[1];
	$synonyms{NET}{$a[0]}=$a[0];
	$synonyms{GAF}{$a[0]}=$a[1];
      }
      close($S);
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
#include "uthash.h"
  struct my_hash {
  int id;                   
  char name[40];             /* key */
  UT_hash_handle hh;         /* makes this structure hashable */
};

int c_loop(int max, int verbose, SV *SVGOlist, SV *SVgos, SV *SVwanted_prots, SV *SVinteractions, SV *SVinter_lists, SV *SVgo_ontos, SV *SVonto_counts, SV *SVgo_counts, SV *SVprot_gos){
  int counter1,counter2=0, inversed=0, both, same, cc1,cc2;
  char *go1[12], *go2[12], pair[30], oo[3];
  
  HV *go_ontos = (HV*)SvRV(SVgo_ontos); // Get a ref to %{$go_ontos{ONTO}}
  AV *GOsSV_ref1 = (AV*)SvRV(SVGOlist); // Get a ref to @GOs
  HV *wanted_prots=(HV*)SvRV(SVwanted_prots); //get a ref to %wanted_prots
  HV *go_prots=(HV*)SvRV(SVgos); //get a ref to %gos (list of prots per go)
  HV *inters=(HV*)SvRV(SVinteractions); //get a ref to %interactions
  HV *inters_lists=(HV*)SvRV(SVinter_lists); //get a ref to %inter_lists
  HV *onto_count=(HV*)SvRV(SVonto_counts); //get a ref to %{$go_counts{ONTOs}}
  HV *go_counts=(HV*)SvRV(SVgo_counts); //get a ref to \%{$go_counts{GOs}}
  HV *prot_gos=(HV*)SvRV(SVprot_gos); //get a ref to \%{$prots{GOs}}

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
    for(counter2 = counter1; counter2 <=max; counter2++) {        
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
	inversed=1; //keeps track of which go is which for printing output
      }     

      /***********************/
      /* get go2 onto	     */
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

      /**************************************************************/
      /* When go1 and go2 are the same, all their prots will        */
      /* be seen twice. Make a hash to keep track		    */
      /**************************************************************/
      struct my_hash *s,*seen_prot_pairs = NULL;
      int hash_counter=0;

      /*********************************************************************/
      /* both counts the number of interactions connecting these two ontos */
      /* Reset to 0 at each iteration.					   */
      /*********************************************************************/
      both=0;

      /********************************************************/
      /* same counts the number of self interactions. 	      */
      /* Reset to 0 at each iteration.			      */
      /********************************************************/
      same=0;

      /*********************************/
      /* Iterate the list of go1 prots */
      /*********************************/
      for(cc1=0;cc1<=av_len(go1_prot_list); cc1++){ 
	/***************************************/
        /* Get the name of the go1 prot	       */
        /***************************************/
	SV **go1_prot_ref=av_fetch(go1_prot_list,cc1,0);
	SV *go1_prot;
	if(go1_prot_ref != NULL){
	    go1_prot=*go1_prot_ref;
	}
	
	/*******************************************************/
        /* Check if this protein is wanted in this onto	       */
        /*******************************************************/
	if(onto_wanted_prots_ref != NULL &&
	   hv_exists_ent(onto_wanted_prots,go1_prot,0)==0){
	  continue;
	}

	/*******************************************************/
	/* Get the list of proteins interacting with go1_prot  */
	/*******************************************************/
	SV **go1_prot_inter_ref=hv_fetch(inters_lists,SvPV(go1_prot,PL_na),strlen(SvPV(go1_prot,PL_na)),0);
	AV *go1_prot_inter_list;
	if(go1_prot_inter_ref != NULL){
	  go1_prot_inter_list =(AV*)SvRV(*go1_prot_inter_ref );
	}

	/*****************************************************/
        /* Iterate the list of go1_prot's interactors	     */
        /*****************************************************/
	for(cc2=0;cc2<=av_len(go1_prot_inter_list); cc2++){ 
	  SV **inter_ref=av_fetch(go1_prot_inter_list,cc2,0);
	  SV *inter_prot;
	  if(inter_ref != NULL){
	    inter_prot=*inter_ref;
	  }
	  /***********************************************************/
          /* Check if this protein is wanted in this ontology	     */
          /***********************************************************/
	  if(hv_exists_ent(onto_wanted_prots,inter_prot,0)==0)
	    continue;
	  
	  /********************************************************************/
          /* Get a reference to the hash containing the annotations	      */
	  /* of inter_prot (%{$prots{GOs}{$inter_prot}}). Keys are gos	      */
          /********************************************************************/
	  SV **prot_gos_ref=hv_fetch(prot_gos, SvPV(inter_prot,PL_na), strlen( SvPV(inter_prot,PL_na)), 0);
	  HV *prot_gos_hash;
	  prot_gos_hash=(HV*)SvRV(*prot_gos_ref); //%{$prots{GOs}{$PROT}}

	  /********************************************/
          /* If inter_prot is annotated to go2	      */
          /********************************************/
	  if(hv_exists(prot_gos_hash,go2,strlen(go2))!=0){
	    /***********************************************************************/
	    /* Sort the two proteins to get a unique hash id and avoid counting	 */
	    /* the same interaction (a-b and b-a) twice for the same go pair.	 */
	    /***********************************************************************/
	    char prot1[20], prot2[20];
	    if(strcmp(SvPV(go1_prot,PL_na), SvPV(inter_prot,PL_na))<0){
	      strcpy(prot1,SvPV(go1_prot,PL_na));
	      strcpy(prot2,SvPV(inter_prot,PL_na));
	    }
	    else{
	      strcpy(prot2,SvPV(go1_prot,PL_na));
	      strcpy(prot1,SvPV(inter_prot,PL_na));
	    }
	    char hash_id[60];
	    strcpy(hash_id,prot1);
	    strcat(hash_id,prot2);
	    
	    /************************************************************************/
            /* Check if this prot pair has been seen before, for these gos.	    */
	    /* If it has skip it, and if not hasn't, add the prot pair to the hash. */
            /************************************************************************/
	    HASH_FIND_STR( seen_prot_pairs, hash_id, s);
	    if (s){
	      continue;
	    }
	    else{
	      //add to hash
	      s = malloc(sizeof(struct my_hash));
	      strcpy(s->name, hash_id);
	      s->id = hash_counter++;
	      HASH_ADD_STR(seen_prot_pairs, name, s );
	    }
	    
	    /*********************************************************/
            /* Keep track of the number of self interactions.	     */
	    /* That way, they can be removed if we so desire.	     */
            /*********************************************************/
	    if(strcmp(prot1,prot2)==0)
	      same++;
	    printf("%s\t%s\n",prot1,prot2);
	    both++;
	  } // end if(hv_exists(prot_gos_hash,go2,strlen(go2))!=0)
	} // end for cc2
      }// end list of go1 prots */

    
      /*****************************************************************/
      /* Get the total number of interactions in this onto combination */
      /*****************************************************************/
      int tot_onto=0;
      if(hv_fetch( onto_count, oo , strlen(oo), 0) != NULL){
	tot_onto=SvIV(*hv_fetch( onto_count, oo , strlen(oo), 0));
      }
      /**********************************************************/
      /* Get the number of interactions annotated to the gopair */
      /* that are between prots annotated to both ontos of oo   */
      /**********************************************************/      
      int go1count=0, go2count=0;
      SV **tmp3=hv_fetch( go1_count_hash, oo , strlen(oo), 0);
      if(tmp3!=NULL)
	go1count=SvIV(*tmp3);
      tmp3=hv_fetch( go2_count_hash, oo , strlen(oo), 0);
      if(tmp3!=NULL)
	go2count=SvIV(*tmp3);
      
      inversed==0 ?
      	printf("%s\t%d\t%d\t%d\t%d\t%d\n",pair,tot_onto,go1count,go2count,both,same) :
      	printf("%s\t%d\t%d\t%d\t%d\t%d\n",pair,tot_onto,go2count,go1count,both,same);
      
      /*************************************/
      /* Free the memory held by the hash  */
      /*************************************/
      /* iterate over hash elements, deleting and freeing them */
      struct my_hash *current_id, *tmp;
      HASH_ITER(hh, seen_prot_pairs, current_id, tmp) {
	HASH_DEL(seen_prot_pairs,current_id);  /* delete; seen_prot_pairs advances to next */
	free(current_id);
      }
      
      
    } //end for counter2
 
  } //end for counter1
  fflush(stdout);
  
}
