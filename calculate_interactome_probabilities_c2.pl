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
# $interactions_o{$oo}{bait+target}= exists if this interaction  bridges $oo. #
# $interactions_g{$go}{$oo}{bait+target}= exists if this interaction  bridges #
#                                         $oo and is annotated to $go.        #
# $interactions_gl{$go}{$oo}= @inters that this bait/target is annotated to.  #
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
getopts('pva:d:G:g:m:n:o:s:',\%opts);
my $gaf_file=$opts{a}||die("Need a gene association file (-a)\n");
my $gen_file=$opts{g}||die("Need a genealogy file (-g)\n");
my $net_file=$opts{n}||die("Need a network file (-n)\n");
my $synfile=$opts{m}||die("Need a map file (-m)\n");
my $subonto=$opts{o}||'PFC';
my $verbose=$opts{v}||0;
my $debug=$opts{d}||0;
my $alt_go_file=$opts{G}||"./data/GO.terms_alt_ids";
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
my %prots=%{$r_prots};
my %gos=%{$r_gos};


###############################
# Load and parse Network file #
###############################
my ($r1,$r2,$r3)=load_network($net_file,1);
my %interactions_o=%{$r1};
my %interactions_g=%{$r2};
my %interactions_gl=%{$r3};
#######################################################
# At this point all information has been collected    #
# except the overlap between each go ($both). So, get #
# this and print				      #
#######################################################
my @GOs=keys(%interactions_g);

################################################################
# Get the number of interactions for each ontology combination #
################################################################
my %onto_count;
foreach my $oo (keys(%interactions_o)){
    $onto_count{$oo}=scalar(keys(%{$interactions_o{$oo}}))||0;
}

#####################################################
# This is the algorithm coded in perl for debugging #
#####################################################
if($opts{p}){
    for (my $i=0; $i<scalar(@GOs);$i++) { 
	my $go1=$GOs[$i];
	print STDERR "$i of $#GOs\r" if $verbose;
	for (my $k=$i; $k<scalar(@GOs);$k++) {
	    
	    my @a=sort {$b lt $a} ($GOs[$i],$GOs[$k]);
	    my $go1=$a[0];
	    my $go2=$a[1];
	    my $pair = join("_", @a);
	    my $oo = join("", sort {$b lt $a} ($go_ontos{ONTO}{$a[0]},$go_ontos{ONTO}{$a[1]}));
	    ## Now go through the prots annotated to each go
	    my ($both,$both2, $same)=(0,0,0);
	    my $tot_onto=$onto_count{$oo}||0;
	    my $go1cnt=scalar(@{$interactions_g{$go1}{$oo}})||0;
	    my $go2cnt=scalar(@{$interactions_g{$go2}{$oo}})||0;
	    if($go1cnt>0){
		my @inters=@{$interactions_gl{$go1}{$oo}};
		foreach my $interaction(@inters){
		  $both++ if defined($interactions_g{$go2}{$oo}{$interaction});
		  @a=split(/_/, $interaction);
		  $same++ if $a[0] eq $a[1];
		}
	    }
	    





	    print "$pair\t$tot_onto\t$go1cnt\t$go2cnt\t$both\t$same\n";
	}
    }    
    exit;
}
else{
  c_loop($#GOs,$verbose,\@GOs,\%interactions_g,\%interactions_o,\%interactions_gl,\%onto_count,\%{$go_ontos{ONTO}});
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
  my (%interactors_o,%interactors_g, %interactors_gl,%h,%i);
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

      my $hash_id=$bait . "+" . $target;
      next if exists($i{$hash_id});
      $i{$hash_id}++;
      

      if ($bait eq $target) {
	  $NONE++ unless defined($prots{GOs}{$bait});
      } 
      else {
	  $NONE++ unless defined($prots{GOs}{$target});
	  $NONE++ unless defined($prots{GOs}{$bait});
      }

      my %seen;
      my @ontos=keys(%want);
  
       
      ## Check which ontology combinations are
      ## bridged by this interaction
      my @redundant_gos=(keys(%{$prots{GOs}{$bait}}),keys(%{$prots{GOs}{$target}}));
      my %k;
      map{$k{$_}++}@redundant_gos;
      my @inter_gos=keys(%k);    
      for (my $i=0; $i<=$#ontos; $i++){
	  my $o1=$ontos[$i];
	  for (my $k=$i; $k<=$#ontos; $k++){
	      my $o2=$ontos[$k];
	      my $oo = join("", sort {$b lt $a} ($o1,$o2));
	      next if $seen{$oo};
	      ## If this interaction bridges $oo
	      if ( (defined($prots{$bait}{ONTOs}{$o1}) && 
		    (defined($prots{$target}{ONTOs}{$o2})) ) ||
		   (defined($prots{$bait}{ONTOs}{$o2}) &&
		    (defined($prots{$target}{ONTOs}{$o1})) ) ) {
		  $interactors_o{$oo}{$hash_id}++;
		  
		  foreach my $go (@inter_gos){
		      $go=$go_ontos{SYNONYM}{$go};
		      $interactors_g{$go}{$oo}{$hash_id}++;
		      ########################################################
		      # There is a bug in Inline.c that causes memory leaks  #
		      # when iterating through a hash. So, get the LIST of   #
		      # each interaction's GOs			             #
		      ########################################################
		      push @{$interactors_gl{$go}{$oo}},$hash_id;
		  }
		  $seen{$oo}++;
	      }
	  }
      }
    }
  }

  close($A);
  if($mode>0 && $verbose){
      print STDERR "\nNONE: $NONE\n";
      map{print STDERR "$_ : " . scalar(keys(%{$interactors_o{$_}})) . "\n"}keys(%interactors_o) if $verbose;
  }
  $mode==0 ? return(\%h) : return(\%interactors_o, \%interactors_g, \%interactors_gl);

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
//c_loop($#GOs,        $verbose,      \@GOs,   \%interactions_g,  \%interactions_o, \%interactions_gl,   \%onto_count, \%{$go_ontos{ONTO}});
void c_loop(int max, int verbose, SV *SVGOlist, SV *SVinter_gos, SV *SVinter_ontos, SV *SVinter_golist,  SV *SVonto_count,SV *SVgo_ontos){
 int counter1,counter2, inversed=0, both, same;
 char *go1, *go2, pair[30], oo[3];  
  
  HV *go_ontos = (HV*)SvRV(SVgo_ontos); // Get a ref to %{$go_ontos{ONTO}}
  AV *GOsSV_ref1 = (AV*)SvRV(SVGOlist); // Get a ref to @GOs
  HV *inter_gos=(HV*)SvRV(SVinter_gos); //get a ref to %interactions_g
  HV *inter_ontos=(HV*)SvRV(SVinter_ontos); //get a ref to %interactions_o
  HV *inter_golist=(HV*)SvRV(SVinter_golist); //get a ref to %interactions_gl
  HV *onto_count=(HV*)SvRV(SVonto_count); //get a ref to %onto_counts

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
    go1=SvPV( *av_fetch(GOsSV_ref1,counter1,0) ,PL_na);
    /************************/
    /* get ontology for go1 */
    /************************/
    char *go1_onto=SvPV(*hv_fetch( go_ontos , go1 , strlen(go1), 0),PL_na);

    /*************************************/
    /* Start iterating GOlist to get go2 */
    /*************************************/
    for(counter2 = counter1; counter2 <=max; counter2++) {
      /******************/
      /* Get go2	*/
      /******************/
      go2=SvPV( *av_fetch(GOsSV_ref1,counter2,0) ,PL_na);

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
      char *go2_onto=SvPV(*hv_fetch( go_ontos, go2 , strlen(go2), 0),PL_na);
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

      /********************************************/
      /* Get these go's $interactions_g{$go} hash */
      /********************************************/
      SV **hashref1,**hashref2;
      HV *go1_inter_hash, *go1_oo_inters, *go2_inter_hash, *go2_oo_inters;
      hashref1 = hv_fetch(inter_gos, go1, strlen(go1),0);
      if(hashref1 != NULL){
	go1_inter_hash=(HV*)SvRV(*hashref1);//%{$interactions_g{$go1}}
	hashref2 = hv_fetch(go1_inter_hash, oo, strlen(oo),0);
	/*****************************************************/
        /* Get the $interactions_g{$go1}{$oo} hash, keys     */
        /* are interacting protein pairs that bridge the     */
        /* two ontologies.				     */
        /*****************************************************/
	if(hashref2 != NULL){
	  go1_oo_inters=(HV*)SvRV(*hashref2);
	}
      }
      SV **hashref4;
      SV **hashref3 = hv_fetch(inter_gos, go2, strlen(go2),0);
      if(hashref3 != NULL){
	go2_inter_hash=(HV*)SvRV(*hashref3);
	hashref4 = hv_fetch(go2_inter_hash, oo, strlen(oo),0);
	if(hashref4 != NULL){
	  go2_oo_inters=(HV*)SvRV(*hashref4);
	}
      }

      /*************************************************/
      /* Get the list of interactions annotated to go1 */
      /*************************************************/
      SV **go1_inter_list_ref=hv_fetch(inter_golist, go1, strlen(go1),0);
      AV *go1_inter_list;
      int go1count=0;
      if(go1_inter_list_ref!= NULL){
	HV *go1_inter_hash=(HV*)SvRV(*go1_inter_list_ref);
	SV **go1_oo_inter_list_ref=hv_fetch(go1_inter_hash,oo,strlen(oo),0);
	if(go1_oo_inter_list_ref!= NULL){
	  go1_inter_list=(AV*)SvRV(*go1_oo_inter_list_ref);
	  go1count=av_len(go1_inter_list)+1;
	}
      }
      /*********************************************************/
      /* Iterate through go1's interactions and count those    */
      /* that are also annotated to go2                        */    
      /*********************************************************/
      both=0;
      same=0;
      int i;
      for (i = 0; i < go1count; i++) {
	SV **SVinteractor = av_fetch(go1_inter_list,i,0);
	SV *svinter;
	if(SVinteractor!=NULL){
	  svinter=*SVinteractor;
	  char *interactor[40];
	  strcpy(interactor,SvPV(svinter,PL_na));
	  if(hashref4 != NULL &&
	     hv_exists_ent(go2_oo_inters,svinter,0)!=0){
	    both++;
	    /******************************************************/
	    /* Count interactions between the same protein	   */
	    /******************************************************/
	    char *inter1=strtok(interactor, "+");
	    char *inter2=strtok(NULL, "+");
	    if(strcmp(inter1,inter2)==0){
	      same++;
	    }
	  }
	}
      }
	
     /**********************************************************/
     /* Get the total number of interactions for this ontology */
     /**********************************************************/
      int tot_onto=0;
      SV **tot_onto_ref = hv_fetch(onto_count, oo, strlen(oo),0);
      if(tot_onto_ref != NULL){
	tot_onto=SvIV(*tot_onto_ref);
      }     

      /*******************************************************/
      /* Get the number of interactions that bridge this oo  */
      /* and that are annotated to go2			     */
      /*******************************************************/
      SV **go2_inter_list_ref=hv_fetch(inter_golist, go2, strlen(go2),0);
      AV *go2_inter_list;
      int go2count=0;
      if(go2_inter_list_ref!= NULL){
	HV *go2_inter_hash=(HV*)SvRV(*go2_inter_list_ref);
	SV **go2_oo_inter_list_ref=hv_fetch(go2_inter_hash,oo,strlen(oo),0);
	if(go2_oo_inter_list_ref!= NULL){
	  go2_inter_list=(AV*)SvRV(*go2_oo_inter_list_ref);
	  go2count=av_len(go2_inter_list)+1;
	}
      }
      
    inversed==0 ?
       	printf("%s\t%d\t%d\t%d\t%d\t%d\n",pair,tot_onto,go1count,go2count,both,same) :
       	printf("%s\t%d\t%d\t%d\t%d\t%d\n",pair,tot_onto,go2count,go1count,both,same);
    } //end for counter2
 
  } //end for counter1
  fflush(stdout);
  return;
 }
