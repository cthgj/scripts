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
getopts('va:d:G:g:o:',\%opts);
my $gaf_file=$opts{a}||die("Need a gene association file (-a)\n");
my $gen_file=$opts{g}||die("Need a genealogy file (-g)\n");
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
@GOs=keys(%{$go_ontos{ONTO}});


print "#     GO\tTotal\tC\tF\tP\tCC\tCF\tCP\tFF\tFP\tPP\n" ;

c_loop($#GOs,$verbose,\@GOs,\%gos,\%{$go_ontos{ONTO}}, \%{$go_counts{ONTOs}}, \%{$go_counts{GOs}}, \%{$prots{GOs}}, \%wanted_proteins);
print STDERR "Finished correctly\n";
exit;
#####################################################
# This is the algorithm coded in perl for debugging #
#####################################################
for (my $i=0; $i<scalar(@GOs);$i++) {
  print STDERR "$i of $#GOs\r" if $verbose;
  for (my $k=$i+1; $k<scalar(@GOs);$k++) {
    my @a=sort {$b lt $a} ($GOs[$i],$GOs[$k]);
    my $pair = join("_", @a);
    my $oo = join("", sort {$b lt $a} ($go_ontos{ONTO}{$a[0]},$go_ontos{ONTO}{$a[1]}));
    ## Now go through the prots annotated to each go
    my $both=0;
    my %seen;
    foreach my $prot (@{$gos{$GOs[$i]}} ) {
      next if $seen{$prot};
      next unless $wanted_proteins{$oo}{$prot};
      $seen{$prot}++;
      if (defined($prots{GOs}{$prot}{$GOs[$k]})) {
	$both++;
      }
    }
    $go_counts{GOs}{$a[0]}{$oo}=0 unless defined($go_counts{GOs}{$a[0]}{$oo});
    $go_counts{GOs}{$a[1]}{$oo}=0 unless defined($go_counts{GOs}{$a[1]}{$oo});
    $go_counts{ONTOs}{$oo}="0" unless defined($go_counts{ONTOs}{$oo});
    print "$pair\t$go_counts{ONTOs}{$oo}\t$go_counts{GOs}{$a[0]}{$oo}\t$go_counts{GOs}{$a[1]}{$oo}\t$both\n";
  }
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
	  # annotations in this ontology	       #
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
	  # If they are different, we want	        #
	  # this protein if it has >=1 DIRECT        #
	  # annotation in each ontology	        #
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
  char go1[12],  oo[3];

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

    int Ccount=0, Fcount=0, Pcount=0, CCcount=0, CFcount=0, CPcount=0, FFcount=0, FPcount=0, PPcount=0, tot=0;


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
      tot=av_len(go1_prot_list);
    }
    /******************************************************************/
    /* Get the $go_counts{GOs}{$go1} hash, keys are onto combinations */
    /******************************************************************/
    HV *go1_count_hash;
   
    SV **tmp1=hv_fetch( go_counts, go1, strlen(go1), 0);
    if(tmp1!=NULL){
      go1_count_hash=(HV*)SvRV(*tmp1);
      char *oo1="CC";
      char *oo2="CF";
      char *oo3="CP";
      char *oo4="PP";
      char *oo5="FP";
      char *oo6="FF";
      /***********************************************************************/
      /* Get the total number of proteins annotated to this onto combination */
      /***********************************************************************/
      int tot_onto=0;
      if(hv_fetch( onto_count, oo1 , strlen(oo1), 0) != NULL){
	CCtot_onto=SvIV(*hv_fetch( onto_count, oo1 , strlen(oo1), 0));
      }
      if(hv_fetch( onto_count, oo2 , strlen(oo2), 0) != NULL){
	CFtot_onto=SvIV(*hv_fetch( onto_count, oo2 , strlen(oo2), 0));
      }
      if(hv_fetch( onto_count, oo3 , strlen(oo3), 0) != NULL){
	CPtot_onto=SvIV(*hv_fetch( onto_count, oo3 , strlen(oo3), 0));
      }
      if(hv_fetch( onto_count, oo4 , strlen(oo4), 0) != NULL){
	PPtot_onto=SvIV(*hv_fetch( onto_count, oo4 , strlen(oo4), 0));
      }
      if(hv_fetch( onto_count, oo5 , strlen(oo5), 0) != NULL){
	FPtot_onto=SvIV(*hv_fetch( onto_count, oo5 , strlen(oo5), 0));
      }
      if(hv_fetch( onto_count, oo6 , strlen(oo6), 0) != NULL){
	FFtot_onto=SvIV(*hv_fetch( onto_count, oo6 , strlen(oo6), 0));
      }
      
      int go1count=0;
      SV **aa=hv_fetch( go1_count_hash, oo1 , strlen(oo1), 0);
      if(aa != NULL)
	CCcount=SvIV(*aa);
      
      aa=hv_fetch( go1_count_hash, oo1 , strlen(oo2), 0);
      if(aa != NULL)
	CFcount=SvIV(*aa);
      
      aa=hv_fetch( go1_count_hash, oo3 , strlen(oo3), 0);
      if(aa != NULL)
	CPcount=SvIV(*aa);
      
      aa=hv_fetch( go1_count_hash, oo4 , strlen(oo4), 0);
      if(aa != NULL)
	PPcount=SvIV(*aa);
      
      aa=hv_fetch( go1_count_hash, oo5 , strlen(oo5), 0);
      if(aa != NULL)
	FPcount=SvIV(*aa);
      
      aa=hv_fetch( go1_count_hash, oo6 , strlen(oo6), 0);
      if(aa != NULL)
	FFcount=SvIV(*aa);
      
    }
    /*************************************************/
    /* If there are no proteins annotated to this go */
    /*************************************************/
    

printf("%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n",go1,tot,Ccount, Fcount, Pcount, CCcount, CFcount, CPcount, FFcount, FPcount, PPcount);





  } //end for counter1
  fflush(stdout);
}
