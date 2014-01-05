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
# $interactions_gc{$go}{$oo}= number of inters annotated to go in this oo.     #
# $interactions_gl{$go}{$oo}= @inters that this bait/target is annotated to.  #
# $go_counts{ONTOs}{$o}= num of interactions annotated to $o                  #
# $go_counts{GOs}{$go}{$oo}= num of interactions between any prot annotatated #
#                            to $go (onto1) and any annotated to onto2        #
###############################################################################

use Getopt::Std;

use IO::File;
use strict;
my %opts;
getopts('epva:dG:g:m:n:o:s:',\%opts);
my $gaf_file=$opts{a}||die("Need a gene association file (-a)\n");
my $gen_file=$opts{g}||die("Need a genealogy file (-g)\n");
my $net_file=$opts{n}||die("Need a network file (-n)\n");
my $synfile=$opts{m}||die("Need a map file (-m)\n");
my $subonto=$opts{o}||'PFC';
my $verbose=$opts{v}||0;
my $debug=$opts{d}||0;
my $alt_go_file=$opts{G}||"./data/GO.terms_alt_ids";
my $exp_codes=$opts{e}||undef; ## Use only good evidence codes
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
my %prots=%{$r_prots};
my %gos=%{$r_gos};


###############################
# Load and parse Network file #
###############################
my ($r1,$r2,$r3,$r4)=load_network($net_file,1);
my %interactions_o=%{$r1};
my %interactions_g=%{$r2};
my %interactions_gl=%{$r3};
my %interactions_gc=%{$r4};
#######################################################
# At this point all information has been collected    #
# except the overlap between each go ($both). So, get #
# this and print				      #
#######################################################


# print "PP \$interactions_o{'PP'}{'GTR1_HUMAN+SUMO1_HUMAN'} : $interactions_o{'PP'}{'GTR1_HUMAN+SUMO1_HUMAN'}\n";
# print "PP \$interactions_g{'GO:0009725'}{'PP'}{'GTR1_HUMAN+SUMO1_HUMAN'} : $interactions_g{'GO:0009725'}{'PP'}{'GTR1_HUMAN+SUMO1_HUMAN'}\n";
# print "PP \$interactions_gc{'GO:0009725'}{'PP'} : $interactions_gc{'GO:0009725'}{'PP'}\n";


################################################################
# Get the number of interactions for each ontology combination #
################################################################


my @GOs=keys(%{$go_ontos{ONTO}});
my @ontos=split(//,$subonto);
my %foo;
foreach my $o1 (@ontos) {
    foreach my $o2 (@ontos) {
	my $oo = join("", sort {$b lt $a} ($o1,$o2));
	$foo{$oo}++;
    }
}
@ontos=keys(%foo);

foreach my $go (@GOs) {
    foreach my $oo (@ontos){
	my $num=0;
	$num=$interactions_gc{$go}{$oo} if defined $interactions_gc{$go}{$oo};
	print "$go\t$oo\t$num\n";
    }
}



#####################################################################
#                          SUBROUTINES                              #
#####################################################################
sub debug{
    print STDERR "@_\n" if $debug;

}
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
      ##$a[1] is the alternate/obsolete term names (if any)
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
  ## Quickly read through the network file to
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
  my (%interactors_o,%interactors_g, %interactors_gl,%interactors_gc,%h,%i);
  my @ontos=keys(%want);
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
      debug("BAIT: $bait, TARGET:$target");
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
		  #$interactors_o{$oo}{$hash_id}++;
		  
		  foreach my $go (@inter_gos){
		      debug("WAS:$go");
		      $go=$go_ontos{SYNONYM}{$go};
		      debug("IS:$go");
		      #$interactors_g{$go}{$oo}{$hash_id}++;
		      $interactors_gc{$go}{$oo}++;
		      ########################################################
		      # There is a bug in Inline.c that causes memory leaks  #
		      # when iterating through a hash. So, get the LIST of   #
		      # each interaction's GOs			             #
		      ########################################################
		      #push @{$interactors_gl{$go}{$oo}},$hash_id;
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
  $mode==0 ? return(\%h) : return(\%interactors_o, \%interactors_g, \%interactors_gl, \%interactors_gc);

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
