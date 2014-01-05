#!/usr/bin/perl 
## Get the numbers necessary to calculate the probabilities of association of each GO pair
#use strict;
use Getopt::Std;
use Inline C;
my %opts;
getopts('hAaOvCi:o:c:s:g:S:d:',\%opts);

my $interval=$opts{i}||500; ## interval for c_loop
my $gaf_file=$ARGV[0]||"gene_association.goa_human";
my $obo_file=$ARGV[1]||"biological_process.obo";
my $gen_file=$ARGV[2];
my $subonto=$opts{s}||'P';
my $over=$opts{O}||undef; ## print OVERrepresented pairs (control)
my $all=$opts{a}||undef;
my $verbose=$opts{v}||undef;
&usage() if $opts{h}; 
my $All_ontos=$opts{A}||undef; ## calculate cross ontology probabilities
my $calculated=$opts{c}||undef; ## useful when I already have some probabilities calculated 
## and I only want the rest. This option will pass a file
## in which the 1st characters up to the 1st space are the
## gopairs I have. 

my $bdir=$opts{d}||'.'; ## base-dir for output files
my $alt_go_file=$opts{g}||"./data/GO.terms_alt_ids";

my $specified_ontology=$opts{o}||undef; ## print only specific ontology combinations, eg FP
my %want;
if($specified_ontology){
    $specified_ontology=join("", sort {$b lt $a} split(//,$specified_ontology));
    map{$want{$_}++}(split(//,$specified_ontology));
}
else{$want{P}=$want{C}=$want{F}=1;}
## Many go pairs have exactly the same stats (eg, 12159_0_0_0_0), 
## to avoid doing the same calculation multiple times in R, 
## I condense them so that instead of having a line like:
##      go1_go2\t$tot\$go1\t$go2\t$both etc
## I have a line like:
##     12159_0_0_0_0\t12159\t0\t0\t0 etc 
## and create a file mapping each go to the appropriate number combination
my $condense=$opts{C}||undef; 

my $species=$opts{S}||undef;
die("Need a species name (-S) for a file suffix if condensing (-C)\n") unless $species;



## Must do this first, as it loads the GO synonyms as well
my %go_ontos=&get_go_ontologies($alt_go_file);
my @GOs=keys(%{$go_ontos{ONTO}});
my ($i,$k)=(0,0);

# timethese (10000, { 
#     'perl' => sub {
	

#     } , 
#     'C' => sub {	print "aa";},
# 	  } );
# die();

# while(my $pair=get_pair($i,$k,@GOs)){
#     print "P: $pair\n";
#     $i++;
#     $k++;
# }

my $aa;

my %papas=&load_genealogy($gen_file);
my %prots=&parse_gaf($gaf_file);
my %go_counts=&count_gos();


################## Main Program ##########################
 @GOs=keys(%{$go_ontos{ONTO}});

## Chose only those gos we want ($specified_ontology)
my @GOs=grep{defined($want{$go_ontos{ONTO}{$_}})}keys(%{$go_ontos{ONTO}});

print STDERR "GOs: " . scalar(@GOs) . "\n" if $verbose; 


#Now, calculate all numbers
my %have;
my (@out,@map);
my $tot=scalar(@GOs);

if($opts{i}){
    for (my $ii=7; $ii<$tot; $ii+=$interval){
	c_loop_lim($ii,$interval,$verbose,\%{$go_ontos{ONTO}},\%{$go_counts{GOs}},\%{$prots{GOs}},\%{$go_counts{ONTOs}},@GOs);
    }
}
else{
# start at 7 because 0==$ii, 1==$interval, 2==verbose, 3==hash_ref1 etc
    c_loop($verbose,\%{$go_ontos{ONTO}},\%{$go_counts{GOs}},\%{$prots{GOs}},\%{$go_counts{ONTOs}},@GOs);
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
	if ($specified_ontology){next unless defined($want{$onto})}
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
		my $oo = join("", sort {$b lt $a} ($os[$n],$os[$k]));
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
		    my $oo = join("", sort {$b lt $a} ($os[$n],$os[$k]));
		    #print "$prot \$go_counts{GOs}{$go}{$oo}++\n";
		    $go_counts{GOs}{$go}{$oo}++ if $this_prot{$oo}>0;
		    #print "$go ($go_ontos{ONTO}{$go}): $prot : $oo : $this_prot{$oo}\n";
		}
	    }
	}
	
    } ## end foreach $prot
    map{print STDERR "$_ : $go_counts{ONTOs}{$_}\n"}keys(%{$go_counts{ONTOs}}) if $verbose;
    
    return(%go_counts);
}


sub usage{
    print STDERR "Not written yet\n"; exit();


}



__END__
__C__
#include <string.h>
#include <stdlib.h>

int get_both( SV* hash_ref1, SV* hash_ref2){
  HV* hash1;
  HV* hash2;
  HE* hash_entry1;
  HE* hash_entry2;
  int num_keys1, i, ii, num_keys2;
  SV* sv_key1;
  SV* sv_key2;
  SV* sv_val1;
  SV* sv_val2;
  hash1 = (HV*)SvRV(hash_ref1);
  hash2 = (HV*)SvRV(hash_ref2);
  num_keys1 = hv_iterinit(hash1); // start hash iteration
  num_keys2 = hv_iterinit(hash2);
  //compare hash values
  int both=0;
  for (i = 0; i < num_keys1; i++) {

    hash_entry1 = hv_iternext(hash1);
    sv_key1 = hv_iterkeysv(hash_entry1);
    sv_val1 = hv_iterval(hash1, hash_entry1);
    ii=0;
    //reset hash2 iteration, otherwise we are already at the end of 
    //the hash (hv_iternext(hash2) is already at the end) because of 
    //the previous time we looped it 
    num_keys2 = hv_iterinit(hash2);
    
    for (ii = 0; ii < num_keys2; ii++) {
      hash_entry2 = hv_iternext(hash2);
      sv_key2 = hv_iterkeysv(hash_entry2);
      sv_val2 = hv_iterval(hash2, hash_entry2);
      if(strcmp(SvPV(sv_key1, PL_na),SvPV(sv_key2, PL_na))==0)
	both++;
    }
  }
    return(both);
}


 
void c_loop(int verbose,SV* hash_ref1, SV* hash_ref2, SV* prots_hash_ref, SV* onto_countref, SV* name1,...){
  HV* go_ontos, *go_counts;
  HE* onto_entry, *hash_entry2;
  char pair[100], go1[12], go2[12], onto[2], oo[3];
  int onto_keys, i, ii, k, num_keys2,go1count,go2count,both,tot_onto,b_go1count,b_go2count,b_both,b_tot_onto;
  float o; //overrepresentation
  SV *go1_onto, *go2_onto, *newhashref, *go1cnt, *go2cnt, **go1_protss , **go2_protss, *go1_prots, *go2_prots;
  HV *go1_count_hash, *go2_count_hash,*prots,*onto_count;
  Inline_Stack_Vars;

  go_ontos = (HV*)SvRV(hash_ref1);
  go_counts = (HV*)SvRV(hash_ref2);
  onto_keys = hv_iterinit(go_ontos); // start hash iteration
  num_keys2 = hv_iterinit(go_counts);
  prots=(HV*)SvRV(prots_hash_ref); // get the %prots{GOS} hash

/****************************************************/
/* This is the main loop that will iterate through  */
/* the list of GO terms and build ALL possible	    */
/* GO term pairs				    */
/****************************************************/
  for (i = 5; i <Inline_Stack_Items; i++) {
    if(verbose==1)
      fprintf(stderr,"%d noliim %d\r",i,Inline_Stack_Items);
    //get go1
    strcpy(go1,SvPV(Inline_Stack_Item(i),PL_na));  
    // get ontology for go1
    go1_onto=SvPV(*hv_fetch( go_ontos , go1 , strlen(go1), 0),PL_na);
    // Get the %go_count{$go1} hash. It will exist
    // if at least one protein is annotated to go1
    if(hv_exists_ent(go_counts,Inline_Stack_Item(i),0)>0){
      newhashref=*hv_fetch( go_counts, go1, strlen(go1), 0);
      go1_count_hash=(HV*)SvRV(newhashref);
    }
    // Get the number of proteins annotated to go1
    go1_protss=hv_fetch( prots, go1 , strlen(go1), 0);
      
    for (k = i+1; k <Inline_Stack_Items; k++) {
      //get go2
      strcpy(go2,SvPV(Inline_Stack_Item(k),PL_na));
      //Get GO pair
      strcpy(pair,go1);
      strcat(pair,"_");
      strcat(pair,go2);
      // get ontology for go2
      go2_onto=SvPV(*hv_fetch( go_ontos , go2 , strlen(go1), 0),PL_na);

            
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

  
      /***************************************/
      /* Get the number of proteins anotated */
      /* to each of the go terms	     */
      /***************************************/
      /****************************************************/
      /* The relevant hash entry in %go_count{$go1} will  */
      /* only exist if at least one prot (of those that   */
      /* have >1 DIRECT annotation from each of the two   */
      /* current ontologies (oo)) is annotated to go1 	  */
      /****************************************************/

      //hv_exists_ent needs an SV*, go1 will not work
      // but Inline_Stack_Item(i) does.
      go1count=0;
      if(hv_exists_ent(go_counts,Inline_Stack_Item(i),0)>0){
	if(hv_exists(go1_count_hash,oo,strlen(oo))>0){
	  go1cnt=*hv_fetch( go1_count_hash, oo , strlen(oo), 0);
	  go1count=SvIV(go1cnt);
	}
      }
      go2count=0;
      if(hv_exists_ent(go_counts,Inline_Stack_Item(k),0)>0){
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
      /* to an SV* and call get_both				     */
      /***************************************************************/
      both=0;
      if(go1_protss != NULL  &&
	 go2_protss != NULL){
	go1_prots=*go1_protss;
	go2_prots=*go2_protss;
	both=get_both(go1_prots,go2_prots);
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
      
      /****************/
      /* Print output */
      /****************/
      printf("%s\t%d\t%d\t%d\t%d\t%d\n",pair,tot_onto,go1count,go2count,both,o);

    } // end for k
    
  } // end for i
}
void c_loop_lim(int start,int interval,int verbose,SV* hash_ref1, SV* hash_ref2, SV* prots_hash_ref, SV* onto_countref, SV* name1,...){
  HV* go_ontos;
  HV* go_counts;
  HE* onto_entry;
  HE* hash_entry2;
  int onto_keys, i, ii, k, num_keys2;
  SV* sv_key1;
  SV* sv_key2;
  SV* sv_val1;
  SV* sv_val2;
  SV* gopair;
  char pair[100], go1[12], go2[12], onto[2], oo[3];
  int go1count,go2count,both,tot_onto,b_go1count,b_go2count,b_both,b_tot_onto;
  float o; //overrepresentation
  SV *go1_onto, *go2_onto, *newhashref, *go1cnt, *go2cnt, **go1_protss , **go2_protss, *go1_prots, *go2_prots;
  HV *go1_count_hash, *go2_count_hash,*prots,*onto_count;
  Inline_Stack_Vars;

  go_ontos = (HV*)SvRV(hash_ref1);
  go_counts = (HV*)SvRV(hash_ref2);
  onto_keys = hv_iterinit(go_ontos); // start hash iteration
  num_keys2 = hv_iterinit(go_counts);
  prots=(HV*)SvRV(prots_hash_ref); // get the %prots{GOS} hash

/****************************************************/
/* This is the main loop that will iterate through  */
/* the list of GO terms and build ALL possible	    */
/* GO term pairs				    */
/****************************************************/
  int lim=start+interval;
  if(lim>Inline_Stack_Items)
    lim = Inline_Stack_Items-1 ;
  for (i = start; i <=lim; i++) {
    if(verbose==1)
      fprintf(stderr,"%d lim %d\r",i,Inline_Stack_Items);
    //get go1
    strcpy(go1,SvPV(Inline_Stack_Item(i),PL_na));  
    // get ontology for go1
    go1_onto=SvPV(*hv_fetch( go_ontos , go1 , strlen(go1), 0),PL_na);
    // Get the %go_count{$go1} hash. It will exist
    // if at least one protein is annotated to go1
    if(hv_exists_ent(go_counts,Inline_Stack_Item(i),0)>0){
      newhashref=*hv_fetch( go_counts, go1, strlen(go1), 0);
      go1_count_hash=(HV*)SvRV(newhashref);
    }
    // Get the number of proteins annotated to go1
    go1_protss=hv_fetch( prots, go1 , strlen(go1), 0);
      
    for (k = i+1; k <Inline_Stack_Items; k++) {
      //get go2
      strcpy(go2,SvPV(Inline_Stack_Item(k),PL_na));
      //Get GO pair
      strcpy(pair,go1);
      strcat(pair,"_");
      strcat(pair,go2);
      // get ontology for go2
      go2_onto=SvPV(*hv_fetch( go_ontos , go2 , strlen(go1), 0),PL_na);

            
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

  
      /***************************************/
      /* Get the number of proteins anotated */
      /* to each of the go terms	     */
      /***************************************/
      /****************************************************/
      /* The relevant hash entry in %go_count{$go1} will  */
      /* only exist if at least one prot (of those that   */
      /* have >1 DIRECT annotation from each of the two   */
      /* current ontologies (oo)) is annotated to go1 	  */
      /****************************************************/

      //hv_exists_ent needs an SV*, go1 will not work
      // but Inline_Stack_Item(i) does.
      go1count=0;
      if(hv_exists_ent(go_counts,Inline_Stack_Item(i),0)>0){
	if(hv_exists(go1_count_hash,oo,strlen(oo))>0){
	  go1cnt=*hv_fetch( go1_count_hash, oo , strlen(oo), 0);
	  go1count=SvIV(go1cnt);
	}
      }
      go2count=0;
      if(hv_exists_ent(go_counts,Inline_Stack_Item(k),0)>0){
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
      /* to an SV* and call get_both				     */
      /***************************************************************/
      both=0;
      if(go1_protss != NULL  &&
	 go2_protss != NULL){
	go1_prots=*go1_protss;
	go2_prots=*go2_protss;
	both=get_both(go1_prots,go2_prots);
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
      
      /****************/
      /* Print output */
      /****************/
      printf("%s\t%d\t%d\t%d\t%d\t%d\n",pair,tot_onto,go1count,go2count,both,o);

    } // end for k
  } // end for i
  Inline_Stack_Void;
  return(0);
}
 
      /*******This prints HASH(XX) style references******/
      /* L1066 is :: HASH(0x557f558) : GO:0046034       */
      /* /\**************************************************\/          */
      /* HV* prots=(HV*)SvRV(prots_hash_ref); */
      /* SV* go1_prots=*hv_fetch( prots, go1 , strlen(go1), 0); */
      /* int anums = hv_iterinit(prots); */
      /* int l; */
      /* for (l = 0; l < anums; l++) { */
      /* 	HE* hashD = hv_iternext(prots); */
      /* 	HE* hashC = hv_iternext(prots); */
      /* 	//	int bb=get_both(hv_iterval(prots,hashD),hv_iterval(prots,hashC)); */
      /* 	printf("L%d is :: %s : %s\n",l,SvPV(hv_iterval(prots,hashD),PL_na),SvPV(hv_iterkeysv(hashD),PL_na )  ) ; */
      /* } */
      /*********************************************/

      /*******This initiliazes get_both succesfully******/
      /**************************************************/         
      /* HV* prots=(HV*)SvRV(prots_hash_ref); */
      /* SV* go1_prots=*hv_fetch( prots, go1 , strlen(go1), 0); */
      /* int anums = hv_iterinit(prots); */
      /* int l; */
      /* for (l = 0; l < anums; l++) { */
      /* 	HE* hashD = hv_iternext(prots); */
      /* 	HE* hashC = hv_iternext(prots); */
      /* 	int bb=get_both(hv_iterval(prots,hashD),hv_iterval(prots,hashC)); */
      /* } */
      /*********************************************/
  
      /* int anums = hv_iterinit(prots); */
      /* int l; */
      /* for (l = 0; l < anums; l++) { */
      /* 	HE* hashD = hv_iternext(prots); */
      /* 	HE* hashC = hv_iternext(prots); */
      /* 	int bb=get_both(hv_iterval(prots,hashD),hv_iterval(prots,hashC)); */
      /* } */
      




	/*******This works for hash iteration******/
	/* SV* newhashref=*hv_fetch( go_counts, go1, strlen(go1), 0); */
	/* HV* go1_count_hash=(HV*)SvRV(newhashref); */
	/* int nums = hv_iterinit(go1_count_hash); */
	/* int l; */
	/* for (l = 0; l < nums; l++) { */
	/*   HE* hashE = hv_iternext(go1_count_hash); */
	/*   printf("L%d/%d/%d is %s(%s):: %s : %s\n",i,k,l,go1,oo,SvPV(hv_iterval(go1_count_hash,hashE),PL_na),SvPV(hv_iterkeysv(hashE),PL_na )   ) ; */
	/* } */
	/*********************************************/


      //make a new SV from pair
      //gopair=newSVpvf(pair);


	//	SV* go1cnt=*hv_fetch( go1_count_hash, oo , strlen(oo), 0);
	//printf("pair: %s %s %s\n",SvPV(gopair,PL_na),oo,go1);	
	//int koko=0;
	//printf("KOKO \n");
	
      


