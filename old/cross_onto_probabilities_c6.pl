#!/usr/bin/perl -w
## Get the numbers necessary to calculate the probabilities of association of each GO pair
#use strict;
#use diagnostics;
use Benchmark ':all';
## Data Structures:
##   $go_ontos{SYNONYM}{$GO:term}= go term synonym
##   $go_ontos{ONTO}{$GO:term}= go term ontology
##   $papas{$GO:term} = @term_ancestor_list
##   $prots{$PROT}{GOs}{$GO:term} = will exist if $PROT is annotated to $GO:term
##   $prots{$PROT}{ONTOs}{$onto} = will exist if $PROT is annotated to $onto
##   @{$prot_annots{$prot}}= list of terms annotating protein
##   @{$gos{$GO:term}}= list of $PROTS annotated to $GO:term  
##   $go_counts{ONTOs}{$onto} = number of proteins with annots of that onto comb
##   $go_counts{GOs}{$GO:term}{$ontology}= number of proteins of that onto annot to term
##   $go_inter_counts{$GO:term}{$ontology}= number of interactions between a protein annotated to
##                                          $GO:TERM and one annotated to any term of other ontology
##   $interactors{$PROT1}{$PROT2}=will exist if $PROT1 interacts with $PROT2
##   $wanted_proteins{$onto}{$PROT}=will exist if this protein is wanted, i.e., if it has  
##                           >=1 DIRECT annotation in one of the desired ontologies and >=1 
##                           INTERACTOR annotated to the other OR if it has  >=1 DIRECT 
##                           annotation in each of the current ontologies
use Getopt::Std;

# use Inline (C => 'DATA',
# 	    CC => 'gcc',
# 	    CCFLAGS => '-O3',
#     );
use Inline C;
use IO::File;
my %opts;
getopts('haOvCiD:o:n:c:s:g:S:d:',\%opts);

print STDERR "\n====================\nREMEMBER::: COUNT SELF-SELF pairs for interactome\n====================\n\n";

my $gaf_file=$ARGV[0]||die("Need a gene association file (ARGV[0])\n");
my $gen_file=$ARGV[1]||die("Need a genealogy file (ARGV[1])\n");
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
my $debug=$opts{D}||0;

my $alt_go_file=$opts{g}||"./data/GO.terms_alt_ids";
my %wanted_proteins;
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


$synonyms{LOADED}=0;
my %go_ontos=&get_go_ontologies($alt_go_file);
my %papas=&load_genealogy($gen_file);
my ($r_prots,$r_gos)=&parse_gaf($gaf_file);
our %prots=%{$r_prots};
our %gos=%{$r_gos};
our (%interactors,%interactions,%go_inter_counts);
if($interactome){
    die("Need a synonyms file (-s) if working with a network\n") unless $synfile;
    my($r1,$r2, $r_counts)=&load_network($net_file); 
    %interactors=%{$r1};
    %interactions=%{$r2};
    %go_inter_counts=%{$r_counts};
}
my %go_counts;
unless($interactome){
    %go_counts=&count_gos();
}


my @ooo=keys(%want);
my $onto1=$ooo[0];
my $onto2=$ooo[1];
$verbose=0 unless $verbose;
if($interactome){
    my %temp_hash;
    foreach my $prot (keys(%prots)){
	$prot=&get_name($prot,"NET",2);
	next unless defined($interactors{$prot});
	foreach my $go (keys(%{$prots{$prot}{GOs}})){
	    $temp_hash{$go}++ if defined($want{$go_ontos{ONTO}{$go}});
	}
    }
    my @GOs=grep{defined($want{$go_ontos{ONTO}{$_}})}keys(%temp_hash);
    print STDERR "GOs: " . scalar(@GOs) . " ($#GOs)\n" if $verbose; 


    # my $cc=scalar(keys(%interactions));
    # my $pl=0;
    # foreach my $interaction (keys(%interactions)){
    # 	$pl++;
    # 	print STDERR "K:$pl of $cc\r";
    # 	my ($p1,$p2)=split(/-/,$interaction);
    # 	my ($o1,$o2)=("a","a");
    # 	my %have;
    # 	foreach my $go1 (keys(%{$prots{$p1}{GOs}})){
    # 	    #next if $go_ontos{ONTO}{$go1} eq $o1;
    # 	    #unless($go1 eq 'GO:0042824' || $go1 eq 'GO:0035672'){next}
    # 	    my $o1=$go_ontos{ONTO}{$go1};
    # 	    foreach my $go2 (keys(%{$prots{$p2}{GOs}})){
    # 		print "$go1 : $go2\n";
    # 		#unless($go2 eq 'GO:0042824' || $go2 eq 'GO:0035672'){next}
    # 		next if $go_ontos{ONTO}{$go2} eq $o2;
    # 		my $o2=$go_ontos{ONTO}{$go2};
    # 		my $oo = join("", sort {$b lt $a} ($o1,$o2));
    # 		$go_counts{GOs}{$go2}{$oo}++ unless $have{$go2};
    # 		$go_counts{GOs}{$go1}{$oo}++ unless $have{$go1};
    # 		$have{$go2}=$have{$go1}=1;
    # 	    }
    # 	}

    #}
    my $cc=scalar(keys(%interactions));
    my @inter_pairs=keys(%interactions);
    my @pps=keys(%interactors);
    my @lll=keys(%wanted_proteins);
    # map{print "ontos $_ : $go_ontos{ONTO}{$_}\n"}keys(%{$go_ontos{ONTO}});
    # map{print "inter $_ : $go_inter_counts{$_}\n"}keys(%go_inter_counts);
    # map{print "cONTO $_ : $go_counts{ONTOs}{$_}\n"}keys(%{$go_counts{ONTOs}});
    # map{print "gosss $_ : $gos{$_}\n"}keys(%gos);
    #map{print "wanted $_ : $wanted{$_}\n"}keys(%interactors);
    
    #print STDERR "$#GOs,$verbose,$subonto,\%{$go_ontos{ONTO}}," . \%go_inter_counts .",\%{$go_counts{ONTOs}}," . \%gos . "," . \%wanted_proteins . "," . \%interactors . "," . \@GOs . " go1: $GOs[0],$GOs[1]\n" if $debug;
    c_loop_inter2($#GOs,$verbose,$subonto,$debug,\%{$go_ontos{ONTO}},\%go_inter_counts,\%{$go_counts{ONTOs}},\%gos,\%wanted_proteins,\%interactors,\@GOs);
###########################################################################################
  exit;
    print "aaaaaaaaaaaaaaaaa malaka\n";
go1:for (my $n=0; $n<=$#GOs;$n++){
    my $go1=$go_ontos{SYNONYM}{$GOs[$n]};
    print STDERR "(p) $n of $#GOs\r" if $verbose; 

    unless($go1 eq 'GO:0042824' || $go1 eq 'GO:0035672'){print STDERR "skipping...\n" if $n==1; next}

    next if $go1 eq 'GO:0008150' ; ## skip P root
    next if $go1 eq 'GO:0005575' ; ## skip C root
    next if $go1 eq 'GO:0003674' ; ## skip F root 
  go2:for (my $k=$n+1; $k<=$#GOs;$k++){
      my $go2=$go_ontos{SYNONYM}{$GOs[$k]};
      unless($go2 eq 'GO:0042824' || $go2 eq 'GO:0035672'){next}

      next if $go2 eq 'GO:0008150' ; ## skip P root
      next if $go2 eq 'GO:0005575' ; ## skip C root
      next if $go2 eq 'GO:0003674' ; ## skip F root 
      my @b=sort {$b lt $a} ($go1,$go2);
      my $pair = join("_",@b);
      $go1=$b[0];
      $go2=$b[1];
      my $oo = join("", sort {$b lt $a} ($go_ontos{ONTO}{$go1},$go_ontos{ONTO}{$go2}));
      
      ## check that the output files do not exist
      if($n==0){
	  if(-e "$bdir/$oo.$species.cond" || -e "$bdir/$oo.$species.map"){
	      die("Output file $bdir/$oo.$species.cond or $bdir/$oo.$species.cond exists\n");
	  }
      }
      $go2='GO:0035672';
      $go1='GO:0042824';
      my $ga=scalar(@{$gos{$go1}})-1;
      my $gb=scalar(@{$gos{$go2}})-1;
      print STDERR "1P$go1:$ga\n";
      print STDERR "2P$go2:$gb\n";

    #my $both=get_both(\%{$prots{GOs}{$go1}},\%{$prots{GOs}{$go2}}); 
      my $both=0;
      my %sskip;
      open(P,">pp");
      foreach my $prot1 (@{$gos{$go1}}){
	  next unless defined($wanted_proteins{$oo}{$prot1});
	  foreach my $prot2 (@{$gos{$go2}}){
	      $sskip{$prot1}{$prot2}++;
	      $sskip{$prot2}{$prot1}++;
	      if($prot1 eq $prot2){next if $sskip{$prot1}{$prot2}>2;}
	      else{next if $sskip{$prot1}{$prot2}>1;}
	      next unless defined($wanted_proteins{$oo}{$prot2});
	      if(defined($interactors{$prot1}{$prot2}) ||
		 defined($interactors{$prot2}{$prot1})){
		  print "P:$prot1-$go1-$prot2-$go2\n";#" :: $wanted_proteins{$oo}{$prot1} :: $wanted_proteins{$oo}{$prot2} :: $oo\n";
		  $both++ ;
	      }
	  }
      }
      close(P);
      #next unless $both>0;
      my $a=$go_counts{GOs}{$go1}{$oo};
      my $b=$go_counts{GOs}{$go2}{$oo};
      print "p:$pair\t$go_counts{ONTOs}{$oo}\t$go_counts{GOs}{$go1}{$oo}\t$go_counts{GOs}{$go2}{$oo}\t$both\n";
  }
      
}

##########################################################################################
}
else{
   my @GOs=grep{defined($want{$go_ontos{ONTO}{$_}})}keys(%gos);
   print STDERR "GOs: " . scalar(@GOs) . "\n" if $verbose; 
   c_loop($#GOs,$verbose,$subonto,$debug,\%{$go_ontos{ONTO}},\%{$go_counts{GOs}},\%{$go_counts{ONTOs}},\%gos,\%wanted_proteins,\@GOs);
}

print STDERR "Exited correctly\n";

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
	    # next if $a[0] eq 'GO:0008150' ; ## skip P root
	    # next if $a[0] eq 'GO:0005575' ; ## skip C root
	    # next if $a[0] eq 'GO:0003674' ; ## skip F root 
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
    close(G);
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
    open(A,"$file")||die("Cannot open $file : $!\n");
    while(<A>){
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
    close(A);
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
    my (%prots,%gos);
    my @ontos=split(//,$subonto);
    while(<GAF>){
	next if /^!/;
#	unless(/GO:0042824/ ||/GO:0035672/){print STDERR "skipping...\n" if $.==1; next}

	print STDERR "Gaf: $.\r" if $verbose;
	chomp;
	my @tmparray=split(/\t/);
	## $tmparray[1]== name, $tmparray[4] == GO $tmparray[8] == ontology
	my $prot=&get_name($tmparray[1],"NET",0);
	# if($interactome){
	#     ##skip the prot unless it is in the network
	#     next unless defined($interactors{$prot});
	# }
	my $onto=$tmparray[8];
	my $go=$go_ontos{SYNONYM}{$tmparray[4]}||$tmparray[4];
	next unless defined($want{$onto});
	next unless $tmparray[3]=~/^$/;  ## skip the NOT annotations
	## This collects each prot's GOs
	#$prots{PROTS}{$prot}{GOs}{$go}++ unless $seen{$prot}{$go};
	$prots{$prot}{GOs}{$go}++ unless $seen{$prot}{$go};
	#print "\$prots{$prot}{GOs}{$go}:$prots{$prot}{GOs}{$go}\n";
	push @{$prot_annots{$prot}}, $go unless $seen{$prot}{$go};
	## This collects each prot's ontos
	$prots{$prot}{ONTOs}{$onto}++ unless $seen{$prot}{$go};
	
	## This collects each GO's prots (if it fulfills the criteria) 
	#$gos{$go}{$prot}++ unless $seen{$prot}{$go};
	push @{$gos{$go}},$prot unless $seen{$prot}{$go};
	#print "$go : $prot :  $seen{$prot}{$go} : @{$gos{$go}}\n";
	$seen{$prot}{$go}++;
	## Now add ancestor GOs
	if(exists($papas{$go})){
	    my @pp=@{$papas{$go}};
	    for(my $i=0; $i<=$#pp; $i++){
		$prots{$prot}{GOs}{$pp[$i]}++;
		#$gos{$pp[$i]}{$prot}++;
		push @{$gos{$pp[$i]}},$prot unless $seen{$prot}{$pp[$i]};
		push @{$prot_annots{$prot}}, $pp[$i] unless $seen{$prot}{$pp[$i]};
		$prots{$prot}{GOs}{$pp[$i]}++ unless $seen{$prot}{$go};
		$seen{$prot}{$pp[$i]}++;
	    }
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
    return(\%prots,\%gos);
}

############################################################
## count_gos                                              ##
############################################################
sub count_gos{
    my %go_counts;
    my $NONE=0;
    my @PROTS;
    if($interactome){@PROTS=keys(%interactors); }
    else{@PROTS=keys(%prots)}
    my $p=0;
    prot:foreach my $prot (@PROTS){
	$p++;
	print STDERR "Counting GOs:Prot $p of $#PROTS\r" if $verbose;
	$prot=&get_name($prot,"NET",2);
	unless ($prots{$prot}{GOs}){
	    $NONE++ ;
	   # print STDERR "$prot\txx\n";
	}
	#my %this_prot;
	#my @os= ('P','F','C');
	my @os=split(//,$subonto);
	for (my $n=0; $n<=$#os; $n++){
	    for (my $k=$n; $k<=$#os;$k++){
		my $oo = join("", sort {$b lt $a} ($os[$n],$os[$k]));
		if ($os[$n] eq $os[$k]){
		    if($interactome){
			## If $prot is annotated to the ontology
			## and has >=1 INTERACTOR annotated to the same
			if (exists($prots{$prot}{ONTOs}{$os[$n]})){
			    inter:foreach my $inter (keys(%{$interactors{$prot}})){
				if(defined($prots{$inter}{ONTOs}{$os[$n]})){
				  $go_counts{ONTOs}{$oo}++ ;
				  $wanted_proteins{$oo}{$prot}=1;
				  next prot;
				  last inter;
				}
			    }
			}
		    }
		    else{
			if(exists($prots{$prot}{ONTOs}{$os[$n]})  &&
			    $prots{$prot}{ONTOs}{$os[$n]}>1){
			    ## Count prots with >1 DIRECT annotations 
			    ## in the current ontologies
			    $go_counts{ONTOs}{$oo}++ ;
			    ## This protein has at least two
			    ## DIRECT annotations in the current onto
			    $wanted_proteins{$oo}{$prot}=1 ;
			}
		    }
		}
		else{
		    if($interactome){
			## If $prot is annotated to one of the ontologies
			## and has >=1 INTERACTOR annotated to the other
			if (exists($prots{$prot}{ONTOs}{$os[$n]})  &&
			    $prots{$prot}{ONTOs}{$os[$n]}>=1){
			  inter:foreach my $inter (keys(%{$interactors{$prot}})){
			      if(defined($prots{$inter}{ONTOs}{$os[$k]})){
				$go_counts{ONTOs}{$oo}++ ;
				$wanted_proteins{$oo}{$prot}=1;
				next prot;
				last inter;
			      }
			  }
			}
			elsif(exists($prots{$prot}{ONTOs}{$os[$k]})  &&
			      $prots{$prot}{ONTOs}{$os[$k]}>=1){
			  inter:foreach my $inter (keys(%{$interactors{$prot}})){
			      if(defined($prots{$inter}{ONTOs}{$os[$n]})){
				$go_counts{ONTOs}{$oo}++ ;
				$wanted_proteins{$oo}{$prot}=1;
				next prot;
				last inter;
			      }
			  }
			}	
		    }
		    else{
			## Count prots with >=1 DIRECT annotation 
			## in each of the current ontologies
			if(defined($prots{$prot}{ONTOs}{$os[$n]}) && 
			   defined($prots{$prot}{ONTOs}{$os[$k]})){
			    $go_counts{ONTOs}{$oo}++;
			    ## This protein has an annotation
			    ## in each of the current ontos
			    $wanted_proteins{$oo}{$prot}=1;
			}
		    }
		}
	    } ## end for my $k
	} ## end for my $n
	

	## Count occurences of each GO term
	## in EACH ontology/ontology combination
	## iff this prot has at least 2 DIRECT annotations
	## in that ontology/ontology combination

	unless($interactome){#&load_network() already does this
	    foreach my $go (keys(%{$prots{$prot}{GOs}})){
		## they should all be their synonym at this point, check just in case
		die("BAD synonym $go : $go_ontos{SYNONYM}{$go}\n") if $go ne $go_ontos{SYNONYM}{$go};
		for (my $n=0; $n<=$#os;$n++){
		    for (my $k=$n; $k<=$#os;$k++){
			my $oo = join("", sort {$b lt $a} ($os[$n],$os[$k]));
			
			## Count the number of INTERACTIONS involving this go and one from 
			## the other ontology. $wanted_proteins{$oo}{$prot}>0 iff
			## $prot is annotated to one of the ontologies and has 
			## >=1 INTERACTOR annotated to the other
			if($interactome){
			    if(exists($wanted_proteins{$oo}{$prot}) &&
			       $wanted_proteins{$oo}{$prot}>0){  
				foreach my $inter(keys(%{$interactors{$prot}})){
				    if(exists($wanted_proteins{$oo}{$inter}) &&
				       $wanted_proteins{$oo}{$inter}>0 &&
				       exists($prots{$inter}{ONTOs}{$go_ontos{ONTO}{$go}})){  
					$go_counts{GOs}{$go}{$oo}++;
					print "\$go_counts{GOs}{$go}{$oo} = $go_counts{GOs}{$go}{$oo}\n";
				    }
				}
			    }
			}
			
			else{
			    ## If annotations, $wanted_proteins{$oo}{$prot}>0 iff 
			    ## this prot has at least 2 DIRECT annotations 
			    ## in this $oo. Elsif interactome, iff $prot is 
			    ## annotated to one of the ontologies and has 
			    ## >=1 INTERACTOR annotated to the other
			    if(exists($wanted_proteins{$oo}{$prot}) &&
			       $wanted_proteins{$oo}{$prot}>0){
				$go_counts{GOs}{$go}{$oo}++;
			    }
			}
		    }
		}
	    }
	}	
    } ## end foreach $prot
    print STDERR "\n"if $verbose;
    map{print STDERR "$_ : $go_counts{ONTOs}{$_}\n"}keys(%{$go_counts{ONTOs}}) if $verbose;
    print STDERR "NONE : $NONE\n";
    return(%go_counts);
}

############################################################
## Load Network                                           ##
############################################################
sub load_network {
    my (%inters,%proteins, %counts);
    my $network_file=shift;
    open(A, "$net_file")|| die("Cannot open $net_file : $!\n");
    while(<A>){
	next if /^\d+$/;
	next if /^!/;
	next if /^\s*$/;
	# print STDERR "Loading Network $.\r" if $verbose;
	chomp;
	my ($bait,$target)=split(/\t/,$_);
	$bait=&get_name($bait,"NET",2);
	$target=&get_name($target,"NET",3);
	## %proteins will hold all the info on all the 
	## network's proteins
	#my $key=join("-",sort {$b lt $a} ($bait,$target));
	my @aa=sort{$b lt $a} $target,$bait;
	$proteins{$aa[0]}{$aa[1]}=1;
	#$proteins{$bait}{$target}=$proteins{$target}{$bait}=1;
	my $ll=$bait . "-" . $target;
	$inters{$ll}++;
    }
    close(A);
    my @ontos=keys(%want);
    my $ll=scalar(keys(%proteins));
    my $c=0;
    my $poss_comb=(scalar(keys(%want))* (scalar(keys(%want))-1)/2) ;
    my %local_seen;
    ## For each interacting protein
    my %both;
    foreach my $p1 (keys(%proteins)){
	$c++;
	print STDERR "$c of $ll (net)\r" if $verbose;
	my %seen_ontos;
	my %local_want;
	
	## Foreach bait's gos, count the num of interactions involving
	## it and any member of each ontology
	my @bait_gos=keys(%{$prots{$p1}{GOs}});
	my @targets=keys(%{$proteins{$p1}});
	#map{$local_seen{$p1}{$_}++;}@targets;
	my @w=keys(%want);
	## Count how many interactions this prot has with each of the 
	## other ontologies
	foreach my $target (@targets){
	    if($debug){print STDOUT "$p1 $target\n";}
	    my %seen_gos;
	    my @target_gos=keys(%{$prots{$target}{GOs}});
	    ##Foreach ontology combination
	    for(my $i=0; $i<scalar(@w);$i++){
		for(my $k=$i; $k<scalar(@w);$k++){
		    if(defined($prots{$p1}{ONTOs}{$w[$i]}) && 
		       defined($prots{$target}{ONTOs}{$w[$k]})){
			my $oo = join("", sort {$b lt $a} ($w[$i],$w[$k]));
			$go_counts{ONTOs}{$oo}++;
			$wanted_proteins{$oo}{$p1}=1;
			$wanted_proteins{$oo}{$target}=1;
			map{
			    $seen_gos{$_}{$oo}++;
			    $counts{$_}{$oo}++; 
			    if($debug){print STDOUT "a\$counts{$_}{$oo}:$counts{$_}{$oo}\n";}
			}@target_gos;
			$local_want{$oo}++;
		    }
		}
	    }
	    unless ($p1 eq $target){
		foreach my $oo (keys(%local_want)){
		    foreach my $bgo (@bait_gos){
			$counts{$bgo}{$oo}++ unless $seen_gos{$bgo}{$oo};
			if($debug){print STDOUT "b\$counts{$bgo}{$oo}:$counts{$bgo}{$oo}\n";}
		    }
		}
	    }
	    # if($p1 eq $target){
	    # 	map{
	    # 	    my $oo=$go_ontos{ONTO}{$_} . $go_ontos{ONTO}{$_};
	    # 	    $counts{$_}{$oo}++; 
	    # 	    print STDOUT "\$counts{$_}{$oo}:$counts{$_}{$oo}\n";
	    # 	}@bait_gos;
	    # }
	    # else{
	    # 	for(my $i=0; $i<scalar(@bait_gos);$i++){
		    


	    # 	    for(my $k=0; $k<scalar(@target_gos);$k++){
	    # 		my $oo = join("", sort {$b lt $a} ($go_ontos{ONTO}{$bait_gos[$i]},$go_ontos{ONTO}{$target_gos[$k]}));
	    # 		$counts{$bait_gos[$i]}{$oo}++;
	    # 		$counts{$target_gos[$k]}{$oo}++;
	    # 		print STDOUT "\$counts{$bait_gos[$i]}{$oo}:$counts{$bait_gos[$i]}{$oo}\n";
	    # 		print STDOUT "\$counts{$target_gos[$k]}{$oo}:$counts{$target_gos[$k]}{$oo}\n";
	    # 	    }
	    # 	}
	    # }
	    
	}
    }
    print STDERR "\n" if $verbose;
    return(\%proteins,\%inters, \%counts)
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
#include "uthash.h"

//c_loop($#GOs,$verbose,$subonto,$debug,\%{$go_ontos{ONTO}},\%{$go_counts{GOs}},\%{$go_counts{ONTOs}},\%gos,\%wanted_proteins,\@GOs);
/* SV* c_loop(int max,int verbose, SV* SVontology, int debug, SV* go_ontos_Oref, SV* go_counts_Gref, SV* go_counts_Oref, SV* gos_hash_ref, SV* wanted_prots_ref, SV* GOsSV){ */
/*   HE* onto_entry, *hash_entry2; */
/*   char pair[100], go1[12], go2[12], onto[2], oo[3], wanted_onto[4], id[22]; */
/*   int onto_keys, counter1,counter2, l, go1count,go2count,both,tot_onto,b_go1count,b_go2count,b_both,b_tot_onto,want; */
/*   float o; //overrepresentation */
/*   SV *go1_onto, *go2_onto, *newhashref, *go1cnt, *go2cnt, **go1_protss , **go2_protss, *go1_prots, *go2_prots, *GOsSV_ref1, **GOsSV_ref2, *GO1, **GOsSV2_ref2, *GO2; */
/*   HV *go_ontos, *go_counts, *go1_prots_hash, *go2_prots_hash, *go1_count_hash, *go2_count_hash,*gos,*onto_count, *wanted_prots; */
/*   go_ontos = (HV*)SvRV(go_ontos_Oref); */
/*   go_counts = (HV*)SvRV(go_counts_Gref); */
/*   gos=(HV*)SvRV(gos_hash_ref); // get the %gos hash */
/*   wanted_prots=(HV*)SvRV(wanted_prots_ref); // get the %wanted_proteins hash */
  
/*   //printf("go1_go2\t\ttot\t#go1\t#go2\t#both\to\n"); */
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
 
/*   for(counter1=0;counter1<max; counter1++){ */
/*     if(verbose==1) */
/*       fprintf(stderr,"%d of %d\r",counter1,max); */
    
/*     //get go1 */
/*     GOsSV_ref2=av_fetch(GOsSV_ref1,counter1,0); */
/*     GO1=*GOsSV_ref2 ; */
/*     strcpy(go1,SvPV(GO1,PL_na)); */
/*     // get ontology for go1 */
/*     go1_onto=SvPV(*hv_fetch( go_ontos , go1 , strlen(go1), 0),PL_na); */
    
/*     /\*******************************************************************\/ */
/*     /\* Check if we want this ontology				       *\/ */
/*     /\* 	                                                               *\/ */
/*     /\* size_t strspn( char *s1, const char *s2) :: returns the length  *\/ */
/*     /\* of the longest substring of s1 that begins at the start of s1   *\/ */
/*     /\* and consists only of  the characters found in s2.               *\/ */
/*     /\*******************************************************************\/ */
/*     want=strspn(go1_onto, wanted_onto); */
/*     if(want==0){ */
/* 	fprintf(stderr, "\nSkipping %s %s\n",go1_onto,go1); */
/* 	continue; */
/*     } */
/*     /\********************************************************\/ */
/*     /\* Get the %go_counts{GOs}{$go1} hash. It will exist    *\/ */
/*     /\* if at least one protein is annotated to go1          *\/ */
/*     /\********************************************************\/ */
/*     if(hv_exists_ent(go_counts,GO1,0)>0){ */
/*       //newhashref=$go_counts{GOs}{$go1} */
/*       newhashref=*hv_fetch( go_counts, go1, strlen(go1), 0); */
/*       go1_count_hash=(HV*)SvRV(newhashref);     */
/*     } */
/*     /\* Get the list of proteins annotated to go1 *\/ */
/*     SV** go1_protlist_rref=hv_fetch( gos, go1 , strlen(go1), 0);  */
/*     AV *go1_protlist; */
/*     if(go1_protlist_rref != NULL){ */
/*       go1_protlist=(AV*)SvRV(*go1_protlist_rref); */
/*      } */
/*     /\*************************************\/ */
/*     /\* Start iterating GOlist to get go2 *\/ */
/*     /\*************************************\/ */
/*     for(counter2 = counter1+1; counter2 <=max; counter2++) { */
/*      //get go2 */
/*       GOsSV2_ref2=av_fetch(GOsSV_ref1,counter2,0); */
/*       GO2=*GOsSV2_ref2 ; */
/*       strcpy(go2,SvPV(GO2,PL_na)); */
/*       //Get GO pair */
/*       strcpy(pair,go1); */
/*       strcat(pair,"_"); */
/*       strcat(pair,go2); */
    
/*       // get ontology for go2 */
/*       go2_onto=SvPV(*hv_fetch( go_ontos, go2 , strlen(go2), 0),PL_na); */
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
/*       if(want==0){ */
/* 	  fprintf(stderr, "\nSkipping %s %s\n",go2,oo); */
/* 	  continue; */
/*       } */
/*       /\*************************************************************\/ */
/*       /\* Are any proteins, of those that have >1 DIRECT annotation *\/ */
/*       /\* in this ontology, annotated to go1?                       *\/ */
/*       /\*************************************************************\/ */
/*       go1count=0;  // initialize go1 count to 0 */
/*       if(newhashref != NULL && */
/* 	 hv_exists(go1_count_hash,oo,strlen(oo))>0){ */
/* 	go1cnt=*hv_fetch( go1_count_hash, oo , strlen(oo), 0); */
/* 	go1count=SvIV(go1cnt); */
/*       } */
/*       /\********************************************************\/ */
/*       /\* Get the %go_counts{GOs}{$go2} hash. It will exist    *\/ */
/*       /\* if at least one protein is annotated to go2          *\/ */
/*       /\********************************************************\/ */
/*       go2count=0;  // initialize go1 count to 0 */
/*       if(hv_exists_ent(go_counts,GO2,0)>0){ */
/* 	//newhashref=$go_counts{GOs}{$go1} */
/* 	newhashref=*hv_fetch( go_counts, go2, strlen(go2), 0); */
/* 	go2_count_hash=(HV*)SvRV(newhashref); */
/* 	/\*************************************************************\/ */
/* 	/\* Are any proteins, of those that have >1 DIRECT annotation *\/ */
/* 	/\* in this ontology, annotated to go2?                       *\/ */
/* 	/\*************************************************************\/ */
/* 	if(newhashref != NULL && */
/* 	   hv_exists(go2_count_hash,oo,strlen(oo))>0){ */
/* 	  go2cnt=*hv_fetch( go2_count_hash, oo , strlen(oo), 0); */
/* 	  go2count=SvIV(go2cnt); */
/* 	} */
/*       } */
/*       /\****************************************************************\/ */
/*       /\* Now, if at least one protein with >1 DIRECT annotation in    *\/ */
/*       /\* each of the ontologies is annotated to each of the gos, then *\/ */
/*       /\* count the prots annotated to both gos                        *\/ */
/*       /\****************************************************************\/ */
/*       both=0; */
/*       if(go1count>0 && go2count>0){ */
/* 	/\* Get the list of proteins annotated to go2 *\/ */
/* 	SV** go2_protlist_rref=hv_fetch( gos, go2 , strlen(go2), 0);  */
/* 	AV *go2_protlist; */
/* 	if(go2_protlist_rref != NULL){ */
/* 	  go2_protlist=(AV*)SvRV(*go2_protlist_rref); */
/* 	} */
/* 	/\*************************************************************\/ */
/* 	/\* Iterate through the list of prots annotated to go1        *\/ */
/* 	/\*************************************************************\/ */
/* 	int cc; */
/* 	for(cc=0;cc<=av_len(go1_protlist);cc++){ */
/* 	  SV **go1_prot_rref=av_fetch(go1_protlist,cc,0); */
/* 	  SV *go1_prot_ref=*go1_prot_rref; */

/* 	  /\**************************************************************\/ */
/*           /\* Skip any proteins that do not have >1 DIRECT annotations   *\/ */
/*           /\* in this ontology combination				*\/ */
/*           /\**************************************************************\/ */
/* 	  SV *onto_wanted_prots_ref=*hv_fetch(wanted_prots,oo, strlen(oo), 0); */
/* 	  HV *onto_wanted_prots=(HV*)SvRV(onto_wanted_prots_ref); */
/* 	  if(onto_wanted_prots_ref != NULL && */
/* 	     hv_exists_ent(onto_wanted_prots,go1_prot_ref,0)==0){ */
/* 	    continue; */
/* 	  } */
/* 	  /\*************************************************************\/ */
/*           /\* Iterate through the list of prots annotated to go2	       *\/ */
/*           /\*************************************************************\/ */
/* 	  int cc2; */
/* 	  for(cc2=0;cc2<=av_len(go2_protlist);cc2++){ */
/* 	    SV **go2_prot_rref=av_fetch(go2_protlist,cc2,0); */
/* 	    SV *go2_prot_ref=*go2_prot_rref; */
/* 	    /\**************************************************************\/ */
/* 	    /\* Skip any proteins that do not have >1 DIRECT annotations   *\/ */
/* 	    /\* in this ontology combination			   	  *\/ */
/* 	    /\**************************************************************\/ */
/* 	    if(onto_wanted_prots_ref != NULL && */
/* 	       hv_exists_ent(onto_wanted_prots,go2_prot_ref,0)==0){ */
/* 	      continue; */
/* 	    } */
/* 	    if(go2_protlist_rref != NULL){ */
/* 	      go2_protlist=(AV*)SvRV(*go2_protlist_rref); */
	      
/* 	      /\****************************************************\/ */
/*               /\* both++ if a protein is annotated to both terms	  *\/ */
/*               /\****************************************************\/ */
/* 	      if(strcmp(SvPV(go1_prot_ref,PL_na),SvPV(go2_prot_ref,PL_na))==0) */
/* 		both++; */
/* 	    } */
/* 	  } */
/* 	} */
/*       } */
/*       /\******************************************************\/ */
/*       /\* Get the total prots for this ontology combination  *\/ */
/*       /\******************************************************\/ */
/*       onto_count=(HV*)SvRV(go_counts_Oref); */
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
/*       /\****************\/ */
/*       /\* Print output *\/ */
/*       /\****************\/ */
/*       printf("%s\t%d\t%d\t%d\t%d\t%d\n",pair,tot_onto,go1count,go2count,both,o); */
/*       if( both > go2count || */
/* 	  both > go1count){ */
/* 	printf("%d and %d\n",counter1,counter2); */
/* 	exit(0); */
/*       } */
/*     } // end for counter2 */
/*   } // end for counter1 */
/*   fflush(stdout); */
/*  } */

//c_loop_inter2(   $#GOs, $verbose,       $subonto,    $debug      \%{$go_ontos{ONTO}},\%go_inter_counts,\%{$go_counts{ONTOs}},\%gos,        \%wanted_proteins,  \%interactors,         \@GOs);
SV* c_loop_inter2(int max,int verbose, SV* SVontology, int debug, SV* go_ontos_Oref, SV* go_counts_Gref, SV* go_counts_Oref, SV* gos_hash_ref, SV* wanted_prots_ref, SV* inters_hash_ref, SV* GOsSV){
  HE* onto_entry, *hash_entry2;
  char pair[100], go1[12], go2[12], onto[2], oo[3], wanted_onto[4], id[22];
  int onto_keys, counter1,counter2, l, go1count,go2count,both,tot_onto=-1,b_go1count,b_go2count,b_both,b_tot_onto,want;
  float o; //overrepresentation
  SV *go1_onto, *go2_onto, *newhashref, *go1cnt, *go2cnt, **go1_protss , **go2_protss, *go1_prots, *go2_prots, *GOsSV_ref1, **GOsSV_ref2, *GO1, **GOsSV2_ref2, *GO2;
  HV *go_ontos, *go_counts, *go1_prots_hash, *go2_prots_hash, *go1_count_hash, *go2_count_hash,*inters,*gos,*onto_count, *wanted_prots;
  go_ontos = (HV*)SvRV(go_ontos_Oref);
  go_counts = (HV*)SvRV(go_counts_Gref); // %go_inter_counts
  gos=(HV*)SvRV(gos_hash_ref); // get the %gos hash
  inters=(HV*)SvRV(inters_hash_ref); // get the %interactors hash
  wanted_prots=(HV*)SvRV(wanted_prots_ref); // get the %wanted_proteins hash
  int bob;
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
  for(counter1=0;counter1<max; counter1++){ 
    if(verbose==1)
      fprintf(stderr,"(c) %d of %d\r",counter1,max);
    //get go1
    GOsSV_ref2=av_fetch(GOsSV_ref1,counter1,0);
    GO1=*GOsSV_ref2 ;
    strcpy(go1,SvPV(GO1,PL_na));
    if (debug>0)
      printf("1GGGGGGGGGGGGGGggggggg :%s\n",go1);
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
    /********************************************************/
    /* Get the %go_counts{GOs}{$go1} hash. It will exist    */
    /* if at least one WANTED protein is annotated to go1   */
    /********************************************************/
    if(hv_exists_ent(go_counts,GO1,0)>0){
      //newhashref=$go_counts{GOs}{$go1}
      newhashref=*hv_fetch( go_counts, go1, strlen(go1), 0);
      go1_count_hash=(HV*)SvRV(newhashref);
    }
    /* Get the list of proteins annotated to go1 */
    SV** go1_protlist_rref=hv_fetch( gos, go1 , strlen(go1), 0); 
    AV *go1_protlist;
    if(go1_protlist_rref != NULL){
      go1_protlist=(AV*)SvRV(*go1_protlist_rref);
     }
   
    /*************************************/
    /* Start iterating GOlist to get go2 */
    /*************************************/
    if (debug>0)
      printf("2GGGGGGGGGGGGGGggggggg :%s: (%d)\n",go1,counter1);

    for(counter2 = counter1+1; counter2 <=max; counter2++) {        
      int inversed=0;
      //get go2
      GOsSV2_ref2=av_fetch(GOsSV_ref1,counter2,0);
      GO2=*GOsSV2_ref2 ;
      strcpy(go2,SvPV(GO2,PL_na));
      if (debug>0)
	printf("3GGGGGGGGGGGGGGggggggg :%s:%s\n",go1,go2);
      //Get GO pair
      if(strcmp(go1,go2)<0){
	strcpy(pair,go1);
	strcat(pair,"_");
	strcat(pair,go2);
      }
      else{
	strcpy(pair,go2);
	strcat(pair,"_");
	strcat(pair,go1);
	inversed=1;
      }
      if(debug>0)
	printf("4GGGGGGGGGGGGGGggggggg :%s:%s (%d/%d)::%d\n",go1,go2,counter1,counter2, tot_onto);
      // get ontology for go2
      go2_onto=SvPV(*hv_fetch( go_ontos, go2 , strlen(go2), 0),PL_na);
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
      if(want==0)
	continue;
      /**************************************************************/
      /* Get the list of wanted prots for this ontology combination */
      /**************************************************************/
      SV **onto_wanted_prots_rref=hv_fetch(wanted_prots,oo, strlen(oo), 0);
      SV *onto_wanted_prots_ref;
      HV *onto_wanted_prots;
      if(onto_wanted_prots_rref != NULL){
	  onto_wanted_prots_ref=*onto_wanted_prots_rref;
	  onto_wanted_prots=(HV*)SvRV(onto_wanted_prots_ref);
      }
      else{continue;}
      /**************************************************************/
      /* How many proteins, of those that have >1 DIRECT annotation */
      /* in this ontology, are annotated to go1?                    */
      /**************************************************************/
      go1count=0;  // initialize go1 count to 0
      // if (exists($go_inter_counts{$go1}))
      if(hv_exists_ent(go_counts,GO1,0)>0 &&
	 hv_exists(go1_count_hash,oo,strlen(oo))>0){
	go1cnt=*hv_fetch( go1_count_hash, oo , strlen(oo), 0);
	go1count=SvIV(go1cnt);
      }
      //printf("go1: %s (%d)\n",go1,go1count);
      /* Get the list of proteins annotated to go2 */
      SV** go2_protlist_rref=hv_fetch( gos, go2 , strlen(go2), 0); 
      AV *go2_protlist;
      if(go2_protlist_rref != NULL){
	go2_protlist=(AV*)SvRV(*go2_protlist_rref);
      }
      if(debug>0)
	printf("C:%d/%d::%d\n",counter1,counter2,av_len(go1_protlist));        

      /********************************************************/
      /* Get the %go_inter_counts{$go2} hash. It will exist    */
      /* if at least one protein is annotated to go2          */
      /********************************************************/
      go2count=0;  // initialize go1 count to 0
      if(hv_exists_ent(go_counts,GO2,0)>0){      
	//newhashref=$go_counts{GOs}{$go1}
	newhashref=*hv_fetch( go_counts, go2, strlen(go2), 0);
	go2_count_hash=(HV*)SvRV(newhashref);
	/*************************************************************/
	/* Are any proteins, of those that have >1  annotation       */
	/* in this ontology, annotated to go2?                       */
	/*************************************************************/
	if(hv_exists(go2_count_hash,oo,strlen(oo))>0){
	  go2cnt=*hv_fetch( go2_count_hash, oo , strlen(oo), 0);
	  go2count=SvIV(go2cnt);
	}
      }      
      /****************************************************************/
      /* Now, count the number of interactions between these 2 gos    */
      /****************************************************************/
      both=0;   
      if(go1count>0 && go2count>0){ 

	/*************************************************************/
	/* Iterate through the list of prots annotated to go1        */
	/*************************************************************/
	int cc;
	for(cc=0;cc<=av_len(go1_protlist);cc++){
	  int seen=0;
	  SV **go1_prot_rref=av_fetch(go1_protlist,cc,0);
	  SV *go1_prot_ref;
	  if(go1_prot_rref != NULL){
	    go1_prot_ref=*go1_prot_rref;
	  }	  
	  /**************************************************************/
          /* Skip any proteins that do not have >=1 INTERACTORS         */
          /* in this ontology combination				*/
          /**************************************************************/
	  //Get prot 1
	  char prot1[14];
	  strcpy(prot1,SvPV(go1_prot_ref,PL_na));
	  /* if(hv_exists(onto_wanted_prots,a,strlen(a))==0){printf("aAasaasdasd\n");} */
	  /* else{prinf("BBBBbbbbbb\n");} */
	  if(onto_wanted_prots_ref != NULL &&
	     hv_exists_ent(onto_wanted_prots,go1_prot_ref,0)==0){
	    continue;
	  }
	  /*************************************************************/
	  /* Iterate through the list of prots annotated to go2        */
	  /*************************************************************/
	  int cc2;
	  for(cc2=0;cc2<=av_len(go2_protlist);cc2++){
	    SV **go2_prot_rref=av_fetch(go2_protlist,cc2,0);
	    SV *go2_prot_ref=*go2_prot_rref;
	    //Get prot 2
	    char prot2[30];
	    strcpy(prot2,SvPV(go2_prot_ref,PL_na));
	    /******************************************************/
	    /* Skip any proteins that do not have >=1 INTERACTORS */
	    /* in this ontology combination			  */
	    /******************************************************/
	    if(onto_wanted_prots_ref != NULL &&
	       hv_exists_ent(onto_wanted_prots,go2_prot_ref,0)==0){
	      continue;
	    }
	    /**********************************************************************************/
            /* When a protein interacts with itself, both $interactors{$prot1}{$prot2}	      */
	    /* and $interactors{$prot2}{$prot1} will exist and the interaction will be	      */
            /* counted twice. Avoid this.						      */
            /**********************************************************************************/
	    if(strcmp(prot1,prot2)==0){
	      if(seen>0){
		continue;
	      }
	      seen++;
	    }

	    /******************************************/
	    /* Get hash of prot1's interactors	      */
	    /******************************************/
	    SV** prot1_interhash_rref=hv_fetch( inters, prot1 , strlen(prot1), 0);
	    HV* prot1_interhash;
	    if( prot1_interhash_rref != NULL){
	      SV* ttemp=*prot1_interhash_rref;
	      prot1_interhash=(HV*)SvRV(ttemp);
	      /**********************************************/
	      /* Check if prot1 interacts with prot2        */
	      /**********************************************/
	      if(hv_exists(prot1_interhash,prot2,strlen(prot2))>0){
		/* inversed==0 ? */
		/*   printf("cc %s-%s-%s-%s\n",prot1,go1,prot2,go2) : */
		/*   printf("cc %s-%s-%s-%s\n",prot2,go2,prot1,go1); */
		both++;
	      }
	    }
	    if(debug>0){
	      printf("AA go1:%s go2:%s p1:%s p2:%s seen:%d, cc:%d cc2:%d both:%d\n",go1,go2,prot1,prot2,seen,cc,cc2,both); 
	    }
	    /******************************************/
	    /* Get hash of prot2's interactors	    */
	    /******************************************/
	    /* SV** prot2_interhash_rref=hv_fetch( inters, prot2 , strlen(prot2), 0); */
	    /* HV* prot2_interhash; */
	    /* if( prot2_interhash_rref != NULL){ */
	    /*   prot2_interhash=(HV*)SvRV(*prot2_interhash_rref); */
	    /* } */
	  }
	}
      }
      /******************************************************/
      /* Get the total prots for this ontology combination */
      /******************************************************/
      onto_count=(HV*)SvRV(go_counts_Oref);
      tot_onto=0;
      if(hv_fetch( onto_count, oo , strlen(oo), 0) != NULL){
	tot_onto=SvIV(*hv_fetch( onto_count, oo , strlen(oo), 0));
	if (debug>0)
	  printf("TOTONTO:%d, oo:%s\n ",tot_onto,oo);
      }
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
      inversed==0 ?
	printf("%s\t%d\t%d\t%d\t%d\n",pair,tot_onto,go1count,go2count,both) :
	printf("%s\t%d\t%d\t%d\t%d\n",pair,tot_onto,go2count,go1count,both);

      /* if(tot_onto<go1count||tot_onto<go2count){ */
      /* 	printf("pair:%s, onto:%s\n",pair,oo); */
      /* 	exit(0); */
      /* } */
	
      } // end for counter2
    } // end for counter1
    fflush(stdout);
 }

//&count_go_interactions1(    $p1      $go,    \%{$go_ontos{ONTO}}, $onto,          \@targets,   \%prot_annots,     \%local_seen );
int count_go_interactions(SV *SVprot,SV* SVgo, int debug, SV *go_ontos_Oref, SV *wanted_onto, SV *SVinters, SV *SVprot_annots, SV *SVseen){
  HV *go_ontos = (HV*)SvRV(go_ontos_Oref);
  HV *prot_annots = (HV*)SvRV(SVprot_annots); //%prot_annots
  HV *seen_prots = (HV*)SvRV(SVseen); //%prot_annots
  AV *inters=(AV*)SvRV(SVinters); 
  int counter1, counter2, counter3;
  char go[12], oo[3], prot[30];
  int interactions=0;
  //get go
  strcpy(go,SvPV(SVgo,PL_na));
  //get prot
  strcpy(prot,SvPV(SVprot,PL_na));
  // get go ontology
  char *go_onto=SvPV(*hv_fetch( go_ontos , go , strlen(go), 0),PL_na);
 printf("---\n");
  /********************************************/
  /* Iterate through this bait's target prots */
  /********************************************/
  for(counter1=0;counter1<=av_len(inters); counter1++){
    int seen=-1;
    char *target[20];
    SV *target_prot=*av_fetch(inters,counter1,0);
    strcpy(target,SvPV(target_prot,PL_na));
    
    /*******************************************************/
    /* Check if we have already seen this target as a bait */
    /* and, if so, skip it				   */
    /*******************************************************/
    SV **seen_bait_rref=hv_fetch( seen_prots, target , strlen(target), 0);
    if(seen_bait_rref!=NULL){
      HV *seen_bait=(HV*)SvRV(*seen_bait_rref); 
      SV **seen_target_rref=hv_fetch( seen_bait,prot  , strlen(prot), 0);
      
      /***********************************************************/
      /* If $local_seen{$bait}{$target} exists, then this 	 */
      /* bait/target pair has already been seen as target/bait.  */
      /* Unless bait==target in which case it has been seen 	 */
      /* if $local_seen{$bait}{$target}>1			 */
      /***********************************************************/
      if(seen_target_rref!=NULL){
	if(strcmp(prot,target)==0){
	  if(SvIV(*seen_target_rref)>1)
	    continue;
	}
	else
	  continue;
      }
      
    }
    /**************************************************/
    /* Now, iterate through each target protein's gos */
    /**************************************************/
    SV **target_annots_rref=hv_fetch( prot_annots, target , strlen(target), 0);
    AV *target_annots;
    if(target_annots_rref !=NULL){
      target_annots=(AV*)SvRV(*target_annots_rref); 
      for(counter2=0;counter2<=av_len(target_annots); counter2++){ 
	if(debug>0)
	  printf(" P/T:%s,%s %s(%s) want:%s seen:%d",prot,target,go,go_onto, SvPV(wanted_onto,PL_na),seen);
	if(seen==1){
	    if(debug>0)
	      printf(" SKIPPED(%s)\n",SvPV(*av_fetch(target_annots,counter2,0), PL_na) );
	    continue;
	}
	char *tgo[20], *t_onto[2];
	SV **tgo_rref=av_fetch(target_annots,counter2,0);
	if (tgo_rref != NULL){
	  SV *tgo_ref=*tgo_rref;
	  strcpy(tgo,SvPV(tgo_ref,PL_na));
	  if(debug>0)
	    printf(", TGO:%s",tgo);

	  // Get target GO ontology
	  strcpy(t_onto,SvPV(*hv_fetch( go_ontos , tgo , strlen(tgo), 0),PL_na));
	  if(debug>0)
	    printf("(%s)",t_onto);
	  /***************************************/
	  /* now sort to get oo, both ontologies */
	  /* concatenated in alphabetical order  */
	  /***************************************/
	  if(strcmp(go_onto,t_onto)<0){
	    strcpy(oo,go_onto);
	    strcat(oo,t_onto);
	  }
	  else{
	    strcpy(oo,t_onto);
	    strcat(oo,go_onto);
	  }
	  // If we want this ontology
	  if(strcmp(oo,SvPV(wanted_onto,PL_na)) == 0){
	    interactions++;
	    seen=1;
	  }
	  if(debug>0)
	    printf("\too:%s, inter:%d, seen2:%d\n",oo,interactions,seen);
	}
      } //end for counter2
    } // end if target_annots_rref !=NULL
  }//end for counter1
  fflush(stdout);
  return(interactions);
} 

  

    
















 
//    c_loop_inter3($#GOs,$cc,$verbose,$subonto,$debug,\%{$go_ontos{ONTO}},\%{$go_counts{GOs}},\%{$go_counts{ONTOs}},\%gos,\%wanted_proteins,\%prot_annot,\@inter_pairs,\@GOs);
/* SV* c_loop_inter3(int max,int inter_max,int verbose, SV* SVontology, int debug, SV* go_ontos_Oref, SV* go_counts_Gref, SV* go_counts_Oref, SV* gos_hash_ref, SV* wanted_prots_ref, SV* prot_annots_ref,SV* inter_pairs_ref, SV* GOsSV) { */
/*   HE* onto_entry, *hash_entry2; */
/*   char pair[100], go1[12], go2[12], onto[2], oo[3], wanted_onto[4], id[22]; */
/*   int onto_keys, counter1,counter2, l, go1count,go2count,both,tot_onto,b_go1count,b_go2count,b_both,b_tot_onto,want; */
/*   float o; //overrepresentation */
/*   SV *go1_onto, *go2_onto, *newhashref, *go1cnt, *go2cnt, **go1_protss , **go2_protss, *go1_prots, *go2_prots, *GOsSV_ref1, **GOsSV_ref2, **pairSV_ref2, *GO1, **GOsSV2_ref2, *GO2, *inter_pairs; */
/*   HV *go_ontos, *go_counts, *go1_prots_hash, *go2_prots_hash, *go1_count_hash, *go2_count_hash,*gos,*onto_count, *wanted_prots; */
/*   go_ontos = (HV*)SvRV(go_ontos_Oref); */
/*   go_counts = (HV*)SvRV(go_counts_Gref); */
/*   gos=(HV*)SvRV(gos_hash_ref); // get the %gos hash */
/*   inter_pairs=(AV*)SvRV(inter_pairs_ref); // get the %interactors hash */
/*   wanted_prots=(HV*)SvRV(wanted_prots_ref); // get the %wanted_proteins hash */
/*   HV *prot_annots=(HV*)SvRV(prot_annots_ref); */

 
/*   struct my_hash { */
/*     int id;                    /\* key *\/ */
/*     char name[20]; */
/*     int interactions; */
/*     UT_hash_handle hh;         /\* makes this structure hashable *\/ */
/*   }; */




/*   struct my_hash *s, *pairs = NULL;  //initialize the pairs hash */

/*   int hash_counter=0; */
/*   for(counter1=0;counter1<inter_max; counter1++){ */
/*     fprintf(stderr,"%d of %d\r",counter1,inter_max); */

/*     //get pair */
/*     pairSV_ref2=av_fetch(inter_pairs,counter1,0); */
/*     SV *pairSV=*pairSV_ref2 ; */
/*     strcpy(pair,SvPV(pairSV,PL_na)); */
/*     char *tmp_pair[100]; */
/*     char *token; */
/*     strcpy (tmp_pair,pair); */
/*     //get go1 */
/*     char *prot1, *prot2; */
/*     prot1=strtok (tmp_pair, "-"); */
/*     prot2=strtok (NULL, "\0"); */
    
/*     //Get go terms for prot1 */
/*     SV** prot1_golist_rref=hv_fetch(prot_annots , prot1 , strlen(prot1), 0);  */
/*     AV *prot1_golist; */
/*     if(prot1_golist_rref != NULL){ */
/*       prot1_golist  =(AV*)SvRV(*prot1_golist_rref); */
/*     } */
/*     int cc1; */
/*     for(cc1=0; cc1<av_len(prot1_golist); cc1++){ */
/*       int inversed=0; */
/*       SV *go1_ref=*av_fetch(prot1_golist,cc1,0); */
/*       //Get go terms for prot2 */
/*       SV** prot2_golist_rref=hv_fetch(prot_annots , prot2 , strlen(prot2), 0);  */
/*       AV *prot2_golist; */
/*       if(prot2_golist_rref != NULL){ */
/* 	prot2_golist  =(AV*)SvRV(*prot2_golist_rref); */
/*       } */
/*       int cc2; */
/*       for(cc2=0; cc2<av_len(prot2_golist); cc2++){ */
/* 	//SV **go2_rref=av_fetch(prot2_golist,cc2,0); */
/* 	SV *go2_ref=*av_fetch(prot2_golist,cc2,0); */
	
/* 	//Get pair */
/* 	if(strcmp(go1,go2)<0){ */
/* 	  strcpy(pair,go1); */
/* 	  strcat(pair,"_"); */
/* 	  strcat(pair,go2); */
/* 	} */
/* 	else{ */
/* 	  strcpy(pair,go2); */
/* 	  strcat(pair,"_"); */
/* 	  strcat(pair,go1); */
/* 	  inversed=1; */
/* 	} */
	
/* 	s = malloc(sizeof(struct my_hash)); */
/*         strcpy(s->name, pair); */
/*         s->id = hash_counter++; */
/*         HASH_ADD_STR( pairs, name, s ); */
	
	

/*       } */
/*     } */



/*     // Foreach go of prot1 */
/*        // Foreach go of prot2 */

/*   } */



/*     fflush(stdout); */
/* } */

