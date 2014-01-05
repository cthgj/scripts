#!/usr/bin/perl -w
use strict;
use Getopt::Std;
use IO::File;
use Switch;
use Term::ANSIColor; 
use Math::Round;

my (%terms_to_GOs,%not_moonlighting,%seen,%moonlighting,%ancestors,%prob,%inter_prob,%offspring,%opts,%proteins,%synonyms,%classes,%precision);
getopts('iACcdhIbgNntVvM:m:a:T:D:S:p:P:G:s:f:',\%opts) || do { print "Invalid option, try 'moonGO.pl -h' for more information\n"; exit(1); };
my (@net_probs);
my $very_verbose=$opts{V}||undef;
my $Aprob_file_read=0; ## is 1 if the annotation prob file for this network has been read, gopair_probability
my $Iprob_file_read=0; ## is 1 if the interactome prob file for this network has been read, interactome_probability
my $have_already_read_terms_file=0;
my $have_already_read_precision_file=0;
my $have_already_read_probs_file=0;
&usage() if $opts{h};
my $DATADIR=$opts{D}||"./data";
my $stats_dir="$DATADIR/gostats";
my $precision_file="$DATADIR/precision.human";
my $inter_stats_dir;
my $no_anc=$opts{C}||undef;
my $check_both=$opts{b}||undef; ## Check BOTH inter and annot probs
my $only_gold=$opts{g}||undef;
my $geneology_file=$opts{G}||"$DATADIR/biological_process.genealogy";
my $go_terms_file=$opts{T}||"$DATADIR/GO.terms_alt_ids";
my $gold_file = "$DATADIR/gold_simple.txt";
my %gold; ## information on the gold standard proteins
my $debug=$opts{d}||undef;
my $help=$opts{h}||undef;
my $check_ancestors=$opts{A}||undef;
my $min_class_cardinality=$opts{M}||0;
my $network_file=$ARGV[0]|| die("Need a network file\n");
my $network_prob_file; ## Network specific probabilities
my $low_annot_prob=$opts{p}||0.9393275;#1.193455e-07;
my $low_inter_prob=$opts{p}||0.02463311;#1.193455e-07;
my $check_interactome_probs=$opts{i}||undef;
my $subonto=$opts{o}||"P";
my $verbose=$opts{v}||undef;
my $no_class=$opts{n}||undef; 
my $synfile;
my $species=$opts{S}||'human';
my $gaf_annotations_file;
## $added_to_prob_file will be 1 if a 
## probability is added to this network's .prob file
my $added_to_prob_file=0; 
my $suffix; ## Used to add extension to a name. Eg, APBB3 => APBB3_MOUSE
my $prob; ## probability significance threshold, either ($low_annot_prob or $low_inter_prob)
          ## see sub parse_input_files_guess_values
my $read_network_prob_file=1;
$opts{N} && do {$read_network_prob_file=undef};
my $TESTS=0;
$verbose=1 if $debug;

my $mode=$opts{m}|| 'i';
$synonyms{LOADED}=0;




############################### Main program ###############################
&parse_input_files_guess_values();


# my $aa=&get_name("P06744","bob");
# my $bb=&get_name("G6PI_HUMAN","bob");
# print "aa :$aa\nbb : $bb\n"; die();

my $interactome_prob_file=$opts{I}||"$stats_dir/$species.interactome.probs";
my $prob_file=$stats_dir . "/$species.prob" || $stats_dir . "/human.prob";
my $annotations_file=$opts{a}||"$DATADIR/$species.clas";
my $date=localtime;
my $LOG = IO::File->new(">> $DATADIR/moonlighting.log")|| die("Cannot open $DATADIR/moonlighting.log : $!\n");
print $LOG "********** $date ******************\n";
my $com=`ps aux | grep $0 | head -1`;
$com=~s/^.+:\d+\s*(.+)\n/$1/;
print $LOG "$com\n";
print $LOG "Network     : $network_file\n";
print $LOG "Annotations : $annotations_file\n" if $annotations_file;
print $LOG "GAF file    : $gaf_annotations_file\n";
print $LOG "Mode        : $opts{m}\n";
print $LOG "Stats dir   : $stats_dir\n";
print $LOG "Species     : $species\n";
print $LOG "Threshold an: $low_annot_prob\n";
print $LOG "Threshold in: $low_inter_prob\n";
print $LOG "Prob        : $prob\n";
print $LOG "\n";
$LOG->close;
 

&load_network($network_file); 
&load_ancestors();
&parse_gaf_file();
&load_annotations($annotations_file) if $annotations_file;
&check_gold(); ## read the list of gold standard prots
my $protnum=scalar(keys(%proteins));
my $c=0; ## counter


# bait:foreach my $bait (keys(%proteins)){
#     my @a=keys(%{$proteins{$bait}{GOs}});
#     print "$bait : @a\n";
# }
# die();

switch($mode) {
    ## Count how many interactions involve at least one pair of distant GOs
    case ['count']{
	my $num=0;
      bait:foreach my $bait (keys(%proteins)){
	  $c++;
	  printf STDERR ("$c of $protnum\r") if $verbose;
	target:foreach my $target (keys(%{$proteins{$bait}{INTERACTORS}})){
	    foreach my $bgo (keys(%{$proteins{$bait}{GOs}})){
		next bgo if $bgo eq 'GO:0008150';
		tgo:foreach my $tgo (keys(%{$proteins{$target}{GOs}})){
		    next tgo if $tgo eq 'GO:0008150';
		    my $P=&gopair_probability($bgo,$tgo);
		    if ($P<=$low_annot_prob){
			$num++;
			next target;
		    }
		}
	    }
	}
      }
	print "$num prots";
    }
    ## Count how many interactions involve UNLIKELY and how many
    ## involve LIKELY pairs
    case ['likeunlike']{
	my $num=0;
	my (%likely,%unlikely);
	print STDERR "Reading Likely..." if $verbose;
	my ($ll,$lu)=(1,1);
	open(L,"$ARGV[1]")||die("Need a list of likely pairs as \$ARGV[1]\n");
	while(<L>){
	    /^(.+?)\s.+\t(.+?)$/;
	    $ll=$2 if $2<$ll;
	    $likely{$1}++
	}
	close(L);
	print STDERR "Done\nReading Unlikely..." if $verbose;
      	open(U,"$ARGV[2]")||die("Need a list of unlikely pairs as \$ARGV[2]\n");
	while(<U>){
	    /^(.+?)\s.+\t(.+?)$/;
	    $lu=$2 if $2<$lu;
	    $unlikely{$1}++
	}
	close(U);
	print STDERR "Done (likely: $ll, unlikely: $lu\n" if $verbose;
	my %seen;
	my $tests=0;
	$likely{PROTS}=0;
	$unlikely{PROTS}=0;
	my %gocount;
      bait:foreach my $bait (keys(%proteins)){
	  $c++;
	  printf STDERR ("$c of $protnum\t $likely{PROTS}/$unlikely{PROTS}/$tests\r") if $verbose;
	target:foreach my $target (keys(%{$proteins{$bait}{INTERACTORS}})){
	    my @pp=sort {$b lt $a} ($bait,$target);
	    my $pair=join("_",@pp);
	    if (defined($seen{$pair})){
		print STDERR "\nSkipped $bait $target\n";
		next target;
	    }
	  bgo:foreach my $bgo (keys(%{$proteins{$bait}{GOs}})){
		next bgo if $bgo eq 'GO:0008150';
	      tgo:foreach my $tgo (keys(%{$proteins{$target}{GOs}})){
		  next tgo if $tgo eq 'GO:0008150';
		  $tests++;
		  my @b = sort {$b lt $a} ($tgo,$bgo);
		  my $gopair=join("_",@b);
		  $gocount{$gopair}++;
		  if(defined($likely{$gopair})){
		      $likely{PROTS}++;
		      die("Crap1 $gopair\n") if defined($unlikely{$gopair});
		  }
		  elsif(defined($unlikely{$gopair})){
		      $unlikely{PROTS}++;
		      die("Crap1 $gopair\n") if defined($likely{$gopair});

		  }
	      }
	    }
	}
      }
	print STDERR "\n" if $verbose;
	map{print "$_\t$gocount{$_}\n"}keys(%gocount);
	print "TOTAL    : $tests\n";
	print "LIKELY   : $likely{PROTS}\n";
	print "UNLIKELY : $unlikely{PROTS}\n";

    }
############################################################################
########### Look for proteins with at least 1 DIRECT annotation  ###########
###########          dissimilar to ALL those of the class        ###########
############################################################################
    case ['a','dis2class'] {
	my %candidates;
	my $class_prob=$low_annot_prob * (10**7);
      bait:foreach my $bait (keys(%proteins)){
	  $c++;
	  printf STDERR ("$c of $protnum\r") if $verbose;
	  my @silly_tmp_array=keys(%{$proteins{$bait}{CLASS}});
	  if(defined($gold{$bait})){push@{$gold{$bait}{MISSED}{$mode}},"UnClassed" if $#silly_tmp_array==-1;}
	class:foreach my $class (keys(%{$proteins{$bait}{CLASS}})){
	    my $is_cand=10;
	    ## Class annotations are redundant
	    my %seen_annots=();
	    my %hits;
	    my @uniqu = grep { ! $seen_annots{$_} ++ } @{$classes{$class}{ANNOT}};
	  c_anot:foreach my $class_annot (@uniqu){
	      next c_anot if $is_cand==0;
	      next c_anot if $class_annot eq 'GO:0008150';
	    bgo:foreach my $bgo (keys(%{$proteins{$bait}{GOs}})){
		next bgo if $is_cand==0;
		next bgo if $bgo eq 'GO:0008150';
		## 'INF' are the class annotations
		next bgo if $proteins{$bait}{GOs}{$bgo} eq "INF";
		## 'UNC' are the 2ndary class annotations
		next bgo if $proteins{$bait}{GOs}{$bgo} eq "UNC";
		#next bgo if $proteins{$bait}{GOs}{$bgo} eq "ANC";
		#next bgo unless $proteins{$bait}{GOs}{$bgo} eq "DIR";
		my $P;
		## If we are checking BOTH annot and inter probs
		if($check_both){
		    my ($Pi,$Pa)=(&interactome_probability($bgo,$class_annot),&gopair_probability($bgo,$class_annot));	
		    if($Pa<=$prob && $Pi<=$prob){
			$is_cand=1;
			## Make sure we keep the LEAST probable of the interactions
			## in question
			if(defined($hits{$bait}{$class})){
			    $hits{$bait}{$class}=~/.+\t:\s(.+?)\t(.+?)$/;
			    my $ka=$1;
			    my $ki=$2;
			    ## If we have already found this bait as a cand,
			    ## keep the best combination of pvals (annot and inter)
			    if($Pi<=$ki && $Pa<=$ka){
				$hits{$bait}{$class}="$bait\t($bgo)  \t:\t$class, $class_annot\t: $Pi\t$Pa";
			    }
			}
			## If this is the first time we find this bait as a cand
			else{
			    $hits{$bait}{$class}="$bait\t($bgo)  \t:\t$class, $class_annot\t: $Pi\t$Pa";
			}
		    }
		    else{
			$is_cand=0;
			## If this prot is in the gold standard,
			## check WHY it was discarded1
			&check_gold($mode,$bait,$bgo,$class_annot,$Pi . "_" . $Pa,$class);
		    }
		}
		else{
		    $check_interactome_probs ?   
			($P=&interactome_probability($bgo,$class_annot)) : 
			($P=&gopair_probability($bgo,$class_annot));		
		    if ($P<=$prob){
			$is_cand=1;
			## Make sure we keep the LEAST probable of the interactions
			## in question
			if(defined($hits{$bait}{$class})){
			    $hits{$bait}{$class}=~/.+\t:\s(.+?)$/;
			    my $koko=$1;
			    $hits{$bait}{$class}="$bait\t($bgo)  \t:\t$class, $class_annot\t: $P" if $P<$koko;
			}
			else{
			    $hits{$bait}{$class}="$bait\t($bgo)  \t:\t$class, $class_annot\t: $P";
			}
		    }
		    else{
			$is_cand=0;
			## If this prot is in the gold standard,
			## check WHY it was discarded1
			&check_gold($mode,$bait,$bgo,$class_annot,$P,$class);
		    }
		}
	    } ## end foreach $bgo
	  }## end foreach $class_annot
	    if($is_cand==1){
		$moonlighting{$bait}{$class}=$hits{$bait}{$class};
	    }
	} ## end foreach class
	  # unless($is_cand==0){
	  # 	  $candidates{$bait}=1;
	  # 	  print "$bait\t$class\n";
	  # }
	  # else{die();
	  # 	  $candidates{$bait}=1;
	  # 	  print "$bait\t$class\n";
	  # 	  next bait; 
	  # }
	  
	  
      } ## end foreach bait
	print STDERR "\n" if $verbose;
	## Now look for improb pairs within the candidates
	$c=0;
	# bait:foreach my $bait (keys(%candidates)){
	#     $c++;
	#     my $is_cand=0;
	#     my $oo=scalar(keys(%candidates));
	#     printf STDERR ("$c of $oo (round 2)\r") if $verbose;
	#   target:foreach my $target (keys(%{$proteins{$bait}{INTERACTORS}})){
	#     tgo:foreach my $tgo (keys(%{$proteins{$target}{GOs}})){
	# 	next tgo if $tgo eq 'GO:0008150';
	# 	my %bait_gos;
	# 	map{$bait_gos{$_}=$_}keys(%{$proteins{$bait}{GOs}});
	#       bgo:foreach my $bgo (keys(%bait_gos)){
	# 	  next bgo if $bgo eq 'GO:0008150';
	# 	  my $P=&gopair_probability($bgo,$tgo,"direct");
	# 	  if ($P<=$low_prob){
	# 	       $moonlighting{$bait}{$target}="$bait ($bgo)  \t:\t$target ($tgo)  \t: $P";
	# 	       next target;
	# 	  }
	# 	   else{
	# 	       $moonlighting{$bait}{$target}=undef;
	# 	       $not_moonlighting{$bait}{$target}++;
	# 	      next target;
	# 	   }
	#       }
	#     }
	#   }
	# }		     


    } ## end case 'a'
############################################################################
########### Look for interactions bridging dissimilar classes    ###########
###########                                                      ###########
############################################################################
   # case ['b','bridge'] {
   #    bait:foreach my $bait (keys(%proteins)){
   # 	  $c++;
   # 	  my $is_cand=0;
   # 	  printf STDERR ("$c of $protnum\r") if $verbose;
   # 	  my @targets=keys(%{$proteins{$bait}{INTERACTORS}});
   # 	  my @bait_classes=keys(%{$proteins{$bait}{CLASS}});
   # 	  if(defined($gold{$bait})){push@{$gold{$bait}{MISSED}{$mode}},"UnClassed" if $#bait_classes==-1;}
   # 	t1:foreach  my $target1(keys(%{$proteins{$bait}{INTERACTORS}})){
   # 	    ## next t1 if it shares a class with bait
   # 	    map{next t1 if defined($proteins{$bait}{CLASS}{$_});}keys(%{$proteins{$target1}{CLASS}});
   # 	    t2:foreach  my $target2(keys(%{$proteins{$bait}{INTERACTORS}})){
   # 		## next t2 if it shares a class with bait
   # 		map{next t2 if defined($proteins{$bait}{CLASS}{$_});}keys(%{$proteins{$target2}{CLASS}});
   # 		## next t2 if t1 shares class with t2
   # 		map{next t2 if defined($proteins{$target1}{CLASS}{$_});}keys(%{$proteins{$target2}{CLASS}});
   # 	      class1:foreach my $class1 (keys(%{$proteins{$target1}{CLASS}})){ 
   # 		  class2:foreach my $class2 (keys(%{$proteins{$target2}{CLASS}})){ 
   # 		      ## Class annotations are redundant
   # 		      my %seen_annots1=();
   # 		      my @uniqu1 = grep { ! $seen_annots1{$_} ++ } @{$classes{$class1}{ANNOT}};
   # 		      my %seen_annots2=();
   # 		      my @uniqu2 = grep { ! $seen_annots2{$_} ++ } @{$classes{$class2}{ANNOT}};
   # 		      c1:foreach my $class1_annot (@uniqu1){
   # 			  next c1 if $class1_annot eq 'GO:0008150';
   # 			  my $prec1=&term_precision($class1_annot);
   # 			  c2:foreach my $class2_annot (@uniqu2){
   # 			      next c2 if $class2_annot eq 'GO:0008150';
   # 			      my $prec2=&term_precision($class2_annot);
   # 			      my $P;
   # 			      my @c= sort {$b lt $a} ($class1,$class2);
   # 			      my $class_name=$c[0] . $c[1];
   # 			      $check_interactome_probs ?   
   # 				  ($P=&interactome_probability($class1_annot,$class2_annot)) : 
   # 				  ($P=&gopair_probability($class1_annot,$class2_annot));
   # 			      if ($P<=$prob){
   # 				  ## Make sure we keep the MOST specific LEAST probable of the interactions
   # 				  ## in question	 
   # 				  if(defined($moonlighting{$bait}{$class_name})){
   # 				      $moonlighting{$bait}{$class_name}=~/^.+?\t.+?\s.+?([\.\d]+)\t.+?\s.+?([\.\d]+).+:\s(.+?)$/;
   # 				      my ($p1,$p2,$pp)=($1,$2,$3);
   # 				      # $moonlighting{$bait}=~/^.+?\t.+?\t.+?\s.+?([\.\d]+)/;
   # 				      # my $p2=$1;
   # 				      # $moonlighting{$bait}=~/.+\t:\s(.+?)$/;
   # 				      if($p1+$p2 > $prec1+$prec2){
   # 					  $moonlighting{$bait}{$class_name}="$bait\t$class1 $class1_annot $prec1\t$class2 $class2_annot $prec2\t: $P"; 
   # 				      }
   # 				      elsif($p1+$p2 == $prec1+$prec2){
   # 					  $moonlighting{$bait}{$class_name}="$bait\t$class1 $class1_annot $prec1\t$class2 $class2_annot $prec2\t: $P" if $pp<=$P;   
   # 				      }
   # 				  }
   # 				  else{
   # 				      $moonlighting{$bait}{$class_name}="$bait\t$class1 $class1_annot $prec1\t$class2 $class2_annot $prec2\t: $P";
   # 				  }
   # 				  #next bait; 
   # 			      }
   # 			      else{
   # 				  my $cc=$class1 . "-" . $class2;
   # 				  &check_gold($mode,$bait,$class1_annot,$class2_annot,$P,$cc);
   # 				  next class2
   # 			      }
   # 			  } ## end c2
   # 		      } ## end c1
   # 		  } ## end class2
   # 		} ## end class1
   # 	      } ## end t2
   # 	  } ## end t2
   #    } ## end foreach bait
   # } ## end case 'b'


############################################################################
########### Look for interactions bridging dissimilar classes    ###########
############################################################################
   case ['B','unique_bridge','b','bridge'] {
     bait:foreach my $bait (keys(%proteins)){
	 my %gg;
	  $c++;
	  my $is_cand=0;
	  printf STDERR ("$c of $protnum\r") if $verbose;
	  my @targets=keys(%{$proteins{$bait}{INTERACTORS}});
	  my @bait_classes=keys(%{$proteins{$bait}{CLASS}});
	  if(defined($gold{$bait})){push@{$gold{$bait}{MISSED}{$mode}},"UnClassed" if $#bait_classes==-1;}
	t1:foreach  my $target1(keys(%{$proteins{$bait}{INTERACTORS}})){
	    ## next t1 if it shares a class with bait
	    map{push@{$gold{$bait}{MISSED}{$mode}},"Shared_Class $target1 $_" if defined($gold{$bait}); next t1 if defined($proteins{$bait}{CLASS}{$_});}keys(%{$proteins{$target1}{CLASS}});
	    t2:foreach  my $target2(keys(%{$proteins{$bait}{INTERACTORS}})){
		## next t2 if it shares a class with bait
		map{push@{$gold{$bait}{MISSED}{$mode}},"Shared_Class $target2 $_" if defined($gold{$bait}); next t2 if defined($proteins{$bait}{CLASS}{$_});}keys(%{$proteins{$target2}{CLASS}});
		## next t2 if t1 shares class with t2
		map{push@{$gold{$bait}{MISSED}{$mode}},"Shared_Class $target1-$target2 $_" if defined($gold{$bait}); next t2 if defined($proteins{$target1}{CLASS}{$_});}keys(%{$proteins{$target2}{CLASS}});
		class1:foreach my $class1 (keys(%{$proteins{$target1}{CLASS}})){ 
		    ## Skip small classes
		    next class1 unless scalar(keys(%{$classes{$class1}{PROTS}}))>$min_class_cardinality;
		  class2:foreach my $class2 (keys(%{$proteins{$target2}{CLASS}})){  
		      ## Skip small classes
		      next class2 unless scalar(keys(%{$classes{$class2}{PROTS}}))>$min_class_cardinality;
		      ## Class annotations are redundant
		     
		      my %seen_annots1=();
		      my @uniqu1 = grep { ! $seen_annots1{$_} ++ } @{$classes{$class1}{ANNOT}};
		      my %seen_annots2=();
		      my @uniqu2 = grep { ! $seen_annots2{$_} ++ } @{$classes{$class2}{ANNOT}};
		    c1:foreach my $class1_annot (@uniqu1){
			  next c1 if $class1_annot eq 'GO:0008150';
			  my $prec1=&term_precision($class1_annot);
			c2:foreach my $class2_annot (@uniqu2){	
			    next c2 if $class2_annot eq 'GO:0008150';
			    my $prec2=&term_precision($class2_annot);
			    my $P;
			    my @c= sort {$b lt $a} ($class1,$class2);
			    my $class_name=$c[0] .  "_" . $c[1];
			    ## If we are checking BOTH annot and inter probs
			    if($check_both){
				my ($Pi,$Pa)=(&interactome_probability($class1_annot,$class2_annot),&gopair_probability($class1_annot,$class2_annot));	
				if($Pa<=$prob && $Pi<=$prob){
				    ## Make sure we keep the MOST specific LEAST probable of the interactions
				    ## in question
				    if(defined($moonlighting{$bait}{$class_name})){
					$moonlighting{$bait}{$class_name}=~/^.+?\t.+?\s.+?([\.\d]+)\t.+?\s.+?([\.\d]+).+:\s*(.+?)\t(.+?)$/;
					##prec1,prec2,inter_prob,annot_prob
					my ($p1,$p2,$pi,$pa)=($1,$2,$3,$4);
					if($p1+$p2 > $prec1+$prec2){
					    $moonlighting{$bait}{$class_name}="$bait\t$class1 " . 
						scalar(keys(%{$classes{$class1}{PROTS}})) . 
						" $class1_annot $prec1\t$class2 " . 
						scalar(keys(%{$classes{$class2}{PROTS}})) . 
						" $class2_annot $prec2\t: $Pi\t$Pa"; 
					}
					## If it is AS precise, keep if less likely
					elsif($p1+$p2 == $prec1+$prec2){
					    $moonlighting{$bait}{$class_name}="$bait\t$class1 " . 
						scalar(keys(%{$classes{$class1}{PROTS}})) . 
						" $class1_annot $prec1\t$class2 " . 
						scalar(keys(%{$classes{$class2}{PROTS}})) . 
						" $class2_annot $prec2\t: $Pi\t$Pa" if $pa<=$P;   
					}
				    }
				    else{
					$moonlighting{$bait}{$class_name}="$bait\t$class1 " . 
					    scalar(keys(%{$classes{$class1}{PROTS}})) . 
					    " $class1_annot $prec1\t$class2 " . 
					    scalar(keys(%{$classes{$class2}{PROTS}})) . 
					    " $class2_annot $prec2\t: $Pi\t$Pa";
				    }
				}
				else{
				    my $cc=$class1 . "-" . $class2;
				    &check_gold($mode,$bait,$class1_annot,$class2_annot,$Pi . "_" . $Pa,$cc);
				    next class2
				}
			    }
			    ## If we are checking EITHER annot OR inter probs
			    else{
				$check_interactome_probs ?   
				    ($P=&interactome_probability($class1_annot,$class2_annot)) : 
				    ($P=&gopair_probability($class1_annot,$class2_annot));
				&debug("Checking: $bait,$class1,$class1_annot,$class2,$class2_annot,$P");
				if ($P<=$prob) {
				    ## We are not interested if the candidate is annotated to GOs that
				    ## are similar to at least one GO of BOTH interactors.
				    # bgo:foreach my $bgo (keys(%{$proteins{$bait}{GOs}})){
				    # 	  next  if $bgo eq 'GO:0008150';
				    # 	tgo:foreach my $tgo (keys(%{$proteins{$target1}{GOs}})){
				    # 	    next  if $tgo eq 'GO:0008150';
				    # 	    my $PP=0;
				    # 	    $check_interactome_probs ?   
				    # 		($PP=&interactome_probability($bgo,$tgo)) : 
				    # 		($PP=&gopair_probability($bgo,$tgo));
				    # 	    if ($PP>$prob){
				    # 		next t1;
				    # 	    }
				    # 	}
				    # 	tgo:foreach my $tgo (keys(%{$proteins{$target1}{GOs}})){
				    # 	    next  if $tgo eq 'GO:0008150';
				    # 	    my $PP=0;
				    # 	     $check_interactome_probs ?   
				    # 		 ($PP=&interactome_probability($bgo,$tgo)) : 
				    # 		 ($PP=&gopair_probability($bgo,$tgo));
				    # 	    if ($PP>$prob){
				    # 		next t2;
				    # 	    }
				    # 	}
				    # }
				    ## Make sure we keep the MOST specific LEAST probable of the interactions
				    ## in question
				    if(defined($moonlighting{$bait}{$class_name})){
#				      print "xxx1 $bait ($target1,$target2,\t$class1\t$class1_annot\t$prec1\t$class2\t$class2_annot\t$prec2\t: $P\n";
					$moonlighting{$bait}{$class_name}=~/^.+?\t.+?\s.+?([\.\d]+)\t.+?\s.+?([\.\d]+).+:\s(.+?)$/;
					my ($p1,$p2,$pp)=($1,$2,$3);
					if($p1+$p2 > $prec1+$prec2){
					    $moonlighting{$bait}{$class_name}="$bait\t$class1 " . 
						scalar(keys(%{$classes{$class1}{PROTS}})) . 
						" $class1_annot $prec1\t$class2 " . 
						scalar(keys(%{$classes{$class2}{PROTS}})) . 
						" $class2_annot $prec2\t: $P"; 
					}
					elsif($p1+$p2 == $prec1+$prec2){
					    $moonlighting{$bait}{$class_name}="$bait\t$class1 " . 
						scalar(keys(%{$classes{$class1}{PROTS}})) . 
						" $class1_annot $prec1\t$class2 " . 
						scalar(keys(%{$classes{$class2}{PROTS}})) . 
						" $class2_annot $prec2\t: $P" if $pp<=$P;   
					}
				    }
				    else{
					$moonlighting{$bait}{$class_name}="$bait\t$class1 " . 
					    scalar(keys(%{$classes{$class1}{PROTS}})) . 
					    " $class1_annot $prec1\t$class2 " . 
					    scalar(keys(%{$classes{$class2}{PROTS}})) . 
					    " $class2_annot $prec2\t: $P";
				    }
				    #next bait; 
				}
				else{
				    my $cc=$class1 . "-" . $class2;
				    &check_gold($mode,$bait,$class1_annot,$class2_annot,$P,$cc);
				    next class2
				}
			    }
			} # end class2_annot
		    } #end class1_annot
		  } #end class2
		} ## end class1 
	    } #end t2
	} ## end t1
     } ## end foreach bait
## Now find those candidates that connect classes 
## that are NOT connected by any other interactions
	 if($mode eq 'B'){
	  cand:foreach my $cand (keys(%moonlighting)){
	      classname:foreach my $classname (keys(%{$moonlighting{$cand}})){
		  $moonlighting{$cand}{$classname}=~/^$cand\t(.+?)\s+(.+?)\s+(.+?)\t(.+?)\s+(.+?)\s+(.+?)\t:\s*(.+)$/||die("Could not match cand : \$moonlighting{$cand}{$classname} : $moonlighting{$cand}{$classname}\n");
		  my ($class1,$class1_annot,$prec1,$class2,$class2_annot,$prec2,$P)=($1,$2,$3,$4,$5);
		  ## Make sure there are no connections
		  ## between any prot in class1 and any in class2
		  my @prots1= keys(%{$classes{$class1}{PROTS}}) ;
		  my @prots2= keys(%{$classes{$class2}{PROTS}}) ;
		p1:foreach my $prot1 (@prots1){
		  p2:foreach my $prot2 (@prots2){
		      if(defined($proteins{$prot1}{INTERACTORS}{$prot2})){
			  push@{$gold{$cand}{MISSED}{$mode}},"$class1\t$class2\tConnected" if defined($gold{$cand});
			  $not_moonlighting{$cand}{$classname}++;
			  next classname;
		      }
		  }
		}
	      } ## end foreach $classname
	  } ## end foreach $cand
	}
   } #end case B/b

############################################################################
###########       Look for improbable interactions between       ###########
########## proteins annotated to different, (dissimilar) classes  ##########
############################################################################
   case ['c','difclass'] {
       print STDERR "\nWARNING: not sure about this one, check source\n\n";
      bait:foreach my $bait (keys(%proteins)){
	  $c++;
	  printf STDERR ("$c of $protnum\r") if $verbose;
	target:foreach my $target (keys(%{$proteins{$bait}{INTERACTORS}})){
	    my $possible_cand=1;
	    ## Skip if we have already seen this
	    if( defined($not_moonlighting{$bait}{$target}) ||
		defined($moonlighting{$bait}{$target}) ||
		defined($moonlighting{$target}{$bait}) ||
		defined($not_moonlighting{$target}{$bait})){
		next target;
	    }
	    my @bait_classes=keys(%{$proteins{$bait}{CLASS}});
	    ##skip target if it shares a class with bait;
	    foreach my $class (@bait_classes){
		next target if defined($proteins{$target}{CLASS}{$class});
	    }
 	    ## So, target shares no class with bait. Now
            ## only check if ALL class annotations are dissimilar
	    # foreach my $bclass (@bait_classes){
	    # 	foreach my $bclass_annot (@{$classes{$bclass}{ANNOT}}){
	    # 	    next if $bclass_annot eq 'GO:0008150';
	    # 	    foreach my $tclass (keys(%{$proteins{$target}{CLASS}})){
	    # 		foreach my $tclass_annot (@{$classes{$tclass}{ANNOT}}){
	    # 		    next if $tclass_annot eq 'GO:0008150';
	    # 		    my $P=&gopair_probability($bclass_annot,$tclass_annot);
	    # 		    ## If the class annotations are disimilar,
	    # 		    ## keep the bait:target pair
	    # 		    if ($P<=$low_prob){
	    # 			$possible_cand=1;
	    # 			$moonlighting{$bait}{$target}="$bait ($bclass:$bclass_annot)  \t:\t$target ($tclass:$tclass_annot)  \t: $P";
	    # 		    }
	    # 		    ## If even one pair of bait/target classes 
	    # 		    ## are similar, next target
	    # 		    else{
	    # 			$possible_cand=0;
	    # 			$not_moonlighting{$bait}{$target}=1;
	    # 			$moonlighting{$bait}{$target}=undef;
	    # 			next target;
	    # 		    }
	    # 		}
	    # 	    }
	    # 	}
	    # }
	    ## If ALL pairs of bait/target classes are dissimilar,
	    ## check the interacting pair
	   # my $is_cand=0;
	    #if($possible_cand==1){

	    	foreach my $bgo (keys(%{$proteins{$bait}{GOs}})){
	    	    next if $bgo eq 'GO:0008150';
	    	    next unless $proteins{$bait}{GOs}{$bgo} eq 'DIR'; 
	    	    foreach my $tgo (keys(%{$proteins{$target}{GOs}})){
	    		next if $tgo eq 'GO:0008150';
	    		next unless $proteins{$target}{GOs}{$tgo} eq 'DIR'; 
	    		my $P;
			## If we are checking BOTH annot and inter probs
			if($check_both){
			    my ($Pi,$Pa)=(&interactome_probability($bgo,$tgo),&gopair_probability($bgo,$tgo));	
			    if($Pa<=$prob && $Pi<=$prob){
				## Make sure we keep the MOST specific LEAST probable of the interactions
				## in question
				if(defined($moonlighting{$bait}{$target})){
				    $moonlighting{$bait}{$target}=~/.+\t:\s*(.+?)\t(.+?)$/;
				    my ($ka, $ki)=($1,$2);
				    $moonlighting{$bait}{$target}="$bait ($bgo)  \t:\t$target ($tgo)  \t: $Pi\t$Pa" if $Pa<$ka;
				}
				else{
				    $moonlighting{$bait}{$target}="$bait ($bgo)  \t:\t$target ($tgo)  \t: $Pi\t$Pa";
				}
			    }
			    else{
				$not_moonlighting{$bait}{$target}=1;
				&check_gold($mode,$bait,$bgo,$tgo,$Pi . "_" . $Pa);
				$moonlighting{$bait}{$target}=undef;
				next target;
			    }
			}
			else{
			    $check_interactome_probs ?   
				($P=&interactome_probability($bgo,$tgo)) : 
				($P=&gopair_probability($bgo,$tgo));
			    
			    if ($P<=$prob){
				## Make sure we keep the LEAST probable of the interactions
				## in question	 
				if(defined($moonlighting{$bait}{$target})){
				    $moonlighting{$bait}{$target}=~/.+\t:\s(.+?)$/;
				    my $koko=$1;
				    $moonlighting{$bait}{$target}="$bait ($bgo)  \t:\t$target ($tgo)  \t: $P" if $P<$koko;
				}
				else{
				    $moonlighting{$bait}{$target}="$bait ($bgo)  \t:\t$target ($tgo)  \t: $P";
				}
			    }
			    else{
				$not_moonlighting{$bait}{$target}=1;
				&check_gold($mode,$bait,$bgo,$tgo,$P);
				$moonlighting{$bait}{$target}=undef;
				next target;
			    }
			}
	    	    }
	    	}
	    #}
	    	    
	   #  my $is_cand=0;
# 	    foreach my $bgo (keys(%{$proteins{$bait}{GOs}})){
# 		next if $bgo eq 'GO:0008150';
# 		next if $proteins{$bait}{GOs}{$bgo} eq 'INF'; 
# 		foreach my $tgo (keys(%{$proteins{$target}{GOs}})){
# 		    next if $tgo eq 'GO:0008150';
# #		    print STDERR "$target : $tgo : $proteins{$bait}{GOs}{$tgo}\n";
# 		    next if $proteins{$target}{GOs}{$tgo} eq 'INF'; 
# 		    my $P=&gopair_probability($bgo,$tgo);
# 		    if ($P<=$low_prob){
# 			## Check ancestor GOs
# 			if($check_ancestors){
# 			    if(&check_ancestor_gos($bgo,$tgo,$bait,$target)==1){
# 				$is_cand=1;
# 			    }
# 			    else{
# 				$is_cand=0;
# 			    }
# 			} 
#                         ## If we are NOT checking ancestors, keep 
# 			else{$is_cand=1;}
# 		    }
		    # else{
		    # 	$not_moonlighting{$bait}{$target}=1;
		    # 	$moonlighting{$bait}{$target}=undef;
		    # 	next target;
		    # }
		    # if($is_cand==1){
		    # 	$moonlighting{$bait}{$target}="$bait ($bgo)  \t:\t$target ($tgo)  \t: $P";
		    # 	next bait;
		    # }
		    # else{
		    # 	$not_moonlighting{$bait}{$target}=1;
		    # 	$moonlighting{$bait}{$target}=undef;
		    # }
		#}
	    #}
	} ## end foreach target
      } ## end foreach bait
   } ## end case "c"

############################################################################
###########       Identify ALL improbable interactions           ###########
############################################################################
    case ['d','direct'] {
     bait:foreach my $bait (keys(%proteins)){
	 $c++;
	 printf STDERR ("$c of $protnum\r") if $verbose;
	 &check_gold($bait);
	 my %bait_gos;
	 map{$bait_gos{$_}=$_}keys(%{$proteins{$bait}{GOs}});
       target:foreach my $target (keys(%{$proteins{$bait}{INTERACTORS}})){
	   if( defined($not_moonlighting{$bait}{$target}) ||
	       defined($moonlighting{$bait}{$target}) ||
	       defined($moonlighting{$target}{$bait}) ||
	       defined($not_moonlighting{$target}{$bait})){
	       next target
	   }
	   #print "$target : " . keys(%{$proteins{$target}{GOs}}) . "\n";
	 tgo:foreach my $tgo (keys(%{$proteins{$target}{GOs}})){
	       ## Useful for counting total number of interactions
	       ## involving dissimilar GOs. In that case, we are only
	       ## interested in DIRECT/INHERITED annotations, no INFERRED
	     if($no_class){next tgo unless $proteins{$target}{GOs}{$tgo} eq "DIR";}
	     
	     next tgo if $tgo eq 'GO:0008150';
	     ## skip target if the bait protein has even one of the target's GOs
	     next target if exists($bait_gos{$tgo});
	   bgo:foreach my $bgo (keys(%bait_gos)){
	       if($no_class){next bgo unless $proteins{$bait}{GOs}{$bgo} eq "DIR";}
	       next bgo if $bgo eq 'GO:0008150';
	       my $P;
	       ## If we are checking BOTH annot and inter probs
		if($check_both){
		    my ($Pi,$Pa)=(&interactome_probability($bgo,$tgo,"direct"),&gopair_probability($bgo,$tgo,"direct"));	
		    ## Make sure we keep the MOST specific LEAST probable of the interactions
		    ## in question	
		    if($Pa<=$prob && $Pi<=$prob){
			if($check_ancestors){
			    if (&check_ancestor_gos($bgo,$tgo,$bait,$target) == 1){
				$moonlighting{$bait}{$target}="$bait ($bgo)  \t:\t$target ($tgo)  \t: $Pi\t$Pa"; 
			    }
			## If one of the tgos is related to one of the bgos, skip target
			    else{$not_moonlighting{$bait}{$target}=1; next target}
			}
			else{
			    ## Make sure we keep the LEAST probable of the interactions
			    ## in question	 
			    if(defined($moonlighting{$bait}{$target})){
				$moonlighting{$bait}{$target}=~/.+\t:\s*(.+?)\t(.+?)$/;
				my ($ka, $ki)=($1,$2);
				$moonlighting{$bait}{$target}="$bait ($bgo)  \t:\t$target ($tgo)  \t: $Pi\t$Pa" if $Pa<$ka;
			    }
			    else{
				$moonlighting{$bait}{$target}="$bait ($bgo)  \t:\t$target ($tgo)  \t: $Pi\t$Pa";
			    }
			}
			next target;
		    }
		    else{
			&check_gold($mode,$bait,$bgo,$tgo,$Pi . "_" . $Pa);
			$not_moonlighting{$bait}{$target}=1;
			next target;
		    }
		}
	       ## If we are NOT checking BOTH inter and annot probs
	       else{
		   $check_interactome_probs ?   
		       ($P=&interactome_probability($bgo,$tgo,"direct")) : 
		       ($P=&gopair_probability($bgo,$tgo,"direct"));
		   if ($P<=$prob){
		       if($check_ancestors){
			   if (&check_ancestor_gos($bgo,$tgo,$bait,$target) == 1){
			       $moonlighting{$bait}{$target}="$bait ($bgo)  \t:\t$target ($tgo)  \t: $P"; 
			   }
			   ## If one of the tgos is related to one of the bgos, skip target
			   else{$not_moonlighting{$bait}{$target}=1; next target}
		       }
		       else{
			   ## Make sure we keep the LEAST probable of the interactions
			   ## in question	 
			   if(defined($moonlighting{$bait}{$target})){
			       $moonlighting{$bait}{$target}=~/.+\t:\s(.+?)$/;
			       my $koko=$1;
			       $moonlighting{$bait}{$target}="$bait ($bgo)  \t:\t$target ($tgo)  \t: $P" if $P<$koko;
			   }
			   else{
			       $moonlighting{$bait}{$target}="$bait ($bgo)  \t:\t$target ($tgo)  \t: $P";
			   }
		       }
		       next target;
		   }
		   else{
		       &check_gold($mode,$bait,$bgo,$tgo,$P);
		       $not_moonlighting{$bait}{$target}=1;
		       next target;
		   }
	       } ## end if/else check_both
	       
	   } ## end foreach my $bgo
	 } ## end foreach my $tgo
       } ## end foreach target
     } ## end foreach bait
   } ## end case "d"
############################################################################
###########       Look for proteins found at the intersection    ###########
###########       of two classes annotated to distant GOs        ###########
############################################################################
    case ['i','inter']{
	die("Need an annotated class file for mode 'i'\n") unless $annotations_file;
      bait:foreach my $bait (keys(%proteins)){
	  $c++;
	  my %seen_classes;
	  my $is_cand=0;
	  printf STDERR ("$c of $protnum\r") if $verbose;
	  if ((scalar(keys(%{$proteins{$bait}{CLASS}})) == 1) && (defined($gold{$bait}))){
	      push@{$gold{$bait}{MISSED}{$mode}},"MonoClassed";
	      next;
	  }
	  my @bait_class_annots=();
	  my @cls=keys(%{$proteins{$bait}{CLASS}});
	c1:for(my $n=0; $n<=$#cls;$n++){
	  c2:for(my $k=$n+1; $k<=$#cls;$k++){
	    c3:foreach my $class1_annotation (@{$classes{$cls[$n]}{ANNOT}}){
		if($class1_annotation=~/0008150/){
		    push@{$gold{$bait}{MISSED}{$mode}},"0008150 $cls[$n]" if defined($gold{$bait}); 
		    next c3;
		}
		my $prec1=&term_precision($class1_annotation);
	      c4:foreach my $class2_annotation (@{$classes{$cls[$k]}{ANNOT}}){		  
		  next c4 if $cls[$k] eq $cls[$n];
		  if($class2_annotation=~/0008150/){
		      push@{$gold{$bait}{MISSED}{$mode}},"0008150 $cls[$k]" if defined($gold{$bait});
		      next c4;
		  }	
		  my $prec2=&term_precision($class1_annotation);		  
		  my $P;
		  ## If we are checking BOTH annot and inter probs
		  if($check_both){
		      my ($Pi,$Pa)=(&interactome_probability($class1_annotation,$class2_annotation),&gopair_probability($class1_annotation,$class2_annotation));	
		      if($Pa<=$prob && $Pi<=$prob){
			  my $a=$cls[$n] . "_" . $cls[$k];
			  my $term1=$class1_annotation;
			  my $term2=$class2_annotation;
			  ## Make sure we keep the MOST specific LEAST probable of the interactions
			  ## in question
			  if(defined($moonlighting{$bait}{$a})){
			      $moonlighting{$bait}{$a}=~/^.+?\t.+?GO:.+?\s([\.\d]+).+?\t.+?GO:.+?\s([\.\d]+).+:\s*(.+?)\t(.+?)$/ || die("AA $moonlighting{$bait}{$a}\n");
			      ##prec1,prec2,inter_prob,annot_prob
			      my ($p1,$p2,$pi,$pa)=($1,$2,$3,$4);
			      ## If this is the most precise case, keep
			      if($p1+$p2 > $prec1+$prec2){
				  $moonlighting{$bait}{$a}="$bait\t$cls[$n] $term1 $prec1 " . 
				      scalar(keys(%{$classes{$cls[$n]}{PROTS}})) . 
				      "\t$cls[$k] $term2 $prec2 " . 
				      scalar(keys(%{$classes{$cls[$k]}{PROTS}})) . 
				      "\t: $Pi\t$Pa";
			      }
			      ## If it is AS precise, keep if less likely
			      elsif($p1+$p2 == $prec1+$prec2){
				  $moonlighting{$bait}{$a}="$bait\t$cls[$n] $term1 $prec1 " . 
				      scalar(keys(%{$classes{$cls[$n]}{PROTS}})) . 
				      "\t$cls[$k] $term2 $prec2 " . 
				      scalar(keys(%{$classes{$cls[$k]}{PROTS}})) . 
				      "\t: $Pi\t$Pa" if $pa<=$Pa;   
			      }
			  }
			  else{
			      $moonlighting{$bait}{$a} ="$bait\t$cls[$n] $term1 $prec1 " . 
				  scalar(keys(%{$classes{$cls[$n]}{PROTS}})) . 
				  "\t$cls[$k] $term2 $prec2 " . 
				  scalar(keys(%{$classes{$cls[$k]}{PROTS}})) . 
				  "\t: $Pi\t$Pa";
			  }
#				    next bait;
		      }
		      else{
			  &check_gold($mode,$bait,$class1_annotation,$class2_annotation,$Pi . "_" . $Pa,0);
			  next c1;
		      }
		  }
		  ## If we are NOT checking BOTH inter and annot probs
		  else{
		      $check_interactome_probs ?   
			  ($P=&interactome_probability($class1_annotation,$class2_annotation)) : 
			  ($P=&gopair_probability($class1_annotation,$class2_annotation));
		      if($P<=$prob){
			  my $a=$cls[$n] . "_" . $cls[$k];
			  # my $term1=&terms_to_GOs($class1_annotation,1);
			  # my $term2=&terms_to_GOs($class2_annotation,1);
			  my $term1=$class1_annotation;
			  my $term2=$class2_annotation;
			  ## Beautify output
			  # while(length($term1)<50){
			  #     $term1.=" ";
			  # }
			  # while(length($term2)<50){
			  #     $term2.=" ";
			  # }
			  
			  ## Make sure we keep the MOST specific LEAST probable of the interactions
			  ## in question
			  if(defined($moonlighting{$bait}{$a})){
			      $moonlighting{$bait}{$a}=~/^.+?\t.+?GO:.+?\s([\.\d]+).+?\t.+?GO:.+?\s([\.\d]+).+:\s(.+?)$/ || die("AA $moonlighting{$bait}{$a}\n");
			      my ($p1,$p2,$pp)=($1,$2,$3);
			      ## If this is the most precise case, keep
			      if($p1+$p2 > $prec1+$prec2){
				  $moonlighting{$bait}{$a}="$bait\t$cls[$n] $term1 $prec1 " . 
				      scalar(keys(%{$classes{$cls[$n]}{PROTS}})) . 
				      "\t$cls[$k] $term2 $prec2 " . 
				      scalar(keys(%{$classes{$cls[$k]}{PROTS}})) . 
				      "\t: $P";
			      }
			      ## If it is AS precise, keep if less likely
			      elsif($p1+$p2 == $prec1+$prec2){
				  $moonlighting{$bait}{$a}="$bait\t$cls[$n] $term1 $prec1 " . 
				      scalar(keys(%{$classes{$cls[$n]}{PROTS}})) . 
				      "\t$cls[$k] $term2 $prec2 " . 
				      scalar(keys(%{$classes{$cls[$k]}{PROTS}})) . 
				      "\t: $P" if $pp<=$P;
			      }
			  }
			  else{
			      $moonlighting{$bait}{$a} ="$bait\t$cls[$n] $term1 $prec1 " . 
				  scalar(keys(%{$classes{$cls[$n]}{PROTS}})) . 
				  "\t$cls[$k] $term2 $prec2 " . 
				  scalar(keys(%{$classes{$cls[$k]}{PROTS}})) . 
				  "\t: $P";
			  }
#				    next bait;
		      }
		      else{
			  &check_gold($mode,$bait,$class1_annotation,$class2_annotation,$P,0);
			  next c1;
		      }
		  } ## end if/else check_both
	      }
	    }
	  }
	}
      }
	
    }


    # case ['ii','iinter']{
    # 	die("Need an annotated class file for mode 'i'\n") unless $annotations_file;
    #   bait:foreach my $bait (keys(%proteins)){
    # 	  $c++;
    # 	  my %seen_classes;
    # 	  my $is_cand=0;
    # 	  printf STDERR ("$c of $protnum\r") if $verbose;
    # 	  next bait unless scalar(keys(%{$proteins{$bait}{CLASS}})) > 1;
    # 	  my @bait_class_annots=();
    # 	  my @classes=keys(%{$proteins{$bait}{CLASS}});
    # 	  class1:foreach my $bait_class1 (@classes){
    # 	     my @class1_annots=@{$classes{$bait_class1}{ANNOT}};
    #           ## Add bait's direct annotations
    # 	      my %k;
    # 	      #map{$k{$_}++ if $proteins{$bait}{GOs}{$_} eq 'DIR'}keys(%{$proteins{$bait}{GOs}});
    # 	      #push @class1_annots,keys(%k);
    # 	   class2:foreach my $bait_class2 (@classes){
    # 		next class2 if $bait_class1 eq $bait_class2;
    # 		#my @class2_annots=@{$classes{$bait_class2}{ANNOT}};
    # 		#print STDERR "\$classes{$bait_class2}{ANNOT} : $classes{$bait_class2}{ANNOT} : @{$classes{$bait_class2}{ANNOT}}\n";
    # 		foreach my $class1_annotation (@class1_annots){
    # 		    aa:foreach my $class2_annotation (@{$classes{$bait_class2}{ANNOT}}){
    # 			next class2 if $seen_classes{$class1_annotation}{$class2_annotation};
    # 			my $P=&gopair_probability($class1_annotation,$class2_annotation,"main_classes",$bait,"NOTARGET");
    # 			if($P>$low_prob){
    # 			    next aa;
    # 			}
    # 			else{
    # 			    # $is_cand=1;
    # 			    # if($check_ancestors){
    # 			    # 	if (&check_ancestor_gos($class1_annotation,$class2_annotation,$bait,"NOTARGET") == 1){			   			    
    # 			    # 	    $is_cand=1;
    # 			    # 	}
    # 			    # 	else{$is_cand=0;}
    # 			    # }
			   
    # 			    # $seen_classes{$class1_annotation}{$class2_annotation} =
    # 			    # 	$seen_classes{$class2_annotation}{$class1_annotation} =	1;
    # 			     #if($is_cand==1){
    # 				my $a=$bait_class1 . "-" . $bait_class2;
    # 				my $term1=&terms_to_GOs($class1_annotation,1);
    # 				my $term2=&terms_to_GOs($class2_annotation,1);
				
    # 				## Beautify output
    # 				while(length($term1)<50){
    # 				    $term1.=" ";
    # 				}
    # 				while(length($term2)<50){
    # 				    $term2.=" ";
    # 				}		 
    # 				$moonlighting{$bait}{$a} ="$bait\t($bait_class1,$bait_class2)   \t$term1\t$term2\t$P" if $is_cand==1;
    # 			    #}
    # 			}
    # 		    }
    # 		}
    # 	    }
    # 	  }
    #   } ## end foreach bait
   # } ## end case "i" 
############################################################################
###########       Look for improbable interactions taking        ###########
###########       only multiclassed proteins as candidates       ###########
############################################################################

    case ['m','multi']{
	my $multi_num;
	if($verbose){
	    foreach my $bait (keys(%proteins)){
		my $a=scalar(keys(%{$proteins{$bait}{CLASS}}));
		$multi_num++ if $a>1;
	    }
	}
	die("Need an annotated class file for mode 'm'\n") unless $annotations_file;
	my $c=0;
      bait:foreach my $bait (keys(%proteins)){
	  ## this will hold the smallest asoc prob of the GOs
	  ## annotating this class
	  my $class_P=1; 
	  ###############
	  my @silly_tmp_array=keys(%{$proteins{$bait}{CLASS}});
	  if(defined($gold{$bait})){
	      push@{$gold{$bait}{MISSED}{$mode}},"UnClassed" if $#silly_tmp_array==-1;
	      push@{$gold{$bait}{MISSED}{$mode}},"MonoClassed" if $#silly_tmp_array==0;
	  }
	  next bait unless scalar(keys(%{$proteins{$bait}{CLASS}})) > 1;
	  $c++;	  
	target:foreach my $target (keys(%{$proteins{$bait}{INTERACTORS}})){
	   ## Skip if we have already seen this
	    if( defined($not_moonlighting{$bait}{$target}) ||
		defined($moonlighting{$bait}{$target}) ||
		defined($moonlighting{$target}{$bait}) ||
		defined($not_moonlighting{$target}{$bait})){
		next target
	    }
	    &debug("$bait ::: $target");
	    bgo:foreach my $bgo (keys(%{$proteins{$bait}{GOs}})){
		next bgo if $bgo eq 'GO:0008150';
		tgo:foreach my $tgo (keys(%{$proteins{$target}{GOs}})){
		    next tgo if $tgo eq 'GO:0008150';
		    &debug("\t $bgo ::: $tgo\n");
		    printf STDERR ("$c of $multi_num\r");
		    my $P;
                    ## If we are checking BOTH annot and inter probs
		    if($check_both){
			my ($Pi,$Pa)=(&interactome_probability($bgo,$tgo),&gopair_probability($bgo,$tgo));
			if($Pa<=$prob && $Pi<=$prob){
			    my $term1=&terms_to_GOs($bgo,1) . ",";
			    my $term2=&terms_to_GOs($tgo,1) . ",";
			    while(length($term1)<50){$term1.=" ";}
			    while(length($term2)<50){$term2.=" ";}
			    ## Make sure we keep the MOST specific LEAST probable of the interactions
			    ## in question
			    if(defined($moonlighting{$bait}{$target})){
				$moonlighting{$bait}{$target}=~/.+\t:\s(.+?)\t(.+?)$/;
				my ($ka, $ki)=($1,$2);
				$moonlighting{$bait}{$target}="$bait ($term1 $bgo,$proteins{$bait}{GOs}{$bgo})  \t:\t$target ($term2 $tgo,$proteins{$target}{GOs}{$tgo})  \t: $Pi\t$Pa" if $Pa<$ka;
			    }
			    else{
				$moonlighting{$bait}{$target} ="$bait ($term1 $bgo,$proteins{$bait}{GOs}{$bgo})  \t:\t$target ($term2 $tgo,$proteins{$target}{GOs}{$tgo})  \t: $Pi\t$Pa"; 
			    }
			    next target;
			}
			# Discard target if ANY of its GOs are related
			# to ANY of the bait's.
			else{
			    &check_gold($mode,$bait,$bgo,$tgo,$Pi . "_" . $Pa);
			    $not_moonlighting{$bait}{$target}=1;
			    next target;
			}
		    }		
		    ## If we are NOT checking BOTH inter and annot probs
		    else{
			$check_interactome_probs ?   
			    ($P=&interactome_probability($bgo,$tgo)) : 
			    ($P=&gopair_probability($bgo,$tgo));
			if ($P<=$prob){
			    if($check_ancestors){
				if (&check_ancestor_gos($bgo,$tgo,$bait,$target) == 1){
				    $moonlighting{$bait}{$target} ="$bait ($bgo,$proteins{$bait}{GOs}{$bgo})  \t:\t$target ($tgo,$proteins{$target}{GOs}{$tgo})  \t: $P";
				    ## Skip to next bait since we have identified this one
				    ## as moonlighting
				    next bait;
				}
			    }
			    else{
				my $term1=&terms_to_GOs($bgo,1) . ",";
				my $term2=&terms_to_GOs($tgo,1) . ",";

				while(length($term1)<50){
				    $term1.=" ";
				}
				while(length($term2)<50){
				    $term2.=" ";
				}
				## Make sure we keep the LEAST probable of the interactions
				## in question	 
				if(defined($moonlighting{$bait}{$target})){
				    $moonlighting{$bait}{$target}=~/.+\t:\s(.+?)$/;
				    my $koko=$1;
				    $moonlighting{$bait}{$target}="$bait ($term1 $bgo,$proteins{$bait}{GOs}{$bgo})  \t:\t$target ($term2 $tgo,$proteins{$target}{GOs}{$tgo})  \t: $P" if $P<$koko;
				}
				else{
				    $moonlighting{$bait}{$target} ="$bait ($term1 $bgo,$proteins{$bait}{GOs}{$bgo})  \t:\t$target ($term2 $tgo,$proteins{$target}{GOs}{$tgo})  \t: $P"; 
				}
				next target;
			    }
			}
			# Discard target if ANY of its GOs are related
			# to ANY of the bait's.
			else{
			    &check_gold($mode,$bait,$bgo,$tgo,$P);
			    $not_moonlighting{$bait}{$target}=1;
			    next target;
			}
		    } ## end if/else check_both
		} ## end foreach tgo
	    } ## end foreach bgo
	} ## end foreach target
      } ## end foreach bait
    } ## end case "m"

    
    case 'o'{
       	my $counter=0;
      bait:foreach my $bait (keys(%proteins)){
	  $counter++;
	  print "bb : $bait : ", scalar(keys(%{$proteins{$bait}{CLASS}})), "\n";
	  next bait unless scalar(keys(%{$proteins{$bait}{CLASS}})) > 1;
	  
      }
    } 
############################################################################
########### If the bait interacts with >=50% with a specific GO, ###########
########### identify interactions with another                   ###########
############################################################################
    case ['p','perc'] {
      bait:foreach my $bait (keys(%proteins)){
	  &check_gold($bait);
	  $c++;
	  printf STDERR ("$c of $protnum\r") if $verbose;
	  my (@target_gos,@target_names);
	  my %bait_info;
	  my @bait_gos=keys(%{$proteins{$bait}{GOs}});
	target:foreach my $target (keys(%{$proteins{$bait}{INTERACTORS}})){
	    map{
		$bait_info{$bait}{$target}{GOs}{$_}++;
	    }keys(%{$proteins{$target}{GOs}});
	    push @{$bait_info{$bait}{TARGETS}},$target;
	}## end foreach $target
	  &calculate_percentages(\%bait_info); ## this will populate %moonlighting
	 
      }## end foreach bait
    } #end case "p"

############################################################################
###########       Look for proteins whose interaction partners   ###########
###########            are annotated ONLY to dissimilar GOs      ###########
############################################################################
 
    case ['P','pairs']{ 
      bait:foreach my $bait (keys(%proteins)){
	  $c++;
	  my $is_cand='nope';
	  printf STDERR ("$c of $protnum\r") if $verbose;
	  my @targets=keys(%{$proteins{$bait}{INTERACTORS}});
	t1:foreach  my $target1(keys(%{$proteins{$bait}{INTERACTORS}})){
	    next t1 if $bait eq $target1;
	  t2:foreach  my $target2(keys(%{$proteins{$bait}{INTERACTORS}})){
	      next t2 if $target1 eq $target2;
	      next t2 if $bait eq $target2;
	    go1:foreach my $go1 (keys(%{$proteins{$target1}{GOs}})){
		next go1 if $go1 eq 'GO:0008150';
	      go2:foreach my $go2 (keys(%{$proteins{$target2}{GOs}})){
		  next go2 if $go2 eq 'GO:0008150';
		  ## If we are checking BOTH annot and inter probs
		  if($check_both){	  
		      my ($Pi,$Pa)=(&interactome_probability($go1,$go2),&gopair_probability($go1,$go2));	
		      if($Pa<=$prob && $Pi<=$prob){
			  my ($last_pa,$last_pi)=(10,10);
			  if ($is_cand=~/\t:\s*(.+?)\t(.+?)$/){
			      $last_pi=$1;
			      $last_pa=$2;
			  }
			  if($Pa<$last_pa and $Pi<$last_pi){
			      $is_cand="($go1)\t$target2 ($go2)\t: $Pi\t$Pa";
			  }
		      }
		      else{
			  &check_gold($mode,$bait,$target1 . "-" . $target2,$go1 . "-" . $go2,$Pi . "_" . $Pa);
			  $is_cand='nope';
			  next t2;
		      }
		  }
		  ## If we are NOT checking BOTH inter and annot probs
		  else{
		      my $P;
		      my $target=$target1 . "_" . $target2; ### for debugging purposes only
		      $check_interactome_probs ?   
			  ($P=&interactome_probability($go1,$go2,"main",$bait,$target)) : 
			  ($P=&gopair_probability($go1,$go2,"main",$bait,$target));
		      if ($P<=$prob){
			  my $last_p=10;
			  if ($is_cand=~/\t:\s*(\d.+)$/){$last_p=$1;}
			  $is_cand="($go1)\t$target2 ($go2)\t: $P" if $P<$last_p;
		      }
		      else{
			  &check_gold($mode,$bait,$target1 . "-" . $target2,$go1 . "-" . $go2,$P);
			  $is_cand='nope';
			  next t2;
		      }
		  }
	      }
	    }
	  }
	    
	    if($is_cand !~ /^nope$/){
		$moonlighting{$bait}="$bait\t$target1 $is_cand";
	    }
	    else{}
	}
      } ## end foreach bait
    } ## end case 'P'
} ## end switch;


########################### Print Results ########################
print STDERR "\n";

my @oo=keys(%moonlighting);
my $pp=0;
my $color="bold blue";
my %output;
my %gfound;
foreach my $bait (keys(%moonlighting)){
    my $is_gold=$gold{$bait}||undef;
    $only_gold && do {next unless $is_gold; };
    my $P;					      
    if(  ($mode eq 'pairs') || ($mode eq 'P')){
	next if defined($not_moonlighting{$bait});
	$moonlighting{$bait}=~/.+:\s*(.+?)$/||die("Could not match mprob : $moonlighting{$bait}\n");
	$P=$1;
	if($check_both){
	    ##l[0] is the inter prob and l[1] is the annot prob
	    my @l=split(/\t/,$P);
	    my $Pi=$l[0]*$TESTS;
	    my $Pa=$l[1]*$TESTS;
	    next unless $Pa<=$prob;
	    next unless $Pi<=$prob;
	    $moonlighting{$bait}=~s/$P$/$Pi\t$Pa/;
	    if($is_gold){ 
		$moonlighting{$bait}=~s/$bait/color("$color").$bait.color("reset")/ige if $opts{c};
		$gfound{$bait}++;
	    }
	    push @{$output{$Pa}},$moonlighting{$bait};
	}
	else{
	    my $P_corr=$P*$TESTS;
	    next unless $P_corr<=$prob;
	    $moonlighting{$bait}=~s/$P$/$P_corr/;
	    if($is_gold){ 
		$moonlighting{$bait}=~s/$bait/color("$color").$bait.color("reset")/ige if $opts{c};
		$gfound{$bait}++;
	    }
	    push @{$output{$P}},$moonlighting{$bait};
	}
    }
    ## For all other modes
    else{
	foreach my $target (keys(%{$moonlighting{$bait}})){
	    next if defined($not_moonlighting{$bait}{$target});
	    $moonlighting{$bait}{$target}=~/.+\t:\s*(.+?\t.+?)$/||die("Could not match mprob $bait:$target: $moonlighting{$bait}{$target}\n");
	    $P=$1;
	    ## If we are checking BOTH annot and inter probs
	    if($check_both){
		##l[0] is the inter prob and l[1] is the annot prob
		my @l=split(/\t/,$P);
		my $Pi=$l[0]*$TESTS;
		my $Pa=$l[1]*$TESTS;
		next unless $Pa<=$prob;
		next unless $Pi<=$prob;
		$moonlighting{$bait}{$target}=~s/$P$/$Pi\t$Pa/;
		if($is_gold){ 
		    $moonlighting{$bait}{$target}=~s/$bait/color("$color").$bait.color("reset")/ige if $opts{c};
		    $gfound{$bait}++;
		}
		push @{$output{$Pa}},$moonlighting{$bait}{$target};
	    }
	    ## If we are NOT checking both inter and annot probs
	    else{
		my $P_corr=$P*$TESTS;
		next unless $P_corr<=$prob;
		$moonlighting{$bait}{$target}=~s/$P$/$P_corr/;
		if($is_gold){ 
		    $moonlighting{$bait}{$target}=~s/$bait/color("$color").$bait.color("reset")/ige if $opts{c};
		    $gfound{$bait}++;
		}
		push @{$output{$P}},$moonlighting{$bait}{$target};
	    }
	}
    }
    $pp++;
}
my @sorted=sort{$a <=> $b} keys(%output);
map{
    map{
	print "$_\n";
    }@{$output{$_}};
}@sorted;

## Print missed gold prots
my $gout;
$check_interactome_probs ? 
    ($gout =$network_file . ".gold." . $species . "." . $mode . ".inter.log") :
    ($gout =$network_file . ".gold." . $species . "." . $mode . ".log");
open(G,">$gout")||die("Could not open $gout : $!\n");
print G "\n==== $date ======\n";
foreach my $g_prot (keys(%gold)){
    my $cc=0; ## just for printing a newline to beautify output
#    $mode eq 'b' && do {next if defined($moonlighting{$g_prot})};
    next if defined($gfound{$g_prot});
    if(($g_prot =~/MISSED/)||($g_prot =~/READ/)){next}
    if ($gold{$g_prot}{NUM}==0){
	print G "$g_prot\tAbsent\n";
	next;
    }
    print G "$g_prot\t";
    map{print G "$g_prot\t" if $cc==1; print G "$_\n"; $cc=1}@{$gold{$g_prot}{MISSED}{$mode}};
    print G "0008150\n" unless $cc==1;
}
my @ff=keys(%gfound);
my $cmnd=`ps aux | grep $$ | grep -v grep`;
$cmnd=~s/^.+(moonG.+?)$/$1/;
print G "\n$cmnd : " . scalar(@ff) . " gold standard proteins were found: @ff\n";
close(G);
#################### Create network specific probability file ####
if($read_network_prob_file){
    unless($added_to_prob_file==0){
	print STDERR "Generating prob file...\n" if $verbose;
	open(NET,">>$network_prob_file");
	map{print NET}@net_probs;
	system("cat $network_prob_file | sort | uniq > lililolo.$$; mv lililolo.$$ $network_prob_file");
	close(NET);
    }
}


################################ END ##############################

############################### SUBROUTINES ################################
sub parse_input_files_guess_values{
    ## Change long mode names to short
    if($mode eq 'dis2class'){$mode = 'a'}
    elsif($mode eq 'bridge'){$mode = 'b'}
    elsif($mode eq 'difclass'){$mode = 'c'}
    elsif($mode eq 'inter'){$mode = 'i'}
    elsif($mode eq 'multi'){$mode = 'm'}
    elsif($mode eq 'perc'){$mode = 'p'}
    elsif($mode eq 'pairs'){$mode = 'P'}
    elsif($mode eq 'direct'){$mode = 'd'}

    ## Guess species, assign variables accordingly
    if(($network_file=~/^dro/) || ($network_file =~/^dm/i)){$species="fly"}
    elsif(($network_file=~/^hs/) || ($network_file=~/^hs/i)){$species="human"}
    elsif(($network_file=~/^mus/) || ($network_file=~/^mm/i)){$species="mouse"}
    elsif(($network_file=~/^scc/) || ($network_file=~/^sc/i)){$species="yeast"}
    elsif(($network_file=~/^ele/) || ($network_file=~/^ce/i)){$species="worm"}

    ## Check ontology
    if($subonto eq "p"){$subonto="P"}
    elsif($subonto eq "c"){$subonto="c"}
    elsif($subonto eq "f"){$subonto="f"}
    elsif($subonto eq "F" || "P" || "C"){}
    else{print STDERR "Subontology (-o) must be one of \"p,c,f\"\n"; exit(1)}
    
    ## Get appropriate files accoriding to species
    if   (($species eq "human") || ($species eq "hum") || ($species eq "hs")){
	$species='human';
	$stats_dir="$DATADIR/gostats/human";
	$synfile=$opts{s}||"$DATADIR/hs.uni2acc.map";
	$prob_file=$stats_dir . "/human.prob";
	$gaf_annotations_file=$opts{f} || "$DATADIR/gene_association.human";
	$suffix='HUMAN';
	$low_annot_prob= 0.9393275 unless $opts{p};
	$low_inter_prob= 0.9157528 unless $opts{p};    
    }
    elsif(($species eq "fly")   || ($species eq "fb")  || ($species eq "dm") || ($species eq "dro")){
	$species='fly';
	$stats_dir="$DATADIR/gostats/fly";
	$synfile=$opts{s}||"$DATADIR/dro.uni2fb.map";
	$prob_file=$stats_dir . "/fly.prob";
	$gaf_annotations_file=$opts{f} || "$DATADIR/gene_association.fly";
	$suffix='DROME';
 	$low_annot_prob= 0.8450468 unless $opts{p};
	$low_inter_prob= 0.9920737 unless $opts{p};    
   }
    elsif(($species eq "worm")  || ($species eq "wb")  || ($species eq "ce") || ($species eq "ele")){
	$species='worm';
	$stats_dir="$DATADIR/gostats/worm";
	$synfile=$opts{s}||"$DATADIR/ele.uni2wb.map";
	$prob_file=$stats_dir . "/worm.prob";
	$gaf_annotations_file=$opts{f} || "$DATADIR/gene_association.worm";
	$suffix='CAEEL';
	$low_annot_prob= 0.8819694 unless $opts{p};
	$low_inter_prob= 0.9900118 unless $opts{p};    

    }
    elsif(($species eq "mouse") || ($species eq "mg")  || ($species eq "mm") || ($species eq "mus")){
	$species='mouse';
	$stats_dir="$DATADIR/gostats/mouse";
	$prob_file=$stats_dir . "/mouse.prob";
	$synfile=$opts{s}||"$DATADIR/mus.uni2mgi.map";
	$gaf_annotations_file=$opts{f} || "$DATADIR/gene_association.mouse";
	$suffix='MOUSE';
	$low_annot_prob= 0.9540235 unless $opts{p};
	$low_inter_prob= 0.9948855 unless $opts{p};    

    }
    elsif(($species eq "yeast") || ($species eq "sgd") || ($species eq "scc")){
	$species='yeast';
	$stats_dir="$DATADIR/gostats/yeast";
	$prob_file=$stats_dir . "/yeast.prob";
	$synfile=$opts{s}||"$DATADIR/scc.uni2sgd.map";
	$gaf_annotations_file=$opts{f} || "$DATADIR/gene_association.yeast";
	$suffix='YEAST';
	$low_annot_prob= 0.7946262 unless $opts{p};
	$low_inter_prob= 0.851656 unless $opts{p};    


    }
    else{print STDERR   <<EndOfMsg;
	 Species must be one of the following:
	 "human" || "hum" || "hs"
	     "fly"   || "fb"  || "dm" || "dro"
	     "worm"  || "wb"  || "ce" || "ele"
	     "mouse" || "mg"  || "mm" || "mus"
	     "yeast" || "sgd" || "scc"
EndOfMsg
	     exit(0);
    }
    ## Choose which threshold to use
    if($check_interactome_probs){
	$prob=$low_inter_prob;
    }
    else{
	$prob=$low_annot_prob;
    }

    
}
###############################################################
sub load_network {
    print STDERR "Loading Network...\n" if $verbose;
    my $A = IO::File->new("< $network_file")|| die("Cannot open $network_file : $!\n");
    while(<$A>){
	next if /^\d+$/;
	next if /^!/;
	next if /^\s*$/;
	my ($bait,$target)=split(/\s+/,$_);
	&debug("$. : $bait $target");
	$bait=&get_name($bait,"load_network");
	$target=&get_name($target,"load_network");
	## %proteins will hold all the info on all the 
	## network's proteins
	#&debug("PROT : \$proteins{$bait}{INTERACTORS}{$target}\n    \$proteins{$target}{INTERACTORS}{$bait}");
	$proteins{$bait}{INTERACTORS}{$target}++;
	$proteins{$target}{INTERACTORS}{$bait}++;
    }
    close(A);
}

############################################################
sub get_name{
    my $name=$_[0];
    my $old_name=$name;
    my $called_by=$_[1]||"null";
    my $reverse=$_[2]||undef;
    if ($synonyms{LOADED}==0){
	if(-e $synfile){
	    open(S,$synfile);
	    while(<S>){
		## If this is a type 1 synonyms file, ie, 
		## network_name\tGAF_name\tAccession
		if(/^(.+)\t(.+)\t(.+)\n/){
		    $synonyms{ACC}{$3}=$3;
		    $synonyms{ACC}{$2}=$3;
		    $synonyms{ACC}{$1}=$3;
		    $synonyms{NAME}{$1}=$2;
		    $synonyms{NAME}{$2}=$2;
		    $synonyms{NAME}{$3}=$2;
		}	
                ## If this is a type 2 synonyms file, ie, 
		## gene name\tlist of synonyms
		elsif(/^(.+)\t(.+\s.+)\n/){
		       my $nn=$1;
		       ## Discard protein if there are naming problems
		       my @syns=split(/\s/,$2);
		       $synonyms{GOOD}{$nn}++;
		       for (my $n=0;$n<scalar(@syns); $n++){
			   $synonyms{$syns[$n]}=$nn;
		       }
		   }
		#if this is a GAF specific synonyms file, i.e.
		# network_name\tGAF_name
		elsif(/^([^\t]+?)\t([^\t]+)\n/){
		    $synonyms{NAME}{$2}=$1; 
		    $synonyms{ACC}{$2}=$2; 
		    $synonyms{NAME}{$1}=$1; 
		    $synonyms{ACC}{$1}=$2; 
		   # print "\$synonyms{NAME}{$2} : $synonyms{NAME}{$2}\n\$synonyms{ACC}{$2} : $synonyms{ACC}{$2}";
		}
		else{
		    die("Bad format synonyms file : $synfile\n");
		}
	    }
	    close(S);
	}
	else{die("Could not open $synfile\n" );
	    # print STDERR "No synonyms file found, parsing annotations... \n NEED TO MODIFY THIS BECAUSE OF THE CHANGE IN SYNONYMS FORMAT" if $verbose;
	#     open(A,"$ARGV[0]")|| die("Cannot open synonyms $ARGV[0]:$!\n");
	#     while(<A>){
	# 	my $silivar="!";
	# 	next if /^$silivar/; ## silly thing cause emacs screws up the syntax highlihting when /!/;
	# 	my @tmparray=split(/\t/);
	# 	map{$synonyms{$_}=$tmparray[2]}split(/\|/,$tmparray[10]);
	#     }
	#     close(A);
	 }	
 	$synonyms{LOADED}=1;
    }
    my $hname;
    $name =~ /_$suffix/i ? ($hname=$name) : ($hname=$name . "_" . $suffix);
    # if(defined($synonyms{ACC}{$name})){
    # 	$name=$synonyms{ACC}{$name};
    # }
    # elsif(defined($synonyms{ACC}{$hname})){
    # 	$name=$synonyms{ACC}{$hname};
    # }
    # elsif(defined($synonyms{GOOD}{$hname})){die("111 : $name\n");
    # 	$name=$synonyms{GOOD}{$hname};
    # }
    if(defined($synonyms{NAME}{$name})){
	$name=$synonyms{NAME}{$name} ;
    }
    elsif(defined($synonyms{NAME}{$hname})){
	$name=$synonyms{NAME}{$hname};
    }
    elsif ($name=~/(.+)_$species/i && defined($synonyms{$1})){die("333 : $name\n");
	$name=$synonyms{$1};
    }
    #else{die("Shit, name problem: $name : $synonyms{NAME}{$name}  $synonyms{ACC}{$name}\n")}
	    
    if ($very_verbose){&debug("Called by : $called_by for $old_name ($synfile), returned $name"); }
#    print "xxx $old_name:$name\n";
    
    return $name;
   
}
############################################################
sub parse_gaf_file{    
    print STDERR "Parsing GAF file...\n" if $verbose;
    open(A,"$gaf_annotations_file")|| die("Cannot open $gaf_annotations_file:$!\n");
    while(<A>){
	next if /^!/;
	chomp;
	my @tmparray=split(/\t/,$_);
	next unless $tmparray[8] eq $subonto;
	next if $tmparray[3] !~/^$/;
	my $name=&get_name($tmparray[1],"parse_gaf_file");
	## If this protein is in the graph
	if(defined($proteins{$name})){
	    
	    $proteins{$name}{GOs}{$tmparray[4]}="DIR";
	    ## Inherit ancestor annotations
	    map{$proteins{$name}{GOs}{$_}="ANC"}@{$ancestors{$tmparray[4]}} unless $no_anc;
	}
    }
    close(A);
}
############################################################
sub load_ancestors{
    open (ANC,"$geneology_file")|| die("cannot open $geneology_file:$!\n");
    while(<ANC>){
	chomp;
	my @terms=split(/\t/);
	my %seen;
	@terms = grep { ! $seen{$_} ++ } @terms;
	my $child=shift(@terms);
	$ancestors{$child}=[@terms];
	map{$offspring{$_}{$child}++}@terms;
    }
    close(ANC);
}
############################################################

sub load_annotations{
    print STDERR "Loading annotations..." if $verbose;
    my ($class,$prot);
    my $file_type=0;
    open(A,"$annotations_file")|| die("Cannot open $_[0]:$!\n");
    my $k=0;
    my (@GOs, @bGOs);
    while(<A>){
	if($.==1){
	    if(/^\[CLASS/){$file_type=1}
	}
	## If this is an annotated classes file
	if($file_type==1){
#	    print if /^\[/;
	    if(/^\[CLASS:\s*(\d+)/){$class=$1}
	    if(/^CA/){
		@GOs=();
		@bGOs=();
		/^CA\s+(.+?)$/;
		my $kk=$1;
		my @ko=split(/\s+/,$kk);
		foreach my $term(@ko){
		    my $T=&terms_to_GOs($term,0);
		    push @GOs,$T unless $T eq 'BAD';
		    ## Inherit ancestor annotations
		    map{push @GOs, $_ unless $_ eq 'BAD';}@{$ancestors{$T}} unless $no_anc;
		    # $foundGOs{$T}++;
		}
	    }
	    if(/^CM/){
	    	/^CM\s+(.+?)\s+(\d+)\/(\d+)/||die("Could not match CM($.): $annotations_file $_\n");
	    	my @ko;
		my @kkk=split(/ /,$1);
	    	my $a=$2*100/$3;
	    	    if ($a>=30){
	    		push @ko, @kkk;
	    	    }
	    	my $kk=$1;
	    	foreach my $term(@ko){
	    	    my $T=&terms_to_GOs($term,0);
	    	    #print STDOUT "pushed $T for $class\n";
	    	    push @bGOs,$T unless $T eq 'BAD';
	    	    # $foundGOs{$T}++;
	    	}
	    }
	    if(/^PN\s*(.+)/){
		my $kk=$1;
		my @a=split(/,\s+/,$kk);
		## get the right names
		my @prots;
		map{push @prots, &get_name($_,'load_annotations');}@a;
		my %seen=();
		foreach my $protein(@prots){
		    next unless defined($proteins{$protein});
		    for (my $n=0; $n<scalar(@GOs); $n++){
			$proteins{$protein}{CLASS}{$class}++;
			$classes{$class}{PROTS}{$protein}++;
			## Add GOs inferred from classes to the prots
			## annotations unless already present
			unless(exists($proteins{$protein}{GOs}{$GOs[$n]})){ 
			    $proteins{$protein}{GOs}{$GOs[$n]}="INF";
			}
			push @{$classes{$class}{ANNOT}}, $GOs[$n] unless $seen{$class}{$GOs[$n]};
			$seen{$class}{$GOs[$n]}++;
		    }
		    ## Also add 2ndary GO annotations
		    for (my $n=0; $n<scalar(@bGOs); $n++){
			unless(exists($proteins{$protein}{GOs}{$bGOs[$n]})){ 
			    $proteins{$protein}{GOs}{$bGOs[$n]}="UNC";
			    #print "$protein $bGOs[$n]\n@bGOs\n"; die();
			}
			push @{$classes{$class}{secondary_ANNOT}}, $bGOs[$n] unless $seen{$class}{$bGOs[$n]};
			$seen{$class}{$bGOs[$n]}++;
		    }
		}		
	    }		
	    
	} # end if($file_type==1)
	else{
	    $k=1 if /^Final\sclasses/i;
	    next if /^Final\sclasses/i;
	    next unless $k>0;
	    next if /^\s*\d+$/;
	    my @a=split(/\s+/);
	    my @b;
	    map{
		my $n=&get_name($_);
		push @b,$n;
	    }@a;
	    $classes{$k}=\@b;
	    
	    ## $proteins{ENPL_HUMAN}{CLASS}{12}=1
	    ## meaning protein ENPL_HUMAN belongs to class 12
	    map{
		if(defined($proteins{$_})){
		    $proteins{$_}{CLASS}{$k}++;
		}
	    }@b;
	    $k++;
	}
    }
    print STDERR "Done\n" if $verbose;

}
############################################################
sub calculate_percentages{
    my %hash=%{(shift)};
    my ($mc,$mcgo,$perc)=0;
    my %gostats;
    my @keys=keys(%hash);
    my $bait=$keys[0];
    %hash=%{$hash{$bait}};    ## hash is now all the info for this bait
    
    ## cycle through ALL target's GOs and count
    foreach my $target (@{$hash{TARGETS}})
    {
    	foreach my $tgo (keys(%{$hash{$target}{GOs}})){ ## tgo is the target's go
    	    $gostats{$tgo}++;
    	}
    }
    foreach my $g (keys(%gostats)){
    	if($gostats{$g} >$mc){
    	    $mc=$gostats{$g}; ##mc number of times mcgo is found
    	    $mcgo=$g; ## mcgo= most common go 
    	}
    }
    $perc=($mc/scalar(@{$hash{TARGETS}}   ))*100;
    ## Skip bait if it does NOT have an mcgo
    if( $perc>=50){
    	# Cycle through target gos, only keep bait
    	# if NO bait gos are related to ANY 
    	# target gos
    	target:foreach my $target (@{$hash{TARGETS}}){
	    my $P;
	    ## skip target if it is annotated to $mcgo
	    next target if defined($hash{$target}{GOs}{$mcgo});
    	    foreach my $tgo (keys(%{$hash{$target}{GOs}})){
		next if $tgo eq 'GO:0008150';
       ## skip target if any of its gos are related to the mcgo
		## If we are checking BOTH annot and inter probs
		if($check_both){
		    my ($Pi,$Pa)=(&interactome_probability($mcgo,$tgo),&gopair_probability($mcgo,$tgo));
		    if($Pa>=$prob && $Pi>=$prob){
			$not_moonlighting{$bait}{$target}=1;
			next target;
		    }
		}
		else{
		    $check_interactome_probs ?   
			($P=&interactome_probability($mcgo,$tgo)) : 
			($P=&gopair_probability($mcgo,$tgo,"gopair_prob"));
		    if($P>=$prob){
			$not_moonlighting{$bait}{$target}=1;
			next target;
		    }
		}
		foreach my $bgo (@{$hash{GOs}}){
		    next if $bgo eq 'GO:0008150';
		    ## If we are checking BOTH annot and inter probs
		    if($check_both){
			my ($Pi,$Pa)=(&interactome_probability($mcgo,$tgo),&gopair_probability($mcgo,$tgo));
			if($Pa>=$prob && $Pi>=$prob){
			    ## skip target if ANY of its gos 
			    ## are close to any bgo
			    $not_moonlighting{$bait}{$target}=1;
			    next target;
			}
			else{
			    ## Check that none of the ancestor GOs
			    ## are related to mcgo either
			    if(&check_ancestor_gos($bgo,$tgo,$bait,$target) == 1){
				## If all is well, keep
				$moonlighting{$bait}{$target}="$bait ($bgo)  \t:\t$target ($tgo)  \t: $Pi\t$Pi";
			    }
			    else{$not_moonlighting{$bait}{$target}=1;}
			}
			
		    }
		    else{
			$check_interactome_probs ?   
			    ($P=&interactome_probability($bgo,$tgo,"gopair_prob")) : 
			    ($P=&gopair_probability($bgo,$tgo,"gopair_prob"));
			if($P>=$prob){
			    ## skip target if ANY of its gos 
			    ## are close to any bgo
			    $not_moonlighting{$bait}{$target}=1;
			    next target;
			}
			else{
			    ## Check that none of the ancestor GOs
			    ## are related to mcgo either
			    if(&check_ancestor_gos($bgo,$tgo,$bait,$target) == 1){
				## If all is well, keep
				$moonlighting{$bait}{$target}="$bait ($bgo)  \t:\t$target ($tgo)  \t: $P";
			    }
			    else{$not_moonlighting{$bait}{$target}=1;}
			}
		    } 
		} ## end foreach my bgo
	    } ## end foreach my tgo
	} #end foreach target
    } ## end if perc >=50
}
############################################################
sub interactome_probability{
    my $go1=shift;
    my $go2=shift;
    my @b = sort {$b lt $a} ($go1,$go2);
    my $gopair=join("_",@b);
    return($inter_prob{$gopair})  if defined($inter_prob{$gopair});
    if (($Iprob_file_read==1) && (not defined($inter_prob{$gopair}))){
	$inter_prob{$gopair}=1;
	return($inter_prob{$gopair});
    }
    if ($go1 eq $go2){  
	$inter_prob{$gopair}=1;
	return($inter_prob{$gopair});
    }
    ## If one term is an ancestor of the other prob=1;
    if((defined($offspring{$go1}{$go2}))||
       (defined($offspring{$go2}{$go1}))){
	$inter_prob{$gopair}=1; 
	return($inter_prob{$gopair});	
    }
    if($read_network_prob_file){
	$network_prob_file = $network_file . ".inter.prob";
	if ((-e "$network_prob_file") && ($Iprob_file_read==0)){
	    open(A,"$network_prob_file")||die("Cannot open $network_prob_file : $!\n");
	    print STDERR "\nReading $network_prob_file..." if $verbose;
	    while(<A>){
		chomp;
		my @tt=split(/\t/);
		my $tempgo=$tt[0];
		$inter_prob{$tempgo}=$tt[1];
	    }
	    close(A);
	    $Iprob_file_read=1;
	    print STDERR "Done\n" if $verbose;
	}
    }
    unless(defined($inter_prob{$gopair})){
	$gopair=~/^GO:(.....)/||die("cannot match gopair : $gopair\n");
	my $file=$1 . ".prob";
	open(A1,"$stats_dir/interactome/$file")|| die("cannot open $stats_dir/interactome/$file : $!\n$gopair\n");
	while(<A1>){
	    if(/$gopair/){
		$added_to_prob_file=1 if $read_network_prob_file;
		chomp;
		my @tt=split(/\t/);
		$inter_prob{$tt[0]}=$tt[6];
		push @net_probs, "$gopair\t" . $inter_prob{$gopair} . "\n" if $read_network_prob_file;
	    }
	}
    }
    defined($inter_prob{$gopair}) ? 
	return($inter_prob{$gopair}):
	return(1);
}
############################################################
sub gopair_probability{
    my $go1=shift;
    my $go2=shift;
    $TESTS++;
    my $source=shift;
    my $bait=shift;
    my $target=shift;
    my $low;
    my @b = sort {$b lt $a} ($go1,$go2);
    my @c = sort {$a lt $b} ($go1,$go2);
    my $bpair=join("_",@c);
    my $gopair=join("_",@b);
    die("$go1,$go2,$source,$bait,$target\n") if $gopair =~/0008150/;
    if (defined($prob{$gopair})){
	return($prob{$gopair}); 
    }
#    print "$gopair\n";
    if ($go1 eq $go2){  
	$prob{$gopair}=1;
	return($prob{$gopair});
    }
    die("aa $gopair: $bpair\n" ) if defined($prob{$bpair});
    ## If one term is an ancestor of the other prob=1;
    if((defined($offspring{$go1}{$go2}))||
       (defined($offspring{$go2}{$go1}))){
	$prob{$gopair}=1;
	return($prob{$gopair});	
    }
    if($read_network_prob_file){
	$network_prob_file = $network_file . ".prob";
	if ( (-e "$network_prob_file") && ($Aprob_file_read==0) ){
	    open(A,"$network_prob_file");
	    print STDERR "\nReading $network_prob_file..." if $verbose;
	    while(<A>){
		chomp;
		my @tt=split(/\t/);
		my $tempgo=$tt[0];
		$prob{$tempgo}=$tt[1];
	    }
	    close(A);
	    $Aprob_file_read=1;
	    print STDERR "Done\n" if $verbose;
	}
    }  

    ## if this prob was not in net prob file
    unless(defined($prob{$gopair})){
	$gopair=~/^GO:(.....)/||die("cannot match gopair : $gopair\n");
	my $file=$1 . ".prob";
	## I have calculated the probabilities of ALL possible
	## GO pairs and split them into files according to the GO
	## number. So, GO:0019952 will be in the file 00199
	if(-e "$stats_dir/annotations/$file.gz") {
	    open(A,"zcat $stats_dir/annotations/$file |")|| die("cannot open $stats_dir/annotations/$file : $!\n$gopair\n");
	}
	elsif(-e "$stats_dir/annotations/$file") {
	    open(A,"$stats_dir/annotations/$file")|| die("cannot open $stats_dir/annotations/$file : $!\n$gopair\n");
	}
#	    open(A,"data/gostats/count_under.HS.BP")|| die("cannot open $stats_dir/$file : $!\n$a1,$a2:$b1,$b2,$go1,$go2,$gopair\n");
	#else{die("bad stats filename: $stats_dir/annotations/$file:$1 $gopair, $bait, $target\n")}
	while(<A>){
	    #   printf STDERR "$. of 37474792\r"  if $verbose;
	    chomp;
	    if(/$gopair/){
		$added_to_prob_file=1 if $read_network_prob_file;
		my @tt=split(/\t/);
		$prob{$gopair}=$tt[6];
		push @net_probs, "$gopair\t" . $prob{$gopair} . "\n" if $read_network_prob_file;
		die("bb  $stats_dir/annotations/$file $gopair: $bpair\n" ) if defined($prob{$bpair});
		return($prob{$gopair});
	    }
	    else{next}
	}
	$have_already_read_probs_file=1;
	close(A);
	
    } ## end unless(defined($prob{$gopair})) 2	
    # } ## end	if($have_already_read_probs_file==0)

    unless(defined($prob{$gopair})){
	## Since I am only giving it the under represented ones,
	## any probs missing will be over represented so we can 
	## safely ignore them.
	# open(SKIPPED, ">>skipped");
	# print SKIPPED "$gopair\n";
	# close(SKIPPED);
	$prob{$gopair}=1;
	push @net_probs, "$gopair\t$prob{$gopair}\n" if $read_network_prob_file;  
	die("cc $gopair: $bpair\n" ) if defined($prob{$bpair});

	return($prob{$gopair});
    }
	## @net_probs will be used to create a prob file for this network
	push @net_probs, "$gopair\t" . $prob{$gopair} . "\n";   
    die("dd $gopair: $bpair\n" ) if defined($prob{$bpair});

	return($prob{$gopair});
}
############################################################

sub check_ancestor_gos{
    my $go1=shift;
    my $go2=shift;
    my $bait=shift;
    my $target=shift;
    &debug("$bait : $target :  $go1: $go2");
    defined($seen{$go1}{$go2}) && do {return($seen{$go1}{$go2})};
    my $good;
    my @go1s=@{$ancestors{$go1}};
    my @go2s=@{$ancestors{$go2}};
    push @go1s,$go1;
    push @go2s,$go2;
    for(my $n=0;$n<=$#go1s; $n++){
	next if $go1s[$n] eq 'GO:0008150';
	for (my $k=0;$k<=$#go2s; $k++){
	    next if $go2s[$k] eq 'GO:0008150';
	    if(&gopair_probability($go1s[$n],$go2s[$k],"anc")>$low_annot_prob){
		$seen{$go1}{$go2}=$seen{$go2}{$go1}=0;
		#print $LOG "Skipped $bait,$target ($go1s[$n],$go2s[$k])\n";
		return(0);
	    }
	}
    }
    $seen{$go1}{$go2}=$seen{$go2}{$go1}=1;
    return(1);
}

############################################################
sub debug{
    if ($debug)
    {
	print STDERR "@_\n";
    }
}
############################################################
sub terms_to_GOs{
    my $term=shift;
    my $mode=shift; ## 0 will return GO:xxx, 1 will return term name
#    unless($term eq 'biological_process'){die("crapiola : -$term-\n");}
    $term=~s/_/ /g;
    if($have_already_read_terms_file==0){
	open(T,"$go_terms_file")|| die("Cannot open terms file : $!\n");
	while(my $line=<T>){
	    next if $line=~/^\!/; 
	    chomp($line);
            ## For some reason, latest GO terms file had
	    ## biological_process instead of biological process.
	    ## So, deal with that:
	    $line=~s/_/ /g; 

	    my @aa=split(/\t/, $line);
	    my @a=@aa;
	    my @terms=($aa[0]);
	    if($line=~/obs$/){pop(@a);}
	    if($aa[1] =~/GO:/){
		push @terms,split(/\s+/,$aa[1]);
	    }
	    if($a[$#a] ne $subonto){
		$terms_to_GOs{TERMS}{$a[$#a-1]}="BAD";
		map{$terms_to_GOs{GOs}{$_}="BAD";}@terms;
	    }
	    else{
		$terms_to_GOs{TERMS}{$a[$#a-1]}=$a[0];
		map{$terms_to_GOs{GOs}{$_}=$a[$#a-1];}@terms;
	    }
	}
	close(T);
	$have_already_read_terms_file=1;
    }

#    &debug("term : $term, id:$terms_to_GOs{TERMS}{$term}, id:$terms_to_GOs{TERMS}{$term} " );
#    print STDERR "term : $term, id:$terms_to_GOs{TERMS}{$term} \n";
    $mode==0 ? 
	return($terms_to_GOs{TERMS}{$term}) :
	return($terms_to_GOs{GOs}{$term}) ;
    
}
############################################################
sub term_precision{
  my $go=shift;
    if($have_already_read_precision_file==0){
	# open(PP,">/tmp/$$")|| die("Could not open precision file /tmp/$$ for writing:$!\n");
	# map{print PP "$_\n";}keys(%found_gos);
	# close(PP);
	# system("term_precision.pl ~/research/GO/biological_process.geneology  /tmp/$$ > /tmp/$$.b");
	open(PO,"$precision_file")|| die "Could not open $precision_file:$!\n";
	while(<PO>){
	    chomp;
	    /(.+)\t(.+)/;
	    $precision{$1}=$2;
	}
	$have_already_read_precision_file=1;
    }
  if(defined($precision{$go})){
      return(nearest(.0001,$precision{$go}));
  }
  else{
      print "MISSING $go\n";
      my $a=`term_precision.pl -n $go $geneology_file | gawk '{print \$NF}'`;
      $precision{$go}=$a;
      return(nearest(.0001,$precision{$go}));

  }
}
############################################################

sub check_gold{
    ## Load gold standard
    unless($gold{READ}){
	open(G,"$gold_file")||die("Could not open $gold_file : $!");
	while(<G>){
	    next if $.==1;
	    my @a=split(/\t/);
	    $gold{$a[0]}{ISO}=$a[1] if defined($proteins{$a[0]});
	    $gold{$a[0]}{GOOD}=$a[2]if defined($proteins{$a[0]});
	    $gold{$a[0]}{NUM}=$a[3]if defined($proteins{$a[0]});
	    $gold{$a[0]}{NUM}=0 unless defined($proteins{$a[0]});
	}
	close(G);
	$gold{READ}=1;
    }
    return unless $_[1]; ## if this is the first time stop here
    my $mode=shift;
    my $bait=shift;
    return(0) unless defined($gold{$bait});
    my ($bgo,$class_annot,$P,$class)=(0,0,0,0);
    ($bgo,$class_annot,$P,$class)=@_;
    $class='' unless $class;
    ## If we are checking BOTH annot and inter probs
    if($check_both){
	##l[0] is the inter prob and l[1] is the annot prob
	my @l=split(/_/,$P);
	push @{$gold{$bait}{MISSED}{$mode}},"$class\t$class_annot\t$bgo\t$l[0]\t$l[1]";

    }
    else{
    push @{$gold{$bait}{MISSED}{$mode}},"$class\t$class_annot\t$bgo\t$P";

    }
    
}

############################################################

sub usage{   
    $0=~/.+\/(.+)/;
    my $name = $1;
    open(HELP, "| more") ;
    print HELP <<EndOfHelp;

USAGE:  
    $name [options] <ANNOTATIONS FILE> <NETWORK FILE>

	$name will take a class annotation file and a network and
	look for interacting proteins whose GO annotations are very dissimilar.
	It will return a list of proteins with interaction probabilities 
	below a given threshold. It also requires certain files that are
	described at the end of this help.

COMMAND-LINE OPTIONS:
    -A : Check ancestor mode: discard candidate GOpairs unless all the combinations
         of their ancestors also have a low probability (modes c,m,i)
    -a : Annotated class file
    -c : Color output
    -C : Ignore ancestor annotations
    -d : Debugging output, VERY verbose
    -D : Data directory, default is current directory.
    -f : Gene Ontology GAF format annotation file (def .\/gene_association.goa_human)
    -h : Print this help and exit
    -g : Only return Gold Standard candidates 							 
    -G : geneology file, output of OntoGenealogy.pl (def .\/biological_process.genealogy)
    -i : Use interactome based probabilities. 
    -m : Mode:
           a|dis2class : look for proteins with annotations dissimilar to those
	                 of their class.
           b|bridge    : look for cases where a protein interacts with 2 partners each 
		         in different and distant classes.
           B|b_unique  : as above, but return only those proteins linking classes that are
			 not otherwise connected. That is, none of the proteins in class1
                         interacts with any of class2 and vice versa.
           c|difclass  : look for interactions between different classes
	   i|inter     : look for proteins found at the intersection of two
	                 classes annotated to distant GOs. 
           m|multi     : take multiclassed ONLY as candidates and look for
                         improbable interactions
           p|perc      : Identify the most common GO (>=50% of all interactors)
                         and look for interactors annotated to a distant GO.
	   P|pairs     : Look for proteins whose interaction partners
                         are annotated to dissimilar GOs
           d|direct    : Direct, find ALL interactions between dissimilar GOs
    -M : Minimum class cardinality for mode B. Classes with fewer members than
         the number given will be ignored. Def=0.
    -n : Ignore annotations inferred from the class file. Useful when counting
         total number of interactions involving dissimilar GOs (d). In that case 
	 we are only interested in the direct/inherited annotations and not 
         the inferred ones.
    -N : Ignore network probability file (NETWORK_NAME.prob)
    -o : Gene Ontology:
           P : biological process (default)
	   C : cellular compartment
	   F : biological function
    -p : minimum probability threshold (def= 1.162632e-10)
    -P : maximum probability threshold (def= 0.5)
    -s : Synonyms file (def .\/synonyms_human)
    -S : Species, (def: HUMAN) this is useful to avoid counting names
         such as CBL and CBL_HUMAN twice. 
    -T : GO terms file (def ./data/GO.terms_ids_obs)
    -v : Verbose output

     
FILES:

There are a number of files required to run $name. Most can 
be specified by command line options. The GO association probabilities 
cannot. The program expects to find a subdirectory .\/gostats containing
files with the various possible go pairs and their probability of being
found together. Given a GO pair such as GO:0032700 and GO:0032707, the 
format of these files is:

                         high tail prob           low tail prob
    0032700	0032707	0.999615355677336	0.000384644322664418


    Annotations File:
         This can be either a GAF file, or an annotated .clas 
	 file of the following format:

	  [CLASS: 22] 	precision : 0.2193	robustesse : 0.740	score : 0.8391
	  CA	 regulation_of_cellular_process
	  P\#	2
	  PN    CCDC6_HUMAN, DAPK1_HUMAN


    Gene Ontology GAF format annotation file:
         This is the standard format GAF file as downloaded from 
	 geneontology.org. For an explanation of the expected format, 
	 see http://www.geneontology.org/GO.format.gaf-2_0.shtml

    Geneology File:
	 This is the output of OntoGenealogy.pl, each line is divided 
	 into tab-separated columns. The first column has the parent 
	 GO term and each of the other column its children:
	 
	 child          parent          parent
	 GO\:0071555	GO:0071554	GO:0008150

    GO terms file:
	 This is a tab separated list of GO term names, their
	 corresponding GO numbers and associated ontology:
	 downloaded from http://www.geneontology.org/doc/GO.terms_and_ids

	 GO\:0000001	mitochondrion inheritance	P	

    Network File:
	 The network to be analysed. One interacting pair per line 
         and, optionally, the number of interactions in the first line.
         eg\:
            27276
	    1433B_HUMAN	1433E_HUMAN
	    
    Synonyms file:
	    Protein synonyms. File is in two tab separated columns, 
	    the first has the uniprot name and the second is a sapce 
	    separated list of its synonyms:
	    
	    1433F_HUMAN	     YWHAH YWHA1 IPI00216319 1433F_HUMAN Q04917 YWHAH

EndOfHelp
close(HELP);

    print STDERR "\n***** @_ *****\n\n" if @_;
    exit(1);

}
