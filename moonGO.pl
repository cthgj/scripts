#!/usr/bin/perl -w
use strict;
use Getopt::Std;
use IO::File;
use Switch;
use Term::ANSIColor; 
use Math::Round;
use DBI();
sub v_say;
sub debug;
require "MY_SUBS.pl";
my (%wanted,%terms_to_GOs,%gos,%not_moonlighting,%seen,%moonlighting,%ancestors,%prob,%inter_prob,%offspring,%opts,%proteins,%synonyms,%classes,%precision);
getopts('BCINVbcdeghilnvxA:D:G:H:L:M:P:S:t:T:a:f:m:p:s:',\%opts) || do { print "Invalid option, try 'moonGO.pl -h' for more information\n"; exit(1); };
usage() unless $ARGV[0];
my (@net_probs);
my $prev_tests=0;
my $no_multi=$opts{l}||undef;
my $cands_list=$opts{L}||undef;
my $very_verbose=$opts{V}||undef;
my $Aprob_file_read=0; ## is 1 if the annotation prob file for this network has been read, gopair_probability
my $Iprob_file_read=0; ## is 1 if the interactome prob file for this network has been read, interactome_probability
my $TESTS=0;
## Global MySQL variables
my ($dbh,$dsn);
my $just_count=$opts{x}||undef; ## Just count the GOs, nothing else
my $gos_seen;
my %go_counts; ## This will hold the gopairs used if -x.
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
$check_both=1 if $opts{B}; ## Check BOTH inter and annot probs, logical AND
my $either=1 unless $opts{B};## Check BOTH inter and annot probs, logical OR
my $only_gold=$opts{g}||undef;
my $geneology_file=$opts{G}||"$DATADIR/biological_process.genealogy";
my $go_terms_file=$opts{T}||"$DATADIR/GO.terms_alt_ids";
my $gold_file = "$DATADIR/gold_simple.txt";
my %gold; ## information on the gold standard proteins
our $debug=$opts{d}||undef;
my $help=$opts{h}||undef;
my $db_port=$opts{P}||3306;
my $db_host=$opts{H}||"10.1.3.30";
$db_host="127.0.0.1" if $db_host =~ /^l/i;
$db_port=3307 if $db_host eq "127.0.0.1";
my $database=$opts{A}||'pogo';
my $exp_codes=$opts{e}||undef; ## Use only "good" evidence codes
my $num_of_tests=$opts{t}||undef;
my ($TESTS1,$TESTS2)=(0,0);
my $min_class_cardinality=$opts{M}||0;
my $network_file=$ARGV[0]|| die("Need a network file\n");
my $network_prob_file; ## Network specific probabilities (annot)
my $network_inter_prob_file; ## Network specific probabilities (inter)
my %foo;
# my $low_annot_prob=$opts{p}||0.9393275;#1.193455e-07;
# my $low_inter_prob=$opts{p}||0.02463311;#1.193455e-07;
my $check_interactome_probs=$opts{i}||undef;
my $subonto=$opts{o}||"P";
our $verbose=$opts{v}||undef;
my $no_class=$opts{n}||undef; 
my $synfile;
my $species=$opts{S}||'human';
my $gaf_annotations_file;
## $added_to_prob_file will be 1 if a 
## probability is added to this network's .prob file
my $added_to_prob_file=0; 
my $added_to_inter_prob_file=0; 
my $suffix; ## Used to add extension to a name. Eg, APBB3 => APBB3_MOUSE
my $prob=$opts{p}||0.05; ## evalue significance threshold
## see sub parse_input_files_guess_values
my $read_network_prob_file;
$opts{N} && do {$read_network_prob_file=1};
#$verbose=1 if $debug;

my $check_all_class_annots_for_i=0;
my $mode=$opts{m}|| 'i';
if ($mode eq 'I') {
    $mode='i';
    $check_all_class_annots_for_i=1;
    
}
$synonyms{LOADED}=0;


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


###############################
# Collect various information #
###############################
&parse_input_files_guess_values();


my $interactome_prob_file=$opts{I}||"$stats_dir/$species.interactome.probs";
my $prob_file=$stats_dir . "/$species.prob" || $stats_dir . "/human.prob";
my $annotations_file=$opts{a}||"$DATADIR/$species.BP.clas";
my $date=localtime;
my $LOG = check_file("$DATADIR/moonlighting.log", "a");
#print "aaa : $date :  $DATADIR : $LOG\n";
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
print $LOG "Threshold   : $prob\n";
print $LOG "\n";



&load_network($network_file); 
&load_ancestors();
&parse_gaf_file();
&load_annotations($annotations_file) if $annotations_file;
&check_gold(); ## read the list of gold standard prots

#######################################################
# If we have been given a list of desired candidates  #
#######################################################
if ($cands_list) {
    my $fh=check_file($cands_list, "r");
    my $file_type=0;
    while (<$fh>) {
	if ($.==1 && /^\#\s*TESTS.+?(\d+)/) {
	    $prev_tests=$1;
	}
	next if /^\s*#/;
	chomp;
	my @a=split(/\t/);	
	$wanted{$a[0]}++;
    }
}

############################################################################
########################### MAIN PROGRAM START #############################
############################################################################

my $protnum=scalar(keys(%proteins));
my $c=0; ## counter
switch($mode) {
  
############################################################################
########### Look for proteins with at least 1 DIRECT annotation  ###########
###########          dissimilar to ALL those of the class        ###########
############################################################################
    case ['a','dis2class'] {
	my %candidates;
      bait:foreach my $bait (keys(%proteins)){
	  $c++;
	  printf STDERR ("$c of $protnum\r") if $verbose;
	  $cands_list && do {
	      next unless defined($wanted{$bait});
	  };
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
		if ($just_count) {
		    my $pair=join("_", sort  {$b lt $a} ($bgo,$class_annot));
		    $go_counts{$pair}=1;
		    $gos_seen++;
		    next bgo;
		}
		my ($test,$Pa,$Pi,$P);
		my $prec1=&term_precision($bgo);
		my $prec2=&term_precision($class_annot);
		## If we are checking BOTH annot and inter probs
		if($check_both){
		    ($Pi,$Pa)=(&interactome_probability($bgo,$class_annot),&gopair_probability($bgo,$class_annot));
		    ## Check if the probabilities are below the threshold
		    $test=&do_test($Pi,$Pa);
		    $P=$Pi . "\t" . $Pa;
		}
		## If we are NOT checking BOTH annot and inter probs
		else{
		    $check_interactome_probs ?   
			($P=&interactome_probability($bgo,$class_annot)) : 
			($P=&gopair_probability($bgo,$class_annot));	
		    $test=&do_test($P);
		    $Pa=$P;
		}		    	
		if ($test==1){
		    $is_cand=1;
		    push @{$gos{$bait}{$class}}, $bgo . "_" . $class_annot;
		    ## Make sure we keep the  least probable of the interactions
		    ## in question
		    if(defined($hits{$bait}{$class})){
			print "AAA \$hits{$bait}{$class} : $hits{$bait}{$class}\n";

			my @jj=split(/\t/,$hits{$bait}{$class});
			$jj[1]=~/\s([^\s]+)$/;
			my $p1=$1;
			$jj[2]=~/\s([^\s]+)$/;
			my $p2=$1;
			my $pp=$jj[3];
			if($p1+$p2 < $prec1+$prec2){
			    $hits{$bait}{$class}="$bait\t$bgo $prec1\t$class $class_annot $prec2\t$P";
			}
			elsif($p1+$p2 == $prec1+$prec2){
			    $hits{$bait}{$class}="$bait\t$bgo $prec1\t$class $class_annot $prec2\t$P" if $Pa<$pp;
			}
		    }
		    else{
			$hits{$bait}{$class}="$bait\t$bgo $prec1\t$class $class_annot $prec2\t$P";
		    }
		}
		else{
		    $is_cand=0;
		    ## If this prot is in the gold standard,
		    ## check WHY it was discarded1
		    &check_gold($mode,$bait,$bgo,$class_annot,$P,$class);
		}    
	    } ## end foreach $bgo
	  }## end foreach $class_annot
	    if($is_cand==1){
		$moonlighting{$bait}{$class}=$hits{$bait}{$class};
	    }
	} ## end foreach class
      } ## end foreach bait
	print STDERR "\n" if $verbose;
	$c=0;
    } ## end case 'a'
############################################################################
########### Look for interactions bridging dissimilar classes    ###########
############################################################################
    case ['B','unique_bridge','b','bridge'] 
    {
 bait:foreach my $bait (keys(%proteins)){
	 my %gg;
	 $c++;
	 my $is_cand=0;
	 printf STDERR ("$c of $protnum\r") if $verbose;
	 $cands_list && do {
	     next unless defined($wanted{$bait});
	 };
	 my @targets=keys(%{$proteins{$bait}{INTERACTORS}});
	  my @bait_classes=keys(%{$proteins{$bait}{CLASS}});
	  if(defined($gold{$bait})){push@{$gold{$bait}{MISSED}{$mode}},"UnClassed" if $#bait_classes==-1;}
	  
	  my @TARGETS=keys(%{$proteins{$bait}{INTERACTORS}});
      t1:for (my $ii=0; $ii<=$#TARGETS; $ii++) {
	      my $target1=$TARGETS[$ii];
	      #t1:foreach  my $target1(keys(%{$proteins{$bait}{INTERACTORS}})){
	      ##########################################
	      # next t1 if it shares a class with bait #
	      ##########################################
	      map{
		  push@{$gold{$bait}{MISSED}{$mode}},"Shared_Class $target1 $_" if defined($gold{$bait}); 
		  next t1 if defined($proteins{$bait}{CLASS}{$_});
	      }keys(%{$proteins{$target1}{CLASS}});
	  t2:for (my $iii=$ii; $iii<=$#TARGETS; $iii++) {
		  my $target2=$TARGETS[$iii];
		  #	  t2:foreach  my $target2(keys(%{$proteins{$bait}{INTERACTORS}})){
		  ##########################################
		  # next t2 if it shares a class with bait #
		  ##########################################
		  map{
		      push@{$gold{$bait}{MISSED}{$mode}},"Shared_Class $target2 $_" if defined($gold{$bait}); 
		      next t2 if defined($proteins{$bait}{CLASS}{$_});
		  }keys(%{$proteins{$target2}{CLASS}});
		  ######################################
		  # next t2 if t1 shares class with t2 #
		  ######################################
		  map{
		      push@{$gold{$bait}{MISSED}{$mode}},"Shared_Class $target1-$target2 $_" if defined($gold{$bait}); 
		      next t2 if defined($proteins{$target1}{CLASS}{$_});
		  }keys(%{$proteins{$target2}{CLASS}});
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
			  my $class_name=join("_", sort  {$b lt $a} ($class1,$class2));
			  if ($mode eq 'B') {
			      ################################################
			      # Check whether these classes are connected by #
			      # any other proteins.			     #
			      ################################################
			      my $conX=class_connected($bait,$class1, $class2);
			      if ($conX>0){
				  $not_moonlighting{$bait}{$class_name}++;
				  next class2 ;
			      }
			  }
		      c1:foreach my $class1_annot (@uniqu1){
			      next c1 if $class1_annot eq 'GO:0008150';
			      my $prec1=&term_precision($class1_annot);
			  c2:foreach my $class2_annot (@uniqu2){
				  next c2 if $class2_annot eq 'GO:0008150';
			if ($just_count) {
			    my $pair=join("_", sort  {$b lt $a} ($class1_annot,$class2_annot));
			    $gos_seen++;
			    $go_counts{$pair}=1;
			    next c2;
			}
			my $prec2=&term_precision($class2_annot);
			my ($test,$Pa,$Pi,$P);
			## If we are checking BOTH annot and inter probs
			if ($check_both) {
			    ($Pi,$Pa)=(&interactome_probability($class1_annot,$class2_annot),&gopair_probability($class1_annot,$class2_annot));
			    ## Check if the probabilities are below the threshold
			    $test=&do_test($Pi,$Pa);
			    $P=$Pi . "\t" . $Pa;
			}
			## If we are NOT checking BOTH annot and inter probs
			else {
			    $check_interactome_probs ?   
				($P=&interactome_probability($class1_annot,$class2_annot)) : 
				    ($P=&gopair_probability($class1_annot,$class2_annot));	
			    $test=&do_test($P);
			}
			if ($test==1) {
			    push @{$gos{$bait}{$class_name}},  $class1_annot . "_" . $class2_annot ;
			    ## Make sure we keep the MOST specific LEAST probable of the interactions
			    ## in question
			    if (defined($moonlighting{$bait}{$class_name})) {
				$moonlighting{$bait}{$class_name}=~/^.+?\t.+?\s.+?([\.\d]+)\t.+?\s.+?([\.\d]+).+\t(.+?)$/;
				my ($p1,$p2,$pp)=($1,$2,$3);
				if ($p1+$p2 < $prec1+$prec2) {
				    $moonlighting{$bait}{$class_name}=save_cand($bait,$class1,$class2,$class1_annot,$class2_annot,$prec1,$prec2, $P);
				} 
				elsif ($p1+$p2 == $prec1+$prec2) {
				    $moonlighting{$bait}{$class_name}=save_cand($bait,$class1,$class2,$class1_annot,$class2_annot,$prec1,$prec2, $P) if $pp<=$P;   
				}
			    } 
			    else {
				$moonlighting{$bait}{$class_name}=save_cand($bait,$class1,$class2,$class1_annot,$class2_annot,$prec1,$prec2, $P);
			    }
			}
			#next bait; 
			else {
			    my $cc=$class1 . "-" . $class2;
			    &check_gold($mode,$bait,$class1_annot,$class2_annot,$P,$cc);
			    next class2
			}
		  
		  
		    }		# end class2_annot
			  }	#end class1_annot
		      }		#end class2
		  }		## end class1 
	      }			#end t2
	  }			## end t1
      }				## end foreach bait
					     
# ## Now find those candidates that connect classes 
# ## that are NOT connected by any other interactions
# 	if($mode eq 'B'){
# 	  cand:foreach my $cand (keys(%moonlighting)){
# 	    classname:foreach my $classname (keys(%{$moonlighting{$cand}})){
# 		$moonlighting{$cand}{$classname}=~/^$cand\t(.+?)\s+(.+?)\s+(.+?)\t(.+?)\s+(.+?)\s+(.+?)\t(.+)$/||die("Could not match cand : \$moonlighting{$cand}{$classname} : $moonlighting{$cand}{$classname}\n");
# 		my ($class1,$class1_annot,$prec1,$class2,$class2_annot,$prec2,$P)=($1,$2,$3,$4,$5);
# 		## Make sure there are no connections
# 		## between any prot in class1 and any in class2
# 		my @prots1= keys(%{$classes{$class1}{PROTS}}) ;
# 		my @prots2= keys(%{$classes{$class2}{PROTS}}) ;
# 	      p1:foreach my $prot1 (@prots1){
# 		p2:foreach my $prot2 (@prots2){
# 		    if(defined($proteins{$prot1}{INTERACTORS}{$prot2})){
# 			print "AAAAA : $prot1 <=> $prot2\n";
# 			push@{$gold{$cand}{MISSED}{$mode}},"$class1\t$class2\tConnected" if defined($gold{$cand});
# 			$not_moonlighting{$cand}{$classname}++;
# 			next classname;
# 		    }
# 		}
# 	      }
# 	    } ## end foreach $classname
# 	  } ## end foreach $cand
# 	}
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
	  $cands_list && do {
	      next unless defined($wanted{$bait});
	  };
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
	    foreach my $bgo (keys(%{$proteins{$bait}{GOs}})){
		next if $bgo eq 'GO:0008150';
		next unless $proteins{$bait}{GOs}{$bgo} eq 'DIR'; 
		foreach my $tgo (keys(%{$proteins{$target}{GOs}})){
		    next if $tgo eq 'GO:0008150';
		    next unless $proteins{$target}{GOs}{$tgo} eq 'DIR'; 
		    if ($just_count) {
			my $pair=join("_", sort  {$b lt $a} ($bgo,$tgo));
			$gos_seen++;
			$go_counts{$pair}=1;
			next;
		    }

		    my ($test,$Pa,$Pi,$P);
		    ## If we are checking BOTH annot and inter probs
		    if($check_both){
			($Pi,$Pa)=(&interactome_probability($bgo,$tgo),&gopair_probability($bgo,$tgo));
			## Check if the probabilities are below the threshold
			$test=&do_test($Pi,$Pa);
			$P=$Pi . "\t" . $Pa;
		    }
		    ## If we are NOT checking BOTH annot and inter probs
		    else{
			$check_interactome_probs ?   
			    ($P=&interactome_probability($bgo,$tgo)) : 
			    ($P=&gopair_probability($bgo,$tgo));	
			$test=&do_test($P);
			$Pa=$P;
		    }
		    if($test==1){
			push @{$gos{$bait}{$target}}, $bgo  . "_" . $tgo;
			## Make sure we keep the LEAST probable of the interactions
			## in question	 
			if(defined($moonlighting{$bait}{$target})){
			    $moonlighting{$bait}{$target}=~/.+\t(.+?)$/;
			    my $koko=$1;
			    $moonlighting{$bait}{$target}="$bait ($bgo)  \t:\t$target ($tgo)  \t$P" if $Pa<$koko;
			}
			else{
			    $moonlighting{$bait}{$target}="$bait ($bgo)  \t:\t$target ($tgo)  \t$P";
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
	} ## end foreach target
      } ## end foreach bait
    } ## end case "c"

############################################################################
###########       Identify improbable annotations                ###########
############################################################################
    # case ['D', 'Direct']{
    #   bait:foreach my $bait (keys(%proteins)){
    # 	  $c++;
    # 	  printf STDERR ("$c of $protnum\r") if $verbose;
    # 	  &check_gold($bait);
    # 	  my %bait_gos;
    # 	  ## Ignore class annotations
    # 	  map{$bait_gos{$_}=$_ if $proteins{$bait}{GOs}{$_} eq 'DIR' }keys(%{$proteins{$bait}{GOs}});
    # 	  my @bgos=keys(%bait_gos);
    # 	bgo:for (my $n=0; $n<=$#bgos; $n++){
    # 	    next bgo if $bgos[$n] eq 'GO:0008150';
    # 	  tgo:for (my $k=$n; $k<=$#bgos; $k++){
    # 	      my ($bgo,$tgo)=($bgos[$n],$bgos[$k]);
    # 	      next tgo if $tgo eq 'GO:0008150';
    # 	      my ($test,$Pa,$Pi,$P);
    # 	      #my $prec1=&term_precision($bgo);
    # 	      #my $prec2=&term_precision($tgo);

    # 	      my $prec1=1;
    # 	      my $prec2=1;
    # 	      ## If we are checking BOTH annot and inter probs
    # 	      if($check_both){
    # 		  ($Pi,$Pa)=(&interactome_probability($bgo,$tgo),&gopair_probability($bgo,$tgo));
    # 		  ## Check if the probabilities are below the threshold
    # 		  $test=&do_test($Pi,$Pa);
    # 		  $P=$Pi . "\t" . $Pa;
    # 	      }
    # 	      ## If we are NOT checking BOTH annot and inter probs
    # 	      else{
    # 		  $check_interactome_probs ?   
    # 		      ($P=&interactome_probability($bgo,$tgo)) : 
    # 		      ($P=&gopair_probability($bgo,$tgo));	
    # 		  $test=&do_test($P);
    # 		  $Pa=$P;
    # 	      }
    # 	      if($test==1){
    # 		   $moonlighting{$bait}="$bait\t$bgo $prec1\t$tgo $prec2\t$P";
    # 		  ## Make sure we keep the MOST specific, LEAST probable of the interactions
    # 		  ## in question	 
    # 		  # if(defined($moonlighting{$bait})){
    # 		  #     $moonlighting{$bait}=~/.+?GO:\d+\s(.+?)\tGO:\d+\s(.+?)\t(.+?)$/;
    # 		  #     my ($p1,$p2,$koko)=($1,$2,$3);
    # 		  #     ## If this is the most precise case, keep
    # 		  #     if($p1+$p2 > $prec1+$prec2){
    # 		  # 	  $moonlighting{$bait}="$bait\t$bgo $prec1\t$tgo $prec2\t$P";
    # 		  #     }
    # 		  #     ## If it is AS precise, keep if less likely
    # 		  #     elsif($p1+$p2 == $prec1+$prec2){
    # 		  # 	  $moonlighting{$bait}="$bait\t$bgo $prec1\t$tgo $prec2\t$P" if $Pa<$koko;
    # 		  #     }
    # 		  # }
    # 		  # else{
    # 		  #     $moonlighting{$bait}="$bait\t$bgo $prec1\t$tgo $prec2\t$P";
    # 		  # }
    # 		  next bait;
    # 	      } ## end test==1
    # 	      # else{
    # 	      # 	  &check_gold($mode,$bait,$bgo,$tgo,$P);
    # 	      # 	  $not_moonlighting{$bait}=1;
    # 	      # }
    # 	  }## end foreach tgo
    # 	}## end foreach bgo

    #   }## end bait
    # }## end case D


############################################################################
###########       Identify ALL improbable interactions           ###########
############################################################################
    case ['d','direct'] {
      bait:foreach my $bait (keys(%proteins)){
	  $c++;
	  printf STDERR ("$c of $protnum\r") if $verbose;
	  $cands_list && do {
	      next unless defined($wanted{$bait});
	  };
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
		  next tgo if $tgo eq 'GO:0008150';
	      ## skip target if the bait protein has even one of the target's GOs
	      next target if exists($bait_gos{$tgo});
	    bgo:foreach my $bgo (keys(%bait_gos)){
		next bgo if $bgo eq 'GO:0008150';
		if ($just_count) {
		    my $pair=join("_", sort  {$b lt $a} ($bgo,$tgo));
		    $gos_seen++;
		    $go_counts{$pair}=1;
		    next bgo;
		}
		my ($test,$Pa,$Pi,$P);
		my $prec1=&term_precision($bgo);
		my $prec2=&term_precision($tgo);
		## If we are checking BOTH annot and inter probs
		if($check_both){
		    ($Pi,$Pa)=(&interactome_probability($bgo,$tgo),&gopair_probability($bgo,$tgo));
		    ## Check if the probabilities are below the threshold
		    $test=&do_test($Pi,$Pa);
		    $P=$Pi . "\t" . $Pa;
		}
		## If we are NOT checking BOTH annot and inter probs
		else{
		    $check_interactome_probs ?   
			($P=&interactome_probability($bgo,$tgo)) : 
			($P=&gopair_probability($bgo,$tgo));	
		    $test=&do_test($P);
		    $Pa=$P;
		}
		if($test==1){
		    ## Useful for counting total number of interactions
		    ## involving dissimilar GOs. In that case, we are only
		    ## interested in DIRECT/INHERITED annotations, no INFERRED
		    if($no_class){
			next bgo if $proteins{$bait}{GOs}{$bgo} eq "INF";
			next tgo if $proteins{$target}{GOs}{$tgo} eq "INF";

		    }
		    
		    push @{$gos{$bait}{$target}},  $bgo . "_" . $tgo ;
		    ## Make sure we keep the MOST specific, LEAST probable of the interactions
		    ## in question	 
		    if(defined($moonlighting{$bait}{$target})){
			$moonlighting{$bait}{$target}=~/.+?GO:\d+\s(.+?)\t.+?GO:\d+\s(.+?)\t(.+?)$/;
			my ($p1,$p2,$koko)=($1,$2,$3);
			## If this is the most precise case, keep
			if($p1+$p2 > $prec1+$prec2){
			    $moonlighting{$bait}{$target}="$bait\t$bgo $prec1\t$target\t$tgo $prec2\t$P";
			}
			## If it is AS precise, keep if less likely
			elsif($p1+$p2 == $prec1+$prec2){
			    $moonlighting{$bait}{$target}="$bait\t$bgo $prec1\t$target\t$tgo $prec2\t$P" if $Pa<$koko;
			}
		    }
		    else{
			$moonlighting{$bait}{$target}="$bait\t$bgo $prec1\t$target\t$tgo $prec2\t$P";
		    }
		    
		    next target;
		}
		else{
		    &check_gold($mode,$bait,$bgo,$tgo,$P);
		    $not_moonlighting{$bait}{$target}=1;
		    next target;
		}
	    } ## end foreach my $bgo
	  } ## end foreach my $tgo
	} ## end foreach target
      } ## end foreach bait
    } ## end case "d"
############################################################################
###########       Look for proteins found at the intersection    ###########
###########       of two classes annotated to distant GOs        ###########
############################################################################
    case ['i','inter'] 
    {    
     die("Need an annotated class file for mode 'i'\n") unless $annotations_file;
      bait:foreach my $bait (keys(%proteins)){
         $c++;
	 my %seen_classes;
	 my $is_cand=0;
	 printf STDERR ("$c of $protnum\r") if $verbose;
	 $cands_list && do {
	     next unless defined($wanted{$bait});
	 };
	  if ((scalar(keys(%{$proteins{$bait}{CLASS}})) == 1) && (defined($gold{$bait}))){
	      push@{$gold{$bait}{MISSED}{$mode}},"MonoClassed";
	      next;
	  }
	      my @bait_class_annots=();
	      ################################################
              # @cls == list of classes that BAIT belongs to #
              ################################################
	      my @cls=keys(%{$proteins{$bait}{CLASS}});
	    c1:for(my $n=0; $n<=$#cls;$n++){ ## bait class 1
	      c2:for(my $k=$n+1; $k<=$#cls;$k++){ ## bait class 2
		c3:foreach my $class1_annotation (@{$classes{$cls[$n]}{ANNOT}}){ ## annot class 1
		    if($class1_annotation=~/0008150/){
			push@{$gold{$bait}{MISSED}{$mode}},"0008150 $cls[$n]" if defined($gold{$bait}); 
			next c3;
		    }
		    my $prec1=&term_precision($class1_annotation);
		  c4:foreach my $class2_annotation (@{$classes{$cls[$k]}{ANNOT}}){ ## annot class 2	  
			 # print STDERR "$n,$k : $class1_annotation : $class2_annotation\n";
		  next c4 if $cls[$k] eq $cls[$n];
		  if($class2_annotation=~/0008150/){
		      push @{$gold{$bait}{MISSED}{$mode}},"0008150 $cls[$k]" if defined($gold{$bait});
		      next c4;
		  }
		  ############################################################################
                  # If we are just counting, the total number of tests, not skipping 	     #
		  # anything (-x)							     #
                  ############################################################################
		  if ($just_count) {
		      my $pair=join("_", sort  {$b lt $a} ($class1_annotation,$class2_annotation));
		      $gos_seen++;
		      $go_counts{$pair}=1;
		      next c4;
		  }
		  my $prec2=&term_precision($class2_annotation);		  
		  my ($test,$Pa,$Pi,$P);
		  ## If we are checking BOTH annot and inter probs
		  if($check_both){
		     
		      ###############################
		      #my $pair=join("_", sort  {$b lt $a} ($class1_annotation,$class2_annotation));
		      #next unless $pair eq 'GO:0016070_GO:0050896';
		      ###############################

		      ($Pi,$Pa)=(&interactome_probability($class1_annotation,$class2_annotation),&gopair_probability($class1_annotation,$class2_annotation));
		      ## Check if the probabilities are below the threshold
		      $test=&do_test($Pi,$Pa);
		      $P=$Pi . "\t" . $Pa;
		  }
		  ## If we are NOT checking BOTH annot and inter probs
		  else{
		      $check_interactome_probs ?   
			  ($P=&interactome_probability($class1_annotation,$class2_annotation)) : 
			  ($P=&gopair_probability($class1_annotation,$class2_annotation));	
		      $test=&do_test($P);
		      #print STDOUT "$class1_annotation :: $class2_annotation :: $P :: $test\n";
		      $Pa=$P;
		  }
                  #####################################
                  # If this prob passes the threshold #
                  #####################################
		  if($test==1){ 
		      $a=join("_",sort($cls[$n], $cls[$k]));
		      $not_moonlighting{$bait}{$a}=undef;
		      push @{$gos{$bait}{$a}},  $class1_annotation . "_" . $class2_annotation;
		      my $term1=$class1_annotation;
		      my $term2=$class2_annotation;
		      
		      my $b=join("_",sort($term1,$term2));
		      #print  "Accepting $bait for $b $a $P\n"; 
		      
		      ## Make sure we keep the MOST specific LEAST probable of the interactions
		      ## in question
		      if(defined($moonlighting{$bait}{$a})){
			  $moonlighting{$bait}{$a}=~/^.+?\t.+?GO:.+?\s([\.\d]+).+?\t.+?GO:.+?\s([\.\d]+)/ || 
			      die("A1A $moonlighting{$bait}{$a}\n");
			  my ($p1,$p2)=($1,$2);
			  $moonlighting{$bait}{$a}=~/.+\t(.+?)$/ || die("AA2 $moonlighting{$bait}{$a}\n");
			  my $pp=$1;
			  ## If this is the most precise case, keep
			  if($p1+$p2 < $prec1+$prec2){
			      $moonlighting{$bait}{$a}=save_cand($bait,$cls[$n],$cls[$k],$term1,$term2,$prec1,$prec2, $P);
			  }
			  ## If it is AS precise, keep if less likely
			  elsif($p1+$p2 == $prec1+$prec2){
			      $moonlighting{$bait}{$a}= save_cand($bait,$cls[$n],$cls[$k],$term1,$term2,$prec1,$prec2, $P) if $pp<=$Pa;
			  }
		      }
		      else{
			  $a=join("_",sort($cls[$n], $cls[$k]));
			  $moonlighting{$bait}{$a}=save_cand($bait,$cls[$n],$cls[$k],$term1,$term2,$prec1,$prec2, $P);
		      }

		  }
		  else{
		      &check_gold($mode,$bait,$class1_annotation,$class2_annotation,$P,0);
		      ######################################################################
                      # Get only those candidates ALL of whose classes' annotations	   #
                      # are disimilar.  						   #
                      ######################################################################
                      if ($check_all_class_annots_for_i) {
			  my $a=join("_",sort($cls[$n], $cls[$k]));
			  $not_moonlighting{$bait}{$a}++;
			  $moonlighting{$bait}{$a}=undef;
		      }
		      $a=join("_",sort($cls[$n], $cls[$k]));
		      my $b=join("_",sort($class1_annotation,$class2_annotation));
		      #print  "Skipping $bait for $b $a $P\n"; 
		      next c2; ## next bait class			  
		  }
	      } ## end c4
	    } ## end c3
	  } ## end c2
	} ## end c1
      } ## end foreach bait
    } ## end case i
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
	      $cands_list && do {
		  next unless defined($wanted{$bait});
	      };
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
		printf STDERR ("$c of $multi_num\r")  if $verbose;
		if ($just_count) {
		    my $pair=join("_", sort  {$b lt $a} ($bgo,$tgo));
		    $gos_seen++;
		    $go_counts{$pair}=1;
		    next tgo;
		}
		my ($test,$Pa,$Pi,$P);
		## If we are checking BOTH annot and inter probs
		if($check_both){
		    ($Pi,$Pa)=(&interactome_probability($bgo,$tgo),&gopair_probability($bgo,$tgo));
		    ## Check if the probbilities are below the threshold
		    $test=&do_test($Pi,$Pa);
		    $P=$Pi . "\t" . $Pa;
		}
		## If we are NOT checking BOTH annot and inter probs
		else{
		    $check_interactome_probs ?   
			($P=&interactome_probability($bgo,$tgo)) : 
			($P=&gopair_probability($bgo,$tgo));	
		    $test=&do_test($P);
		    $Pa=$P;
		}
		if($test==1){
		    push @{$gos{$bait}{$target}},  $bgo . "_" . $tgo;
		    
		    my $term1=&terms_to_GOs($bgo,1) . ",";
		    my $term2=&terms_to_GOs($tgo,1) . ",";
		    while (length($term1)<50) {
			$term1.=" ";
		    }
		    while (length($term2)<50) {
			$term2.=" ";
		    }
		    ## Make sure we keep the LEAST probable of the interactions
		    ## in question	 
		    if (defined($moonlighting{$bait}{$target})) {
			$moonlighting{$bait}{$target}=~/.+\t(.+?)$/;
			my $koko=$1;
			$moonlighting{$bait}{$target}="$bait ($term1 $bgo,$proteins{$bait}{GOs}{$bgo})  \t:\t$target ($term2 $tgo,$proteins{$target}{GOs}{$tgo})  \t$P" if $Pa<$koko;
		    } else {
			$moonlighting{$bait}{$target} ="$bait ($term1 $bgo,$proteins{$bait}{GOs}{$bgo})  \t:\t$target ($term2 $tgo,$proteins{$target}{GOs}{$tgo})  \t$P"; 
		    }
		    next target;
		    
		}
		# Discard target if ANY of its GOs are related
		# to ANY of the bait's.
		else{
		    &check_gold($mode,$bait,$bgo,$tgo,$P);
		    $not_moonlighting{$bait}{$target}=1;
		    next target;
		}
	    } ## end foreach tgo
	  } ## end foreach bgo
	} ## end foreach target
      } ## end foreach bait
    } ## end case "m"
############################################################################
########### If the bait interacts with >=50% with a specific GO, ###########
########### identify interactions with another                   ###########
############################################################################
    case ['p','perc'] {
      bait:foreach my $bait (keys(%proteins)){
	  &check_gold($bait);
	  $c++;
	  printf STDERR ("$c of $protnum\r") if $verbose;
	  $cands_list && do {
	      next unless defined($wanted{$bait});
	  };
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
	  $cands_list && do {
	      next unless defined($wanted{$bait});
	  };
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
		  if ($just_count) {
		      my $pair=join("_", sort  {$b lt $a} ($go1,$go2));
		      $gos_seen++;
		      $go_counts{$pair}=1;
		      next go2;
		  }
		  my ($test,$Pa,$Pi,$P);
		  if($check_both){
		      ($Pi,$Pa)=(&interactome_probability($go1,$go2),&gopair_probability($go1,$go2));
		      ## Check if the probbilities are below the threshold
		      $test=&do_test($Pi,$Pa);
		      $P=$Pi . "\t" . $Pa;
		  }
		  ## If we are NOT checking BOTH annot and inter probs
		  else{
		      $check_interactome_probs ?   
			  ($P=&interactome_probability($go1,$go2)) : 
			  ($P=&gopair_probability($go1,$go2));	
		      $test=&do_test($P);
		      $Pa=$P;
		  }
		  if($test==1){
		      push @{$gos{$bait}}, $go1 . "_" . $go2;
		      my $last_p=10;
		      if ($is_cand=~/.+\t(.+?)$/){$last_p=$1;}
		      $is_cand="($go1)\t$target2 ($go2)\t$P" if $Pa<$last_p;
		  }
		  else{
		      &check_gold($mode,$bait,$target1 . "-" . $target2,$go1 . "-" . $go2,$P);
		      $is_cand='nope';
		      next t2;
		  }
		  

		  if($is_cand !~ /^nope$/){
		      $moonlighting{$bait}="$bait\t$target1 $is_cand";
		  }
	      }
	    }
	  }
	}
      } ## end foreach bait
    } ## end case 'P'


############################################################################
############################################################################
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
} ## end switch;

############################################################################
########################### MAIN PROGRAM END ###############################
############################################################################

########################### Print Results ########################

#########################################
# If we are just counting the GO pairs  #
#########################################
if ($just_count) {
    print STDERR "\n" if $verbose;
    print "Mode $mode:\t" .  scalar keys(%go_counts) . " unique GO pairs and $gos_seen tests performed\n";
    exit(0);
}




my %header=(
	    "i" => "Candidate\tClass1\tClass1 Card.\tgo1\tgo1 prec.\tClass2\tClass2 Card.\tgo2\tgo2 prec.\tprobability",
	   );
$header{"I"}=$header{"b"}=$header{"B"}=$header{i};
if ($check_both) {
    $header{"i"}="Candidate\tClass1\tClass1 Card.\tgo1\tgo1 prec.\tClass2\tClass2 Card.\tgo2\tgo2 prec.\tInter. e-val\tAssoc. e-val";
    $header{"I"}=$header{"b"}=$header{"B"}=$header{i};

}
my @oo=keys(%moonlighting);
my $nn=scalar(@oo);
my $pp=1;
my $color="bold blue";
my %output;
my %gfound;
# If we are checking both probabilities, the number of tests performed will
# be doubled. Divide $TESTS/2 to avoid this artificial stringency
$TESTS=$TESTS/2 if $check_both;
print STDERR "\n" . scalar keys(%moonlighting) .  " candidates found (before MT correction).\n" if $verbose;

##################################################
# If we have been given a candidate list and if  #
# that list set a value for $TESTS, use that	 #
##################################################
$cands_list && do {
    $TESTS=$prev_tests if $prev_tests>0;
	
};
if ($no_multi){
     print STDERR "\nSkipping multiple testing correction.\n" if $verbose;
     $TESTS=1;
}
else {
    print STDERR "\nMultiple testing correction...\n" if $verbose;
}
if (@oo) {
    print STDOUT "# TESTS : $TESTS\n"; 
    print STDOUT "# $header{$mode}\n";
}
else {
    exit(0);
}
    print STDERR "# TESTS : $TESTS\nAnnot. : $TESTS1\nInter : $TESTS2\n" if $verbose; 

my %cands_kept;
    $TESTS=$num_of_tests if $num_of_tests;
bait:foreach my $bait (@oo){
    print STDERR "$pp of $nn\r" if $verbose;
    my $is_gold=$gold{$bait}||undef;
    $only_gold && do {next unless $is_gold; };
    my $P;					      
    if(  ($mode eq 'pairs') || ($mode eq 'P') || ($mode eq 'D') || ($mode eq 'Direct') ){
	next if defined($not_moonlighting{$bait});
	## Now that we have the number of $TESTS performed,
	## go through each of the go pairs annotating this cad/target pair
	## and only keep cand if ALL the CORRECTED pair probabilities
	## pass the threshold 
	my ($result,$Pi_corr,$Pa_corr,$P_corr);
	foreach my $gop (@{$gos{$bait}}){
	    my ($go1,$go2)=split(/_/,$gop);
	    if($check_both){
		my ($Pi,$Pa)=(&interactome_probability($go1,$go2,'nocount'),&gopair_probability($go1,$go2,'nocount'));
		($result,$Pi_corr,$Pa_corr)=&do_test($Pi,$Pa,$TESTS);
		next bait unless $result==1;
	    }
	    else{
		my $PP;
		$check_interactome_probs ?   
		    ($PP=&interactome_probability($go1,$go2,'nocount')) : 
		    ($PP=&gopair_probability($go1,$go2,'nocount'));
		($result,$P_corr)=&do_test($PP,$TESTS);
		next bait unless $result==1;
	    }
	}
	if($check_both){
	    $moonlighting{$bait}=~/.+\t(.+?\t.+?)$/||die("Could not match mprob : $moonlighting{$bait}\n");
	    $P=$1;
	    $moonlighting{$bait}=~s/$P$/$Pi_corr\t$Pa_corr/;
	    if($is_gold){ 
		$moonlighting{$bait}=~s/$bait/color("$color").$bait.color("reset")/ige if $opts{c};
		$gfound{$bait}++;
	    }
	    push @{$output{$Pa_corr}},$moonlighting{$bait};
	    $cands_kept{$bait}++;
	}
	else{
	    $moonlighting{$bait}=~/.+\t(.+?)$/||die("Could not match mprob1 : $moonlighting{$bait}\n");
	    $P=$1;
	    my $P_corr=$P*$TESTS;
	    $P_corr=1 if $P_corr>1;
	    print "$moonlighting{$bait}\n";
	    next unless $P_corr<=$prob;
	    $moonlighting{$bait}=~s/$P$/$P_corr/;
	    if($is_gold){ 
		$moonlighting{$bait}=~s/$bait/color("$color").$bait.color("reset")/ige if $opts{c};
		$gfound{$bait}++;
	    }
	    push @{$output{$P_corr}},$moonlighting{$bait};
	    $cands_kept{$bait}++;
	}
    }
    ## For all other modes
    else{
	target:foreach my $target (keys(%{$moonlighting{$bait}})){
	    next if defined($not_moonlighting{$bait}{$target});
	     ##########################################################################
             # Now that we have the number of $TESTS performed,			      #
	     # go through each of the go pairs annotating this cand/target pair	      #
	     # and only keep cand if ALL the CORRECTED pair probabilities	      #
	     # pass the threshold 						      #
             ##########################################################################
	    my ($result,$Pi_corr,$Pa_corr,$P_corr);
		foreach my $gop (@{$gos{$bait}{$target}}){
		    my ($go1,$go2)=split(/_/,$gop);
		    if($check_both){
			my ($Pi,$Pa)=(&interactome_probability($go1,$go2,'nocount'),&gopair_probability($go1,$go2,'nocount'));
			($result,$Pi_corr,$Pa_corr)=&do_test($Pi,$Pa,$TESTS);
			$Pi_corr=1 if $Pi_corr>1;
			$Pa_corr=1 if $Pa_corr>1;
			next target unless $result==1;
		    }
		    else{
			my $PP;
			$check_interactome_probs ?   
			    ($PP=&interactome_probability($go1,$go2,'nocount')) : 
			    ($PP=&gopair_probability($go1,$go2,'nocount'));
			($result,$P_corr)=&do_test($PP,$TESTS);
			$P_corr=1 if $P_corr>1;
			next target unless $result==1;
		    }
		}
	    ## an 2.33483243394472e-11, in 2.33483243394472e-11
	    ## If we are checking BOTH annot and inter probs
	    if($check_both){
		$moonlighting{$bait}{$target}=~/.+\t(.+?\t.+?)$/ || 
		    die("Could not match mprob2 $bait:$target: $moonlighting{$bait}{$target}\n");
		$P=$1;
		my ($Pi,$Pa)=split(/\t/,$P);
		## Correct the probs
		my $Pi_corr=$Pi*$TESTS;
		my $Pa_corr=$Pa*$TESTS;
		die("Something wrong with multi testing\n") if $Pa_corr>$prob;
		die("Something wrong with multi testing\n") if $Pi_corr>$prob;
		my $PPP=$Pi_corr . "\t" . $Pa_corr;
		$moonlighting{$bait}{$target}=~s/$P$/$PPP/;
		if($is_gold){ 
		    $moonlighting{$bait}{$target}=~s/$bait/color("$color").$bait.color("reset")/ige if $opts{c};
		    $gfound{$bait}++;
		}
		push @{$output{$Pa_corr}},$moonlighting{$bait}{$target};
		$cands_kept{$bait}++;
	    }
	    ## If we are NOT checking both inter and annot probs
	    else{
		$moonlighting{$bait}{$target}=~/.+\t\s*(.+?)$/ || 
		    die("Could not match mprob3 $bait:$target: $moonlighting{$bait}{$target}\n");
		$P=$1;
		my $P_corr=$P*$TESTS;
		$moonlighting{$bait}{$target}=~s/$P$/$P_corr/;
		if($is_gold){ 
		    $moonlighting{$bait}{$target}=~s/$bait/color("$color").$bait.color("reset")/ige if $opts{c};
		    $gfound{$bait}++;
		}
		push @{$output{$P_corr}},$moonlighting{$bait}{$target};
		$cands_kept{$bait}++;
	    }
	}
    }
    $pp++;
}

my @sorted=sort{$a <=> $b} keys(%output);
$c=0;
map{
    map{
	print "$_\n";
    }@{$output{$_}};
}@sorted;
print STDERR "\n" . scalar keys(%cands_kept) .  " candidates found (after MT correction).\n" if $verbose;

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
chomp($cmnd);
print G "\n$cmnd : " . scalar(@ff) . " gold standard proteins were found: @ff\n";
close(G);
#################### Create network specific probability file ####
if($read_network_prob_file){
    unless($added_to_prob_file==0){
	print STDERR "Generating prob file...\n" if $verbose;
	my $NN = $network_file . ".prob";
	open(NET,">>$NN");
	map{print NET}@net_probs;
	system("cat $NN | sort | uniq > lililolo.$$; mv lililolo.$$ $NN");
	close(NET);
    }
    unless($added_to_inter_prob_file==0){
	print STDERR "Generating Inter prob file...\n" if $verbose;
	my $NN = $network_file . ".inter.prob";
	open(NET,">>$NN");
	map{print NET}@net_probs;
	system("cat $NN | sort | uniq > lililolo.$$; mv lililolo.$$ $NN");
	close(NET);
    }
}


################################ END ##############################

############################### SUBROUTINES ################################
sub parse_input_files_guess_values{
    ## Change long mode names to short
    if($mode eq 'dis2class'){$mode = 'a'}
    elsif($mode eq 'bridge'){$mode = 'b'}
    elsif($mode eq 'b_unique'){$mode = 'B'}
    elsif($mode eq 'difclass'){$mode = 'c'}
    elsif($mode eq 'inter'){$mode = 'i'}
    elsif($mode eq 'multi'){$mode = 'm'}
    elsif($mode eq 'perc'){$mode = 'p'}
    elsif($mode eq 'pairs'){$mode = 'P'}
    elsif($mode eq 'direct'){$mode = 'd'}

    ## Guess species, assign variables accordingly
    $species=guess_species($network_file) unless $opts{S};
    unless (defined ($species)) {
	if(($network_file=~/^dro/) || 
	   ($network_file =~/^dm/i) || 
	   ($network_file =~/^fly/i)){$species="fly"}
	elsif(($network_file=~/^hs/) || 
	      ($network_file=~/^hs/i) || 
	      ($network_file =~/^hum/i)){$species="human"}
	elsif(($network_file=~/^mus/) || 
	      ($network_file=~/^mm/i)|| 
	      ($network_file =~/^mou/i)){$species="mouse"}
	elsif(($network_file=~/^scc/) || 
	      ($network_file=~/^sc/i)|| ($network_file =~/^yea/i)){$species="yeast"}
	elsif(($network_file=~/^ele/) || 
	      ($network_file=~/^ce/i)|| ($network_file =~/^wor/i)){$species="worm"}
    }
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
	$synfile=$opts{s}||"$DATADIR/human.map";
	$prob_file=$stats_dir . "/human.prob";
	$gaf_annotations_file=$opts{f} || "$DATADIR/gene_association.human";
	$suffix='HUMAN';
	# $low_annot_prob= 0.9393275 unless $opts{p};
	# $low_inter_prob= 0.9157528 unless $opts{p};    
    }
    elsif(($species eq "fly")   || ($species eq "fb")  || ($species eq "dm") || ($species eq "dro")){
	$species='fly';
	$stats_dir="$DATADIR/gostats/fly";
	$synfile=$opts{s}||"$DATADIR/fly.map";
	$prob_file=$stats_dir . "/fly.prob";
	$gaf_annotations_file=$opts{f} || "$DATADIR/gene_association.fly";
	$suffix='DROME';
 	# $low_annot_prob= 0.8450468 unless $opts{p};
	# $low_inter_prob= 0.9920737 unless $opts{p};    
    }
    elsif(($species eq "worm")  || ($species eq "wb")  || ($species eq "ce") || ($species eq "ele")){
	$species='worm';
	$stats_dir="$DATADIR/gostats/worm";
	$synfile=$opts{s}||"$DATADIR/worm.map";
	$prob_file=$stats_dir . "/worm.prob";
	$gaf_annotations_file=$opts{f} || "$DATADIR/gene_association.worm";
	$suffix='CAEEL';
	# $low_annot_prob= 0.8819694 unless $opts{p};
	# $low_inter_prob= 0.9900118 unless $opts{p};    

    }
    elsif(($species eq "mouse") || ($species eq "mg")  || ($species eq "mm") || ($species eq "mus")){
	$species='mouse';
	$stats_dir="$DATADIR/gostats/mouse";
	$prob_file=$stats_dir . "/mouse.prob";
	$synfile=$opts{s}||"$DATADIR/mouse.map";
	$gaf_annotations_file=$opts{f} || "$DATADIR/gene_association.mouse";
	$suffix='MOUSE';
	# $low_annot_prob= 0.9540235 unless $opts{p};
	# $low_inter_prob= 0.9948855 unless $opts{p};    

    }
    elsif(($species eq "yeast") || ($species eq "sgd") || ($species eq "scc")){
	$species='yeast';
	$stats_dir="$DATADIR/gostats/yeast";
	$prob_file=$stats_dir . "/yeast.prob";
	$synfile=$opts{s}||"$DATADIR/yeast.map";
	$gaf_annotations_file=$opts{f} || "$DATADIR/gene_association.yeast";
	$suffix='YEAST';
	# $low_annot_prob= 0.7946262 unless $opts{p};
	# $low_inter_prob= 0.851656 unless $opts{p};    


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
    # if($check_interactome_probs){
    # 	$prob=$low_inter_prob;
    # }
    # else{
    # 	$prob=$low_annot_prob;
    # }

    
}
###############################################################
sub load_network {
    print STDERR "Loading Network $network_file...\n" if $verbose;
    my $A = IO::File->new("< $network_file")|| die("Cannot open network $network_file : $!\n");
    while(<$A>){
	next if /^\d+$/;
	next if /^!/;
	next if /^\s*$/;
	my ($bait,$target)=split(/\s+/,$_);
	$bait=&get_name($bait,"load_network");
	$target=&get_name($target,"load_network");
	## %proteins will hold all the info on all the 
	## network's proteins
	$proteins{$bait}{INTERACTORS}{$target}++;
	$proteins{$target}{INTERACTORS}{$bait}++;
    }
    close($A);
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
	else{die("Could not open $synfile\n" );	}	
 	$synonyms{LOADED}=1;
    }
    my $hname;
    $name =~ /_$suffix/i ? ($hname=$name) : ($hname=$name . "_" . $suffix);
    if(defined($synonyms{NAME}{$name})){
	$name=$synonyms{NAME}{$name} ;
    }
    elsif(defined($synonyms{NAME}{$hname})){
	$name=$synonyms{NAME}{$hname};
    }
    elsif ($name=~/(.+)_$species/i && defined($synonyms{$1})){die("333 : $name\n");
							      $name=$synonyms{$1};
    }
    if ($very_verbose){&debug("Called by : $called_by for $old_name ($synfile), returned $name"); }
    return $name;
}
############################################################
sub parse_gaf_file{    
    print STDERR "Parsing GAF file $gaf_annotations_file...\n" if $verbose;
    open(A,"$gaf_annotations_file")|| die("Cannot open $gaf_annotations_file:$!\n");
    while(<A>){
	next if /^!/;
	chomp;
	my @tmparray=split(/\t/,$_);
	next unless $tmparray[8] eq $subonto;
	next if $tmparray[3] !~/^$/;
	#########################################
	# If we only want "good" evidence codes #
	#########################################
	if ($exp_codes) {
	    next unless defined($good_codes{$tmparray[6]});
	}
	my $name=&get_name($tmparray[1],"parse_gaf_file");
	## If this protein is in the graph
	if(exists($proteins{$name})){
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
    print STDERR "Loading annotations $annotations_file...\n" if $verbose;
    my ($class,$prot);
    my $file_type=0;
    open(A,"$annotations_file")|| die("Cannot open $_[0]:$!\n");
    my $k=0;
    my (@GOs, @bGOs);
    while(<A>){
	if($.==1){
	    if(/^\[CLASS/ || /^\#\s*Inter/){$file_type=1}
	}
	## If this is an annotated classes file
	if($file_type==1){
	    if(/^\[CLASS:\s*(\d+)/){$class=$1}
	    if(/^CA/){
		@GOs=();
		@bGOs=();
		## Benoit decided to change the effin format again, remove the quotes he added
		s/[\[\]]//g;
		s/\',/ /g;
		s/\'//g;
		/^CA\s+(.+?)$/;
		my $kk=$1;
		my @ko=split(/\s+/,$kk);
		
		foreach my $term(@ko){
                    my $T=terms_to_GOs($term,0)||die("$term could not be matched\n");
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
	    	    my $T=&terms_to_GOs($term,0)||die("Unknown term $term\n");
	    	    push @bGOs,$T unless $T eq 'BAD';
	    	}
	    }
	    if(/^PN\s*(.+)/){
		my $kk=$1;
		my @a=split(/,\s*/,$kk);
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
			#print "BBBBB @bGOs :$bGOs[0]:" . scalar(@bGOs) . "\n";
			unless(exists($proteins{$protein}{GOs}{$bGOs[$n]})){ 
			    $proteins{$protein}{GOs}{$bGOs[$n]}="UNC";
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
sub do_test{
    my $result=0;
    my $tests=1;
    if($check_both){
	my ($Pi,$Pa)=@_;
	$tests=$_[2] if $_[2];
	if($either){
	    $result=1 if $Pa*$tests<=$prob or $Pi*$tests<=$prob;
	    return($result,$Pi*$tests,$Pa*$tests) if($_[2]);
	}
        else{
	    $result=1 if $Pa*$tests<=$prob and $Pi*$tests<=$prob;
	    return($result,$Pi*$tests,$Pa*$tests) if($_[2]);
	}
    }
    else{
	my $P=$_[0];
	$tests=$_[1] if $_[1];
	$result=1 if $P*$tests<=$prob;
	return($result,$P*$tests) if $_[1];
    }
    return($result);
}
############################################################
sub save_cand{
    my ($bait, $cls1, $cls2, $go1, $go2, $p1, $p2, $P)=@_;
    my $card1=scalar(keys(%{$classes{$cls1}{PROTS}}));
    my $card2=scalar(keys(%{$classes{$cls2}{PROTS}}));
    my $cc="$bait\t$cls1 $card1 $go1 $p1\t$cls2 $card2 $go2 $p2\t$P";
    return($cc);
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
	  next target if defined($hash{$target}{GOs}{$mcgo});
	  next target if $mcgo eq 'GO:0008150';
	  foreach my $tgo (keys(%{$hash{$target}{GOs}})){
	      next if $tgo eq 'GO:0008150';
	      if ($just_count) {
		  my $pair=join("_", sort  {$b lt $a} ($mcgo,$tgo));
		  $go_counts{$pair}=1;
		  next;
	      }
	      my ($test,$Pa,$Pi,$P);
	      ## If we are checking BOTH annot and inter probs
	      if($check_both){
		  ($Pi,$Pa)=(&interactome_probability($mcgo,$tgo),&gopair_probability($mcgo,$tgo));
		  ## Check if the probbilities are below the threshold
		  $test=&do_test($Pi,$Pa);
	      }
	      ## If we are NOT checking BOTH annot and inter probs
	      else{
		  $check_interactome_probs ?   
		      ($P=&interactome_probability($mcgo,$tgo)) : 
		      ($P=&gopair_probability($mcgo,$tgo));	
		  $test=&do_test($P);
	      }
	      ## skip target if any of its gos are related to the mcgo
	      if($test==0){
		  $not_moonlighting{$bait}{$target}=1;
		  next target;
	      }
	      push @{$gos{$bait}{$target}},  $mcgo . "_" . $tgo ;
	      foreach my $bgo (@{$hash{GOs}}){
		  next if $bgo eq 'GO:0008150';
		  if($check_both){
		      ($Pi,$Pa)=(&interactome_probability(),&gopair_probability());
		      ## Check if the probbilities are below the threshold
		      $test=&do_test($Pi,$Pa);
		  }
		  ## If we are NOT checking BOTH annot and inter probs
		  else{
		      $check_interactome_probs ?   
			  ($P=&interactome_probability()) : 
			  ($P=&gopair_probability());	
		      $test=&do_test($P);
		  }
		  ## skip target if ANY of its gos 
		  ## are close to any bgo
		  if($test==0){
		      $not_moonlighting{$bait}{$target}=1;
		      next target;
		  }
		  else{
		      $moonlighting{$bait}{$target}="$bait ($bgo)  \t:\t$target ($tgo)  \t$Pi\t$Pi";
		  }
		   
	      } ## end foreach my bgo
	  } ## end foreach my tgo
      } #end foreach target
    } ## end if perc >=50
}
############################################################
## Check whether two classes are connected by >1protein 
## (for mode B)
sub class_connected{
    my ($cand, $class1, $class2)=@_;
    my @prots1= keys(%{$classes{$class1}{PROTS}}) ;
    my @prots2= keys(%{$classes{$class2}{PROTS}}) ;
    my $con=0;
 p1:foreach my $prot1 (@prots1){
	return($con) if $con>0;	
    p2:foreach my $prot2 (@prots2){
	    return($con) if $con>0;	
	    if(defined($proteins{$prot1}{INTERACTORS}{$prot2})){
		$con++;
	    }
	}
    }
return($con);	
}
############################################################
sub interactome_probability{
    my $go1=shift;
    my $go2=shift;
    ## Count the number of tests performed unless we are at the end and
    ## are re-checking for multiple testing correction.
    unless($_[2] && $_[2] eq 'nocount'){$TESTS++;  $TESTS2++; }

    ############################
    # Declare database details #
    ############################
    my $host=$db_host;
    my $port = $db_port;
    if ($exp_codes) {
	$database="pogo_exp";
    }
    my $user="root";
    my $pw="yenapas";
    my $table="inter_" . $species;    
    ###############################################
    # If this is the first run, connect to the DB #
    ###############################################
    if (scalar keys (%inter_prob)==0) {
	$dsn = "DBI:mysql:database=$database;host=$host;port=$port;"; 
	##Select database 
	$dbh=DBI->connect($dsn,$user,$pw);
	debug("DB: $host\t$port\t$database\tTable: $table");
    }
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
	my $NN = $network_file . ".inter.prob";
	if ((-e "$NN") && ($Iprob_file_read==0)){
	    open(A,"$NN")||die("Cannot open $NN : $!\n");
	    print STDERR  "\nReading $NN..." if $verbose;
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
	$added_to_inter_prob_file=1 if $read_network_prob_file;
	## Define mysql query
	my $query="SELECT P_low from $table where gopair='$gopair'";  
	debug("QQ : $query;");
	## Run
	my $result=$dbh->prepare("$query");                        
	## Run query
	$result->execute;  
	my $ref = $result->fetchrow_hashref();
	## If a gopair has no prob value, most likely one of its
	## constituent terms does not annotate any proteins DIRECTLY.
	## We can therefore safely ignore the pair by setting P=1.
	$inter_prob{$gopair}=$ref->{'P_low'};
	$inter_prob{$gopair}= 1 unless defined($ref->{'P_low'});
	debug("QQ res : $inter_prob{$gopair}");
    }
    
    defined($inter_prob{$gopair}) ? 
	return($inter_prob{$gopair}):
	return(1);
}
############################################################
sub gopair_probability{
    my $go1=shift;
    my $go2=shift;
    ## Count the number of tests performed unless we are at the end and
    ## are re-checking for multiple testing correction.
    unless($_[0] && $_[0] eq 'nocount'){$TESTS++; $TESTS1++; }
    my $source=shift;
    my $bait=shift;
    my $target=shift;
    my $low;
    my @b = sort {$b lt $a} ($go1,$go2);
    my @c = sort {$a lt $b} ($go1,$go2);
    my $bpair=join("_",@c);
    my $gopair=join("_",@b);
    die("$go1,$go2,$source,$bait,$target\n") if $gopair =~/0008150/;

    ############################
    # Declare database details #
    ############################
    my $host=$db_host;
    if ($exp_codes) {
	$database="pogo_exp";
    }
    my $user="root";
    my $pw="yenapas";
    my $table="$species";    
    my $port=$db_port;

    ###############################################
    # If this is the first run, connect to the DB #
    ###############################################
    if (scalar keys (%prob)==0) {
	$dsn = "DBI:mysql:database=$database;host=$host;port=$port;"; 
	##Select database 
	$dbh=DBI->connect($dsn,$user,$pw);
	debug("DB: $host\t$port\t$database\tTable: $table");
    }

    if (defined($prob{$gopair})){
	return($prob{$gopair}); 
    }
    if ($go1 eq $go2){  
	$prob{$gopair}=1;
	return($prob{$gopair});
    }
    die("bad gopair $gopair: $bpair\n" ) if defined($prob{$bpair});
    ## If one term is an ancestor of the other prob=1;
    if((defined($offspring{$go1}{$go2}))||
       (defined($offspring{$go2}{$go1}))){
	$prob{$gopair}=1;
	return($prob{$gopair});	
    }
    if($read_network_prob_file){
	my $NN = $network_file . ".prob";
	if ( (-e "$NN") && ($Aprob_file_read==0) ){
	    open(A,"$NN")|| die "Could not open prob file $NN . $!\n";
	    print STDERR "\nReading $NN..." if $verbose;    
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
	#print STDERR "Fetching $gopair from DB : ";
	$added_to_prob_file=1;
	## Define mysql query
	my $query="SELECT P_low from $table where gopair='$gopair'";  
	debug("QQ :$query;\n");
	## Run
	my $result=$dbh->prepare("$query");                        
	## Run query
	$result->execute;                                          
	my $ref = $result->fetchrow_hashref();
	## If a gopair has no prob value, most likely one of its
	## constituent terms does not annotate any proteins DIRECTLY.
	## We can therefore safely ignore the pair by setting P=1.
	$prob{$gopair}=$ref->{'P_low'};
	$prob{$gopair}= 1 unless defined($ref->{'P_low'});
	debug("QQ res : $prob{$gopair}");
	#print STDERR "$query;\n" unless defined($prob{$gopair});
    }


    # unless(defined($prob{$gopair})){
    # 	## Since I am only giving it the under represented ones,
    # 	## any probs missing will be over represented so we can 
    # 	## safely ignore them.
    # 	$prob{$gopair}=1;
    # 	push @net_probs, "$gopair\t$prob{$gopair}\n" if $read_network_prob_file;  
    # 	die("cc $gopair: $bpair\n" ) if defined($prob{$bpair});

    # 	return($prob{$gopair});
    # }
    ## @net_probs will be used to create a prob file for this network
    push @net_probs, "$gopair\t" . $prob{$gopair} . "\n";   
    die("dd $gopair: $bpair\n" ) if defined($prob{$bpair});

    return($prob{$gopair});
}


############################################################
# sub debug{
#     if ($debug)
#     {
# 	print STDERR "@_\n";
#     }
# }
############################################################
sub terms_to_GOs{
    my $term=shift;
    $term="biological_process" if $term eq 'biological_process_unknown';
    
    my $mode=shift; ## 0 will return GO:xxx, 1 will return term name
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
    $mode==0 ? 
	return($terms_to_GOs{TERMS}{$term}) :
	return($terms_to_GOs{GOs}{$term}) ;
    
}
############################################################
sub term_precision{
    my $go=shift;
    if($have_already_read_precision_file==0){
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
	print STDERR "MISSING $go : term_precision.pl -n $go -a $go_terms_file $geneology_file | gawk '{print \$NF}'\n";
	my $a=`term_precision.pl -n $go -a $go_terms_file $geneology_file | gawk '{print \$NF}'`;
	$precision{$go}=$a;
	return(nearest(.0001,$precision{$go}));

    }
}
############################################################
sub count_tests{
    # my $TESTS=0;
    # my $mode=shift;
    # switch($mode) {
    #     case ['a','dis2class'] {
    # 	    foreach my $bait (keys(%proteins)){
    # 	      class:foreach my $class (keys(%{$proteins{$bait}{CLASS}})){
    # 		  my @uniqu = grep { ! $seen_annots{$_} ++ } @{$classes{$class}{ANNOT}};
    # 		c_anot:foreach my $class_annot (@uniqu){
    # 		    next c_anot if $is_cand==0;
    # 		    next c_anot if $class_annot eq 'GO:0008150';  
    # 		  bgo:foreach my $bgo (keys(%{$proteins{$bait}{GOs}})){
    # 		      next bgo if $proteins{$bait}{GOs}{$bgo} eq "INF";
    # 		      next bgo if $proteins{$bait}{GOs}{$bgo} eq "UNC";
    # 		      $TESTS++;
    # 		  }
    # 		}
    # 	      }
    # 	    }
    # 	} ## end a
    # } ## end switch
    return($TESTS);
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
	my @l=split(/\t/,$P);
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
#    open(HELP, "| less") ;
#    print HELP <<EndOfHelp;
    print STDERR <<EndOfHelp;

USAGE:  
    $name [options] <ANNOTATIONS FILE> <NETWORK FILE>

	$name will take a class annotation file and a network and
	look for interacting proteins whose GO annotations are very dissimilar.
	It will return a list of proteins with interaction probabilities 
	below a given threshold. It also requires certain files that are
	described at the end of this help.

COMMAND-LINE OPTIONS:
    -A : Database name, def: pogo.
    -a : Annotated class file
    -b : Check BOTH interactiome AND annotation probs (logical OR).
    -B : Check BOTH interactiome AND annotation probs (logical AND).
    -c : Color output
    -C : Ignore ancestor annotations, keep only DIRECT annotations for each protein.
    -d : Debugging output, VERY verbose
    -D : Data directory, default is .\/data.
    -f : Gene Ontology GAF format annotation file (def .\/gene_association.goa_human)
    -h : Print this help and exit
    -g : Only return Gold Standard candidates 							 
    -G : geneology file, output of OntoGenealogy.pl (def .\/<DATA DIRECTORY>\/biological_process.genealogy)
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
           I|Inter     : like 'i' but will only print those candidates that are found at the
                         intersecion of two classes that are ONLY annotated to disimilar GOs.
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
    -l : Do not correct for multiple testing.
    -L : File containing the candidates we want. Any proteins not in this list will be ignored. The 
         file can either be a list (one name per line) or the output of a previous $name run. In the 
         latter case, the multiple testing correction will use the number of tests from the header of the file.
    -N : Read (and update) the network probability file (NETWORK_NAME.prob)
    -o : Gene Ontology:
           P : biological process (default)
	   C : cellular compartment
	   F : biological function
    -p : minimum e-value threshold (def= 0.05)
    -P : Database port for connecting to PoGO (def: 3306).
    -H : Host name/ip for the PoGO database (def: 10.1.3.30). Setting "-H l" sets the host 
         to 127.0.0.1 and port to 3307. Useful when tunneling, eg:
         ssh -fN -p 24222 cchapple\@139.124.66.43 -L 3307:10.1.3.30:3306
    -s : Synonyms file (def .\/<DATA DIRECTORY>\/<SPECIES>.map)
    -S : Species, (def: HUMAN). 
    -t : Number of tests for multiple testing correction. By default, this is the actual number of tests performed
         during the current run. Use this option to set it artificially. 
    -T : GO terms file (def .\/<DATA DIRECTORY>\/GO.terms_ids_obs)
    -v : Verbose output
    -x : Simply count the number of GO pairs tested for this mode. Normal heuristics (skipping to
         the next candidate atthe first GO pair that fails to pass the threshold) are ignored, therefore
         this is the actual number of tests performed for this mode.
         

     
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
#close(HELP);

    print STDERR "\n***** @_ *****\n\n" if @_;
    exit(1);

}
