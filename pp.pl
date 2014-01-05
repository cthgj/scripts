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
			next c3;
		    }
		    my $prec1=&term_precision($class1_annotation);
		  c4:foreach my $class2_annotation (@{$classes{$cls[$k]}{ANNOT}}){ ## annot class 2	  
			 # print STDERR "$n,$k : $class1_annotation : $class2_annotation\n";
		  next c4 if $cls[$k] eq $cls[$n];
		  if($class2_annotation=~/0008150/){
		      next c4;
		  }
		  ############################################################################
                  # If we are just counting, the total number of tests, not skipping 	     #
		  # anything (-x)							     #
                  ############################################################################
		      my $pair=join("_", sort  {$b lt $a} ($class1_annotation,$class2_annotation));
		      $gos_seen++;
		      $go_counts{$pair}=1;
		      next c4;

	      } ## end c4
	    } ## end c3
	  } ## end c2
	} ## end c1
      } ## end foreach bait
    } ## end case i
} ## end switch;

############################################################################
########################### MAIN PROGRAM END ###############################
############################################################################

########################### Print Results ########################

#########################################
# If we are just counting the GO pairs  #
#########################################
    print STDERR "\n" if $verbose;
 print "$_\n" for keys(%go_counts);

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
