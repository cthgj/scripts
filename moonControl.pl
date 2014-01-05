#!/usr/bin/perl -w
use strict;
use Getopt::Std;
use IO::File;
use Switch;
use Term::ANSIColor; 

my %read;
my (%terms_to_GOs,%not_moonlighting,%seen,%moonlighting,%ancestors,%prob,%inter_prob,%offspring,%opts,%proteins,%synonyms,%classes,%precision);
getopts('iACcdRhIgNntVvr:M:m:a:T:D:S:p:P:G:s:f:',\%opts) || do { print "Invalid option, try 'moonGO.pl -h' for more information\n"; exit(1); };
my (@net_probs);
my $very_verbose=$opts{V}||undef;
my $prob_file_read=0; ## is 1 if the prob file for this network has been read, gopair_probability

my $have_already_read_terms_file=0;
my $have_already_read_precision_file=0;
my $have_already_read_probs_file=0;
&usage() if $opts{h};
my $DATADIR=$opts{D}||"./data";
my $stats_dir="$DATADIR/gostats";
my $precision_file="$DATADIR/precision.human";
my $inter_stats_dir;
my $no_anc=$opts{C}||undef;
my $only_gold=$opts{g}||undef;
my $geneology_file=$opts{G}||"$DATADIR/biological_process.genealogy";
my $go_terms_file=$opts{T}||"$DATADIR/GO.terms_alt_ids";
my $gold_file = "$DATADIR/gold_simple.txt";
my %gold; ## information on the gold standard proteins
my $debug=$opts{d}||undef;
my $help=$opts{h}||undef;
my $dont_check_ancestors=$opts{A}||undef;
my $min_class_cardinality=$opts{M}||0;
my $network_file=$ARGV[0]|| die("Need a network file\n");
my $network_prob_file; ## Network specific probabilities
my $low_annot_prob=$opts{p}||0.9393275;#1.193455e-07;
my $low_inter_prob=$opts{p}||0.02463311;#1.193455e-07;
my $check_interactome_probs=$opts{i}||undef;
my $subonto=$opts{o}||"P";
my $random=$opts{r}||undef;
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
    ## Bin probability ranges and count how many interactions involve 
    ## each bin
    case ['bin']{
	my $num=0;
	my (%data, %seen);
	my @bins=(0,0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.10, 0.11, 0.12, 0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19, 0.20, 0.21, 0.22, 0.23, 0.24, 0.25, 0.26, 0.27, 0.28, 0.29, 0.30, 0.31, 0.32, 0.33, 0.34, 0.35, 0.36, 0.37, 0.38, 0.39, 0.40, 0.41, 0.42, 0.43, 0.44, 0.45, 0.46, 0.47, 0.48, 0.49, 0.50, 0.51, 0.52, 0.53, 0.54, 0.55, 0.56, 0.57, 0.58, 0.59, 0.60, 0.61, 0.62, 0.63, 0.64, 0.65, 0.66, 0.67, 0.68, 0.69, 0.70, 0.71, 0.72, 0.73, 0.74, 0.75, 0.76, 0.77, 0.78, 0.79, 0.80, 0.81, 0.82, 0.83, 0.84, 0.85, 0.86, 0.87, 0.88, 0.89, 0.90, 0.91, 0.92, 0.93, 0.94, 0.95, 0.96, 0.97, 0.98, 0.99, 1.00,1.01);
      bait:foreach my $bait (keys(%proteins)){
	  $c++;
	  printf STDERR ("$c of $protnum\r") if $verbose;
	target:foreach my $target (keys(%{$proteins{$bait}{INTERACTORS}})){
	    my @pp=sort {$b lt $a} ($bait,$target);
	    my $pair=join("_",@pp);
	    next target if defined($seen{$pair});
	    $seen{$pair}++;
	  bgo:foreach my $bgo (keys(%{$proteins{$bait}{GOs}})){
		next bgo if $bgo eq 'GO:0008150';
		tgo:foreach my $tgo (keys(%{$proteins{$target}{GOs}})){
		    next tgo if $tgo eq 'GO:0008150';
		    my @b = sort {$b lt $a} ($tgo,$bgo);
		    my $gopair=join("_",@b);
		    next tgo if defined($seen{$gopair});
		    $seen{$gopair}++;
		    # my @b = sort {$b lt $a} ($tgo,$bgo);
		    # my $gopair=join("_",@b);
		    # $data{$gopair}++;
		    my $P=&gopair_probability_whole_file($bgo,$tgo);
		    i:for(my $i=0; $i<$#bins; $i++){
		    	if($P<$bins[$i+1] and $P>=$bins[$i]){
		    	    $data{$bins[$i]}++;
		    	    last i;
		    	}
		    }
		}
	    }
	}
      }
	map{print "$_\t$data{$_}\n"}keys(%data);
    }

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
} ## end switch;


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
	$network_prob_file = $network_file . ".inter.prob";
	$prob=$low_inter_prob;
    }
    else{
	$network_prob_file = $network_file . ".prob";
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
	    	my $a=$2*100/$3;
	    	    if ($a>=30){
	    		push @ko, $1;
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
	    ## skip target if it is annotated to $mcgo
	    next target if defined($hash{$target}{GOs}{$mcgo});
    	    foreach my $tgo (keys(%{$hash{$target}{GOs}})){
		next if $tgo eq 'GO:0008150';
		## skip target if any of its gos are related to the mcgo
		my $P;
		$check_interactome_probs ?   
		    ($P=&interactome_probability($mcgo,$tgo)) : 
		    ($P=&gopair_probability($mcgo,$tgo,"gopair_prob"));
		if($P>=$prob){
		    $not_moonlighting{$bait}{$target}=1;
		    next target;
		}
		foreach my $bgo (@{$hash{GOs}}){
		    next if $bgo eq 'GO:0008150';
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
    if (($prob_file_read==1) && (not defined($inter_prob{$gopair}))){
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
	if ((-e "$network_prob_file") && ($prob_file_read==0)){
	    open(A,"$network_prob_file")||die("Cannot open $network_prob_file : $!\n");
	    print STDERR "\nReading $network_prob_file..." if $verbose;
	    while(<A>){
		chomp;
		my @tt=split(/\t/);
		my $tempgo=$tt[0];
		$inter_prob{$tempgo}=$tt[1];
	    }
	    close(A);
	    $prob_file_read=1;
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
	if ( (-e "$network_prob_file") && ($prob_file_read==0) ){
	    open(A,"$network_prob_file");
	    print STDERR "\nReading $network_prob_file..." if $verbose;
	    while(<A>){
		chomp;
		my @tt=split(/\t/);
		my $tempgo=$tt[0];
		$prob{$tempgo}=$tt[1];
	    }
	    close(A);
	    $prob_file_read=1;
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
sub gopair_probability_whole_file{
    my $go1=shift;
    my $go2=shift;
    my @b = sort {$b lt $a} ($go1,$go2);
    my $gopair=join("_",@b);
    if ($go1 eq $go2){  
	$prob{$gopair}=1;
	return($prob{$gopair});
    }
    ## If one term is an ancestor of the other prob=1;
    unless($dont_check_ancestors){
	if((defined($offspring{$go1}{$go2}))||
	   (defined($offspring{$go2}{$go1}))){
	    $prob{$gopair}=1;
	    return($prob{$gopair});	
	}
    }
    ## Set a random prob foreach gopair
    if($opts{R}){
	$prob{$gopair}=rand;
	return($prob{$gopair});	
    }
    

    if($have_already_read_probs_file==0){
	print STDERR "Reading Probs...";
	$random ? open(P,"$random"):
	    open(P,"cat $stats_dir/annotations/*prob |");
	while(<P>){
	    print STDERR "." if $. % 1000000==0;
	    my @tt=split(/\t/);
	    $prob{$tt[0]}=$tt[$#tt];
	}
	$have_already_read_probs_file=1;
	close(P);
	print STDERR "Done...";
    }
    if (defined($prob{$gopair})){
	return($prob{$gopair}); 
    }
    else{
	$prob{$gopair}=1;    
	return($prob{$gopair}); 
    }

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
    push @{$gold{$bait}{MISSED}{$mode}},"$class\t$class_annot\t$bgo\t$P";
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
