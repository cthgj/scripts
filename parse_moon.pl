#!/usr/bin/perl -w
# use lib "/home/cchapple/lib";
# use moon_cand;
##
## SSH TUNNEL : ssh -fN -p 24222 cchapple@139.124.66.43 -L 3307:10.1.3.30:3306

use strict;
use Getopt::Std;
use Switch;
use Term::ANSIColor; 
use Math::Round;
use DBI();
require "MY_SUBS.pl";
usage() unless $ARGV[0];
my (%opts, %exclude,%cc_prob,%terms_to_GOs, %cc_annots,%found,%gold,%found_gos,%precision, %pairs);
my $have_already_read_terms_file=0;
my $have_already_read_CC_file=0;
## Global MySQL variables
my ($dbh,$dsn);

my $tests=1;
getopts('dMvTgnhGVXPS:t:D:L:c:l:x:s:N:w:f:C:m:H:p:r:R:o:B:',\%opts);
&usage() if $opts{h};
my $excel_style=$opts{X}||undef;
my $species=$opts{S}||undef;
my $db_port=$opts{t}||3306;
my $db_host=$opts{H}||"10.1.3.30";
$db_host="127.0.0.1" if $db_host =~ /^l/i;
$db_port=3307 if $db_host eq "127.0.0.1";
my $count_pairs=$opts{P}||undef;
$species=guess_species($ARGV[0]) unless defined($species);
my $make_networks=$opts{n}||undef;
my $very_verbose=$opts{V}||undef;
my $verbose=$opts{v}||undef;
my $color="bold blue";
my $only_gold=$opts{G}||undef;
my $wanted_modes=$opts{m}||'all';
my $count_methods=$opts{M}||undef;
my $min_prob=$opts{p}||1;
my $cc_prob=$opts{l}||0.05;
my $min_card=$opts{c}||1; ## minimum class cardinality for mode 'i'
my $found_by=$opts{f}||undef;
my $count_gos=$opts{g}||undef;
my $BASEDIR=$opts{B}||$ENV{"HOME"} . "/research/new_moon";
my $DATADIR=$opts{D}||"$BASEDIR/data";
my $class_file=$opts{C}||"$DATADIR/$species.BP.clas";
my $cc_annots_file=$opts{L}||"$DATADIR/$species.CC.clas";
my $cc_probs_file="$DATADIR/$species.CC.prob";
my $network_file=$opts{N}||"$DATADIR/../$species.gr";
my $go_terms_file=$opts{T}||"$DATADIR/GO.terms_alt_ids";
my @MODES;
my $least_likely=1 if $opts{p};
my $gold_file=$opts{F}||"$DATADIR/gold_simple.txt";
my $wanted=$opts{w}||undef;
my $excluded=$opts{x}||undef;
my @files=@ARGV;
my %candidates;
my %fcands;
my $min_prec=0;
my $min_prec2=0;
my $min_sum=$opts{s}||0;
my $total_cands=0; ## counter
our $debug=$opts{d}||undef;
my $check_cc=1;
$check_cc=0 unless $opts{l};

if ($opts{g} && $opts{P}) {
    die "Options -g and -P are mutually exclusive.\n";
}
$count_gos=1 if $count_pairs;

if($opts{r}){
    $min_prec=$opts{r};
}
if($opts{R}){
    $min_prec2=$opts{R};
}
#if($opts{r} && $opts{R}){die "-r means that EITHER go's precision must be below the threshold, -R that BOTH must be. The two options cannot be used together\n"}



&check_gold(); 
&get_class_CC(0);
my (%want, %want_gos, %want_pairs);
if($wanted){
    if(-e $wanted){
	open(W,"$wanted")||die("Could not open file $wanted:$!\n");
	while(<W>){
	    chomp;
	    next if /^Name\t/;
	    /^([^\s]+)\s*/;
	    my $a=$1;
	    $want{$a}++;
	    if (/GO:\d+_/) {
		$want_pairs{$a}++;		
	    }
	    elsif (/GO:\d+/) {
		$want_gos{$a}++;		
	    }
	}
    }
    else{
	my @s=split(/,/,$wanted);
	map{
	    $want{$_}++;
	    if (/GO:\d+_/) {
		$want_pairs{$_}++;		
	    }
	    elsif (/GO:\d+/) {
		$want_gos{$_}++;		
	    }
	}@s;
    }
    die("Error, need a list of names, either a file or a comma separated list with option -w\n") unless scalar(keys(%want))>0;
}
close(W);
if($excluded){
    if(-e $excluded){
	open(X,"$excluded")||die("Could not open file $excluded:$!\n");
	while(<X>){
	    chomp;
	    next if /^Name\t/;
	    /^(.+?)\s/;
	    $exclude{$1}++;
	}
    }
    else{
	my @s=split(/,/,$excluded);
	map{$exclude{$_}++;}@s;
    }
    die("Error, need a list of names, either a file or a comma separated list with option -x\n") unless scalar(keys(%exclude))>0;
}
close(X);
##############################################
if($count_gos){count_gos(\@ARGV)}
else{
    foreach my $file (@ARGV){
	$file=~/.+\.([BabiImpPcd])\./;
	my $mode=$1;
	#$species=guess_species($file);
	if($wanted_modes ne 'all'){
	    next unless $wanted_modes=~/$mode/;
	}
	push @MODES,$mode;
	open(A,"$file")||die("Cannot open $file : $!\n");
	while(<A>){
	    next if /^\s*$/;
	    if (/^\#.*TESTS.*?(\d+)/) {
		$tests=$1;
		next;
	    }
	    next if /^\#/;
	    $total_cands++;
	    chomp;
	    my @a=split(/\t/);
	    if ($wanted || $excluded) {
		my $keep=want_or_not($_);
		next unless $keep==1;
	    }

	    my @a1=split(/\s+/,$a[1]);
	    my @a2=split(/\s+/,$a[2]);
#	  #  for(my $i=0;$i<=$#a;$i++){print "$i:$a[$i]\n"}die();
	    ## Only count results that pass the thresholds
	    if ($a[$#a]<=$min_prob &&  ($a1[3]>=$min_prec || $a2[3]>=$min_prec)  &&
		($a1[3]>=$min_prec2 && $a2[3]>=$min_prec2) &&
		($a1[3] + $a2[3] >=$min_sum) && not defined($exclude{$a[0]})) {
		$found{$a[0]}{$mode}=1 if $a[$#a]<=$min_prob;
		$found{$a[2]}{$mode}=1 if $a[$#a]<=$min_prob and $mode eq 'd';
	    }
	    # else {
	    debug("$a[$#a]<=$min_prob &&  $a1[3]>=$min_prec  &&  $a2[3]>=$min_prec $a1[3]>=$min_prec2 && $a2[3]>=$min_prec2 $a1[3] + $a2[3] >=$min_sum");
	    # }
	    ## Check what type of file this is
	    switch($mode){
		# case /[i]/i{
		#     next unless $a[$#a]<=$min_prob;
		#     $fcands{$a[0]}{NUMS}++;
		#     $fcands{$a[0]}{MODE}{$mode}++ ;
		#     $a[2]=~s/\s*$//;
		#     $a[$#a]=~s/\s*$//;
		#     ## I use this to deal with cases where the same cand
		#     ## appears in a file >1 time
		#     $fcands{$a[0]}{$mode}{$fcands{$a[0]}{MODE}{$mode}}{CLASS1}=$a1[0];
		#     $fcands{$a[0]}{$mode}{$fcands{$a[0]}{MODE}{$mode}}{CLASS2}=$a2[0];
		#     ## Class Cardinality
		#     $fcands{$a[0]}{$mode}{$fcands{$a[0]}{MODE}{$mode}}{CARD1}=$a1[3];
		#     $fcands{$a[0]}{$mode}{$fcands{$a[0]}{MODE}{$mode}}{CARD2}=$a2[3];
			
		#     $fcands{$a[0]}{$mode}{$fcands{$a[0]}{MODE}{$mode}}{go1}= $a1[1] ;
		#     $fcands{$a[0]}{$mode}{$fcands{$a[0]}{MODE}{$mode}}{go2}= $a2[1] ;
		#     $found_gos{$a1[1]}++;
		#     $found_gos{$a2[1]}++;
		    
		#     $a[$#a]=~s/\s*:\s*//;
		#     $fcands{$a[0]}{$mode}{$fcands{$a[0]}{MODE}{$mode}}{SCORE}=nearest(0.000001,$a[$#a]);
		#     $fcands{$a[0]}{$mode}{$fcands{$a[0]}{MODE}{$mode}}{SCORE}=$a[$#a] if $fcands{$a[0]}{$mode}{$fcands{$a[0]}{MODE}{$mode}}{SCORE}==0;	
		#     #$a[3]=~/(\d+\.\d{0,3})/;
		#     # my $o=$1;
		#     # my $b="";
		#     # if($a[$#a]=~/\d+(e\-\d+)$/){ $b=$1;}
		#     # ## catch probs like 0.0000000002
		#     # if($a[4]!~/e/ && $b eq ''){
		#     # 	$fcands{$a[0]}{$mode}{$fcands{$a[0]}{MODE}{$mode}}{SCORE}=$a[4];
		#     # }
		#     # else{
		#     # 	$fcands{$a[0]}{$mode}{$fcands{$a[0]}{MODE}{$mode}}{SCORE}=$o . $b;
		#     # }
		# }
		case /[pP]/{
		    next unless $a[$#a]<=$min_prob;
		    $fcands{$a[0]}{NUMS}++;
		    $fcands{$a[0]}{MODE}{$mode}++;
		    $a[1]=~/^(.+)\s\((GO:\d+)/;
		    ## $fcands{$a[0]}{MODE}{$mode} is just a counter
		    $fcands{$a[0]}{$mode}{$fcands{$a[0]}{MODE}{$mode}}{INTER1}=$1;
		    $fcands{$a[0]}{$mode}{$fcands{$a[0]}{MODE}{$mode}}{go1}=$2 ;
		    $found_gos{$2}++;
		    $a[2]=~/^(.+)\s\((GO:\d+)/;
		    $fcands{$a[0]}{$mode}{$fcands{$a[0]}{MODE}{$mode}}{INTER2}=$1;
		    $fcands{$a[0]}{$mode}{$fcands{$a[0]}{MODE}{$mode}}{go1}=$2;
		    $found_gos{$2}++;

		    $a[$#a]=~s/\s*:\s*//;
		    $fcands{$a[0]}{$mode}{$fcands{$a[0]}{MODE}{$mode}}{SCORE}=nearest(0.000001,$a[$#a]);
		    $fcands{$a[0]}{$mode}{$fcands{$a[0]}{MODE}{$mode}}{SCORE}=$a[$#a] if $fcands{$a[0]}{$mode}{$fcands{$a[0]}{MODE}{$mode}}{SCORE}==0;	
		}
		case /[BbiI]/
		{
		 # for(my $i=0;$i<=$#a;$i++){print "$i:$a[$i]\n"}die("$_\n");
		 next unless  $a[$#a]<=$min_prob;
		 $fcands{$a[0]}{NUMS}++;
		 $fcands{$a[0]}{MODE}{$mode}++;
		 # $a[1]=~/(\d+)\s::\s(GO:\d+)/;
		 $fcands{$a[0]}{$mode}{$fcands{$a[0]}{MODE}{$mode}}{CLASS1}=$a1[0];
		 $fcands{$a[0]}{$mode}{$fcands{$a[0]}{MODE}{$mode}}{go1}=$a1[2];
		 $found_gos{$a1[2]}++;
		 #$a[3]=~/(\d+)\s::\s(GO:\d+)/;
		 $fcands{$a[0]}{$mode}{$fcands{$a[0]}{MODE}{$mode}}{CLASS2}=$a2[0];
		 $fcands{$a[0]}{$mode}{$fcands{$a[0]}{MODE}{$mode}}{go2}=$a2[2];
		 $found_gos{$a2[2]}++;
		 $a[$#a]=~s/\s*:\s*//;
		 #$fcands{$a[0]}{$mode}{$fcands{$a[0]}{MODE}{$mode}}{SCORE}=nearest(0.000001,$a[3]);
		 $fcands{$a[0]}{$mode}{$fcands{$a[0]}{MODE}{$mode}}{SCORE}=nearest(0.000001,$a[$#a]);
		 ############## $fcands{$a[0]}{MODE}{$mode} is just a counter ########
		 $fcands{$a[0]}{$mode}{$fcands{$a[0]}{MODE}{$mode}}{SCORE}=$a[$#a] if $fcands{$a[0]}{$mode}{$fcands{$a[0]}{MODE}{$mode}}{SCORE}==0;	
		 $fcands{$a[0]}{$mode}{$fcands{$a[0]}{MODE}{$mode}}{CARD1}=$a1[1];
		 $fcands{$a[0]}{$mode}{$fcands{$a[0]}{MODE}{$mode}}{CARD2}=$a2[1];
		}
		case /[d]/{
		    next unless $a[$#a]<=$min_prob;
		    $fcands{$a[0]}{NUMS}++;
		    $fcands{$a[0]}{MODE}{$mode}++;
		    $fcands{$a[0]}{$mode}{$fcands{$a[0]}{MODE}{$mode}}{SCORE}=$a[$#a];
		    $fcands{$a[0]}{$mode}{$fcands{$a[0]}{MODE}{$mode}}{INTER1}=$a[2];

		    ## for case d both prots are candidates
		    $fcands{$a[2]}{NUMS}++;
		    $fcands{$a[2]}{MODE}{$mode}++;
		    $fcands{$a[2]}{$mode}{$fcands{$a[2]}{MODE}{$mode}}{SCORE}=$a[$#a];
		    $fcands{$a[2]}{$mode}{$fcands{$a[2]}{MODE}{$mode}}{INTER1}=$a[0];
		    
		    $a[$#a]=~s/\s*:\s*//;
		    $fcands{$a[0]}{$mode}{$fcands{$a[0]}{MODE}{$mode}}{SCORE}=nearest(0.000001,$a[$#a]);
		    $fcands{$a[2]}{$mode}{$fcands{$a[0]}{MODE}{$mode}}{SCORE}=nearest(0.000001,$a[$#a]);
		    $fcands{$a[0]}{$mode}{$fcands{$a[0]}{MODE}{$mode}}{SCORE}=$a[$#a] if $fcands{$a[0]}{$mode}{$fcands{$a[0]}{MODE}{$mode}}{SCORE}==0;	
		    $fcands{$a[2]}{$mode}{$fcands{$a[2]}{MODE}{$mode}}{SCORE}=$a[$#a] if $fcands{$a[2]}{$mode}{$fcands{$a[2]}{MODE}{$mode}}{SCORE}==0;

		    $a[1]=~/^(.+)\s/;
		    $fcands{$a[0]}{$mode}{$fcands{$a[0]}{MODE}{$mode}}{go1}= $1 ;
		    $fcands{$a[2]}{$mode}{$fcands{$a[2]}{MODE}{$mode}}{go2}= $1 ;
		    $found_gos{$1}++;
		    $a[3]=~/^(.+)\s/;
		    $fcands{$a[0]}{$mode}{$fcands{$a[0]}{MODE}{$mode}}{go2}= $1 ;
		    $fcands{$a[2]}{$mode}{$fcands{$a[2]}{MODE}{$mode}}{go1}= $1 ;
		    $found_gos{$1}++;
		}
		case /[cm]/{
		    die("Not ready for mode $mode yet\n");
		    next unless $a[$#a]<=$min_prob;
		    my @b=split(/\s+/,$a[0]);
		    $fcands{$b[0]}++;
		    @b=split(/\s+/,$a[2]);
		    $fcands{$b[0]}{NUMS}++;
		    $fcands{$b[0]}{MODE}{$mode}++;
		    $fcands{$b[0]}{$mode}{$fcands{$a[0]}{MODE}{$mode}}{SCORE}=$a[$#a];
		    $a[$#a]=~s/\s*:\s*//;
		    $fcands{$a[0]}{$mode}{$fcands{$a[0]}{MODE}{$mode}}{SCORE}=nearest(0.000001,$a[$#a]);
		    $fcands{$a[0]}{$mode}{$fcands{$a[0]}{MODE}{$mode}}{SCORE}=$a[$#a] if $fcands{$a[0]}{$mode}{$fcands{$a[0]}{MODE}{$mode}}{SCORE}==0;	
		}
		case /a/{
		    next unless $a[$#a]<=$min_prob;
		    my $name=$a[0];
		    $fcands{$name}{NUMS}++;
		    $fcands{$name}{MODE}{$mode}++;
		    $a[1]=~/(GO:\d+)/;
		    $fcands{$name}{$mode}{$fcands{$a[0]}{MODE}{$mode}}{go1}=$1;
		    $found_gos{$1}++;
		    $a[3]=~/(\d+),\s*(GO:\d+)\s*/;
		    $fcands{$name}{$mode}{$fcands{$a[0]}{MODE}{$mode}}{CLASS1}=$1;
		    $fcands{$name}{$mode}{$fcands{$a[0]}{MODE}{$mode}}{go2}=$2;
		    $found_gos{$2}++;
		    $a[$#a]=~s/\s*:\s*//;
		    $fcands{$a[0]}{$mode}{$fcands{$a[0]}{MODE}{$mode}}{SCORE}=nearest(0.000001,$a[$#a]);
		    $fcands{$a[0]}{$mode}{$fcands{$a[0]}{MODE}{$mode}}{SCORE}=$a[$#a] if $fcands{$a[0]}{$mode}{$fcands{$a[0]}{MODE}{$mode}}{SCORE}==0;	
		}
	    } ## end switch
	}## end while(<A>)
   } ## end foreach my $file

    ## get precision values for all GOs found
    &term_precision(0,\%found_gos); 
    my $left=$total_cands;
    ## If we just want those candidates that were found by
    ## at least $found_by methods
    if($found_by){	    
	my $total=0;
	my @modes=sort {uc($a) cmp uc($b)} @MODES;
	$"="\t";
	print STDERR "\t\tTOT\t@modes\n";
	cand:foreach my $cand (keys(%fcands)){
	    $left--;
	    print STDERR "$left of $total_cands\r" if $verbose;
	    if($only_gold){
		next unless defined($gold{$cand}) ;
	    }
	    next cand unless scalar(keys(%{$found{$cand}}))>=$found_by;
	    $total++;

	    my $cases;
	    foreach my $M (@MODES){
		#my @aa=keys(%{$found{$cand}{$M}});
		$cases+=$fcands{$cand}{MODE}{$M} if defined($fcands{$cand}{MODE}{$M}) ;
	    }
	    my $gcand=$cand;
	    $gcand=color("$color").$cand.color("reset") unless $excel_style;
	    my ($mm,$spacer)=(1,"");
	    if(defined($gold{$cand})){
		$mm=15-length("$cand");
		for (my $i=0;$i<=$mm;$i++){
		    $spacer.=" ";
		}
		print "$gcand$spacer$cases\t";
	    }
	    else{
		$mm=15-length("$cand");
		for (my $i=0;$i<=$mm;$i++){
		    $spacer.=" ";
		}
		print "$cand$spacer$cases\t";
	    }

	    my @modes=sort {uc($a) cmp uc($b)} keys(%{$fcands{$cand}{MODE}});
	    foreach my $M (sort {uc($a) cmp uc($b)} @MODES){
		if(exists($found{$cand}{$M})){
		    print "$fcands{$cand}{MODE}{$M}\t";
		}
		else{
		    print "0\t";
		}
	    }#map{print "$_ " }@modes;
	    print "\n";
	    if($make_networks){
		foreach my $mode (@modes){
		    next unless $mode =~/[aBbi]/;
		    next unless defined($found{$cand}{$mode});
		    for (my $i=1; $i<=$fcands{$cand}{MODE}{$mode}; $i++){
			&make_networks($cand,$mode,$fcands{$cand}{$mode}{$i}{CLASS1},$fcands{$cand}{$mode}{$i}{CLASS2});
		    }
		}
	    }
	}
	  # if($only_gold){
	  #     $total=0;
	    #     map {$total++ if defined($gold{$_})}keys(%fcands);
	  # }
	print "TOTAL : $total\n";
#	  map{print "xx $_\n"}keys(%fcands);
    }
    ## Just print what methods found each candidate 
    elsif($count_methods){
	foreach my $cand (keys(%fcands)){
	    $left--;
	    print STDERR "$left of $total_cands\r" if $verbose;
	    if($only_gold){
		next unless defined($gold{$cand}) ;
	    }
	    if($wanted && not defined($want{$cand})){ next;}
	    if($excluded &&  defined($exclude{$cand})){ next;}
	    my $gcand=$cand;
	    unless ($excel_style){
		defined($gold{$cand}) ? 
		    ($gcand=color("$color").$cand.color("reset")) : 
			($gcand=$cand);
	    }
	    print "$gcand\t" . scalar(keys(%{$fcands{$cand}{MODE}})) . "\t";
	    foreach my $M (sort {uc($a) cmp uc($b)} @MODES){
		defined($found{$cand}{$M}) ? print "$M\t" : print "\t";
	    }
#	map{print "$_ " }keys(%{$fcands{$cand}{MODE}});
	    print "\n";
	}
    }
     ## If we want the least likely candidates for a given method(s)
    else{
	if ($excel_style) {
	    print "Name\tMode\tClass1 (card.)\tClass1: BP (GO1, prec)\tClass1:CC\tClass2 (card.)\tClass2: BP (GO2, prec)\tClass2:CC\tCCprob.\tBP prob.\n"	}
	else {
	    $very_verbose ? 
		print "##Mode\tClass1 (card.)\tClass1: BP (GO1, prec)\tClass1:CC\tClass2 (card.)\tClass2: BP (GO2, prec)\tClass2:CC\tCCprob.\tBP prob.\n":
		print "##Mode\tClass1:BP Annotation,Class1:CC Annotation\tClass2:BP_Annotation,Class2:CC_Annotation\tCC_probability\tBP_probability\t(Class1: GO, precision)\t(Class2: GO, precision)\n";
	}
	my $total=0;
	my $cases=0; ## How many cases (not cands) found
	my @modes;
	if($wanted_modes ne 'all'){
	    @modes=split(//,$wanted_modes);
	}
	else{@modes=@MODES;}
	cand:foreach my $cand (keys(%fcands)){
	    my $kk=0; ## used to count later, see $total
	    $left--;
	    print STDERR "$left of $total_cands\r" if $verbose;
	    if($only_gold){
		next unless defined($gold{$cand}) ;
	    }
	    my $gcand=$cand;
	    
	   unless ($excel_style){
	       defined($gold{$cand}) ? 
		   ($gcand=color("$color").$cand.color("reset")) : 
		       ($gcand=color("bold").$cand.color("reset"));
	   }
	    my $printed=0;
	    foreach my $m ((keys(%{$fcands{$cand}{MODE}}))){
		## for each instance of this cand in this file
		my @to_be_printed;
	      ii:for (my $i=1; $i<=$fcands{$cand}{MODE}{$m}; $i++){
		      my $ccprob=2;
		      next unless $fcands{$cand}{$m}{$i}{SCORE}<=$min_prob;
		      # ## Skip if we do NOT want this GO
		      # if($excluded &&  
		      # 	 defined($exclude{$fcands{$cand}{$m}{$i}{go1}})||
		      # 	 defined($exclude{$fcands{$cand}{$m}{$i}{go2}})	){ next;}
		      my $prec1=&term_precision(1,$fcands{$cand}{$m}{$i}{go1})||1;
		      my $prec2=&term_precision(1,$fcands{$cand}{$m}{$i}{go2})||1;
		      my (@class1_CC_annots,@class2_CC_annots);
		      unless($m eq 'P' or $m eq 'd'){
			  @class1_CC_annots=@{&get_class_CC($fcands{$cand}{$m}{$i}{CLASS1})};
			  @class2_CC_annots=@{&get_class_CC($fcands{$cand}{$m}{$i}{CLASS2})} unless $m eq 'a';
		      }
		      ## skip classes with fewer than desired members
		      if($m eq 'i' || $m eq 'I' || $m eq 'B'){
			  next ii unless $fcands{$cand}{$m}{$i}{CARD1}>=$min_card;
			  next ii unless $fcands{$cand}{$m}{$i}{CARD2}>=$min_card;
		      }
		      for (my $k=0;$k<=$#class1_CC_annots;$k++){
			  for (my $kk=0;$kk<=$#class2_CC_annots;$kk++){
			      my $cc=&gopair_prob($class1_CC_annots[$k],$class2_CC_annots[$kk], $tests);
			      if($check_cc and not ($m eq 'a' || $m eq 'P' || $m eq 'd')){		      
				next ii unless $cc<=$cc_prob;
			    }
			      $ccprob=$cc if $cc<=$ccprob;
			  }
		      }
		      $ccprob="N/A" if $ccprob==2;
		    # if($check_cc and not ($m eq 'a' || $m eq 'P' || $m eq 'd')){		      
		    # 	for (my $k=0;$k<=$#class1_CC_annots;$k++){
		    # 	    for (my $kk=0;$kk<=$#class2_CC_annots;$kk++){
		    # 		my $cc=&gopair_prob($class1_CC_annots[$k],$class2_CC_annots[$kk], $tests);
		    # 		die();
		    # 		next ii unless $cc<=$cc_prob;
		    # 		$ccprob=$cc if $cc<=$ccprob;
		    # 		$ccprob="N/A" if $ccprob==2;
		    # 	    }
		    # 	}
		    # }
		    next unless $prec1+$prec2>=$min_sum;
		    if($opts{R}){			
			next unless $prec1>=$min_prec2;
			next unless $prec2>=$min_prec2;
		    }
		      if ($opts{r}) {
			  unless( $prec1>=$min_prec ||
				  $prec2>=$min_prec){
			      next;
			  }
		      }
		      $kk=1;
		      if ($printed==0) {
			  push @to_be_printed, "$gcand\n" unless $excel_style;
			  $printed=1;
		      }
		      $cases++;
		      push @to_be_printed, &print_results($cand,$m,$i,$prec1,$prec2,$ccprob,'bob');
		     
		      if (($opts{n}) && ($m eq 'i' || $m eq 'I' || $m eq 'B' || $m eq 'P' || $m eq 'b')) { 
			  #for (my $i=1; $i<=$fcands{$cand}{MODE}{$m}; $i++){
			  &make_networks($cand,$m,$fcands{$cand}{$m}{$i}{CLASS1},$fcands{$cand}{$m}{$i}{CLASS2}) if defined($found{$cand}{$m});
			  #}
		      } elsif (($opts{n}) && ($m eq 'a')) {
			  #for (my $i=1; $i<=$fcands{$cand}{MODE}{$m}; $i++){
			  &make_networks($cand,$m,$fcands{$cand}{$m}{$i}{CLASS1}) if defined($found{$cand}{$m});
			  #}
		      }

		  }		# end for $i
		######################################
                # Print the formatted results	     #
                ######################################
		map{
		    $excel_style ? 
			print "$gcand\t$_" : 
			    print ; 
		}@to_be_printed;
	    }			## end foreach $m
	    $total++ if $kk==1;
	}			# end foreach cand
	print STDERR "\n" if $verbose;
	print STDERR "TOTAL : $cases cases for $total candidates \n";
    }
}
############################################################
sub want_or_not{
 #   my @line=@_;
    my @gos=(/(GO:\d+)/g);
    my $pair=join("_",sort {$b lt $a} (@gos));
    my @names=(/([^\s]+_\w{5})/g);
    my $wnt=0;
    my $ex=0;
    map{$ex++ if defined ($exclude{$_})} (@gos, $pair, @names);
    map{$wnt++ if defined ($want{$_})}(@gos, $pair, @names);

    if ($ex>0) {return(0)}
    if ($wnt>0) {return(1)}
    # if ($ex>0) {return(0)}
    # if ($ex>0) {return(0)}
}

############################################################
sub make_networks{
    my $cand=shift;
    my $mode=shift;
    if(($mode eq 'i')|| ($mode eq 'I') ||($mode eq 'B') || ($mode eq 'b') || ($mode eq 'a')){
	my $class1=shift;
	my $class2;
	$mode eq 'a' ? ($class2='none') : ($class2=shift);
	my %names;
	open(A,"moon_get_class_prots.pl $class1 $class_file |");
	while(<A>){
	    chomp;
	    next if $cand eq $_;
	    $names{$_}{$class1}++;
	}
	close(A);
	unless ($mode eq 'a'){
	    open(B,"moon_get_class_prots.pl $class2 $class_file |")||die "Cannot run moon_get_class_prots.pl: $!\n";
	    while(<B>){
		chomp;
		next if $cand eq $_;
		$names{$_}{$class2}++;
	    }
	    close(B);
	}
	my $i=$class1 . "_"  . $class2;
	my $attr_file=$cand . "." . $mode . ".$i" . ".class.attr";
        my $attr_file1=$cand . "." . $mode . ".$i" . ".true_class.attr";
        my $attr_file2=$cand . "." . $mode . ".$i" . ".names.attr";
	open(C,"> $attr_file")||die("Could not open $attr_file for writing:$!\n");
	open(C1,"> $attr_file1")||die("Could not open $attr_file1 for writing:$!\n");
	open(C2,"> $attr_file2")||die("Could not open $attr_file2 for writing:$!\n");
	print C "Class (class=String)\n$cand" . "_$cand = 0\n";
	print C1 "True_Class (class=String)\n";
	print C2 "Label\n$cand" . "_$cand = $cand\n";
	print C1 "$cand" . "_$cand = ($class1" . "::$class2)\n" if $mode eq 'i';
	foreach my $name (keys(%names)){
	    my $NetName=$cand . "_" . $name;
	    print C2 "$NetName = $name\n";
	    if (defined($names{$name}{$class1}) && defined($names{$name}{$class2})){
		print C "$NetName = AB\n";
		print C1 "$NetName = ($class1" . "::$class2)\n";
	    }
	    elsif (defined($names{$name}{$class1})){
		print C "$NetName = A\n";
		print C1 "$NetName = ($class1)\n";
	    }
	    elsif(defined($names{$name}{$class2})){
		print C "$NetName = B\n";
		print C1 "$NetName = ($class2)\n";
	    }
	    else{die("Huh!!??\n$NetName,$class1,$class2,$class1,$class2\n");}
	}
	close(C);
	my $subnet_file=$cand . "." .  $mode . ".$i.sif";
	my $subnet_png=$cand . "." .  $mode . ".$i.png";
	open(S,">$subnet_file")||die "Could not open $subnet_file:$!\n";
	open(N,"$network_file")||die("Could not open $network_file:$!\n");
	while(<N>){
	    next if /^\d+$/;
	    next if /^\s+$/;
	    chomp;
	    my @a=split(/\t/);
	    if(((($a[0] eq $cand) || ($a[1] eq $cand)) && (defined($names{$a[0]}) ||defined($names{$a[1]})))||
	       (defined($names{$a[0]}) && defined($names{$a[1]}))){
		print S "$cand" . "_$a[0] pp $cand" . "_$a[1]\n";
	    }
	}
	close(S);
	close(N);
	
	my $cyto_script=$cand . "." . $mode  . ".$i" . ".cyto.txt";
	open(C3,">$cyto_script")||die "Could not open $cyto_script:$!\n";
	my $cur_dir=`pwd`;
	chomp($cur_dir);
	print C3 "session open file=\"$BASEDIR/results/candidates.cys\"\n";
	print C3 "network import file=\"$cur_dir/$subnet_file\"\n";
	print C3 "node import attributes file=\"$cur_dir/$attr_file\"\n";
	print C3 "node import attributes file=\"$cur_dir/$attr_file1\"\n";
	print C3 "node import attributes file=\"$cur_dir/$attr_file2\"\n";
	print C3 "layout force-directed\n";
	print C3 "network view export file=\"$cur_dir/$subnet_png\" zoom=\"5\"\n"; 
	close(C3);
	

    } 
}


############################################################
sub count_gos{
    my @files=@{$_[0]};
    my (%prots,%gos);
    my $total_cands=0;
    my $cases=0;
    foreach my $file (@files){
	$file=~/.+\.([BabIimpPcd])\./;
	my $mode=$1;
	
	open(A,"$file")||die("Cannot open $file : $!\n");
	while(<A>){
	    next if /^\s*$/;
	    next if /^\#/;
	    chomp;
	    /^(.+?)\s/;
	    my $name=$1;
	#    $gos{PROTS}{$name}++;
	    if($mode =~ /[ib]/i){
		my @a=split(/\t/);
		map{s/\s*$//}@a;
		my @a1=split(/\s+/,$a[1]);
		my @a2=split(/\s+/,$a[2]);
		if ($a[$#a]<=$min_prob &&  $a1[3]>=$min_prec  &&  $a2[3]>=$min_prec &&
		    ($a1[3]>=$min_prec2 && $a2[3]>=$min_prec2) && 
		    ($a1[3] + $a2[3] >=$min_sum) && !defined($exclude{$a[0]}) && !defined($exclude{$a1[2]})
		    && !defined($exclude{$a2[2]})) {
		    if ($count_pairs) {
			my $gopair=join("_",sort {$b lt $a} ($a1[2],$a2[2]));
			$pairs{$gopair}++;
		    }
		    $gos{$a1[2]}++;
		    $gos{$a2[2]}++;
		    $prots{$name}=1;
		    ## For many prots, more than one case is found that makes them candidates.
		    ## In a bridge case, for example, a protein can link >1 dissimilar classes
		    ## So, when calculating the percentage, I need the total number of CASES.
		    $cases++;
		}
	    }
	    else{die("Not sure if this works for modes other than i,b,I,B. Current mode is $mode\n");
		my @a=/(GO:\d+)/g;
		$gos{$a[0]}++;
		$gos{$a[1]}++;
	    }
	}
    }
    my @GOs=keys(%gos);
    ## get precision values for all GOs found
    &term_precision(0,\%gos);
    my $total_gos=scalar(@GOs);
    my $total_prots=scalar(keys(%prots));
    my $total_pairs=scalar keys(%pairs);
    my %output;
    if ($count_pairs) {
	print "## GO pair\tCands\tCases\tGOs\tFreq\t%\tPrec1\tPrec2\tTerm1\tTerm2\tP-value\n";
	my @go_pairs=keys(%pairs);
	foreach my $go_pair (@go_pairs){
	    next if $go_pair eq 'TOTAL_PROTS';
	    next if $go_pair eq 'PROTS';
	    my $p=sprintf("%.2f",(($pairs{$go_pair})*100)/($cases));
	    my $koko="$go_pair\t$total_prots\t$cases\t$total_gos\t$pairs{$go_pair}\t$p";
	    my ($go1,$go2)=split(/_/,$go_pair);
	    my $prob=gopair_prob($go1,$go2);
	    $koko.="\t"  . &term_precision(1,$go1) . "\t"  . &term_precision(1,$go2) . 
		"\t" . &terms_to_GOs($go1,1,"PROTS") . "\t" . &terms_to_GOs($go2,1,"P") 
		    . "\t" . $prob;
	    push @{$output{$p}}, $koko;
	}
    }
    else {
	print  "## GO term\tCands\tCases\tGOs\tFreq\t%\tPrec\tTerm name\n";
	foreach my $go (@GOs){
	    next if $go eq 'TOTAL_PROTS';
	    next if $go eq 'PROTS';
	    #	print "$go\t$total_prots\t$cases\t$total_gos\t$gos{$go}\t";
	    ## Each case is associated to 2 GOs, so the number of slots
	    ## a go can take is $cases*2. That is, to be 100% a go has
	    ## to be present $case*2 times.
	    my $p=sprintf("%.2f",(($gos{$go})*100)/($cases*2));
	    my $koko="$go\t$total_prots\t$cases\t$total_gos\t$gos{$go}\t$p";
	    $koko.="\t"  . &term_precision(1,$go) . "\t" . &terms_to_GOs($go,1,"P");
	    push @{$output{$p}}, $koko;
	   
	}
    } 
    my @sorted=sort{$a <=> $b} keys(%output);
    map{
	map{
	    print "$_\n";
	}@{$output{$_}}
    }@sorted;
    $verbose && do {
	$count_pairs ? 
	    print STDERR "$total_pairs GO pairs were found\n" :
		print STDERR "$total_gos GOs were found\n";
    }

}

############################################################
sub terms_to_GOs{
    my $term=shift;
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
	    $terms_to_GOs{TERMS}{$a[$#a-1]}=$a[0];
	    #print "\$terms_to_GOs{TERMS}{$a[$#a-1]}=$a[0];\n";
	    map{
		$terms_to_GOs{GOs}{$_}=$a[$#a-1];
	    }@terms;
	    #}
	    
	}
	close(T);
	$terms_to_GOs{TERMS}{"cellular compartment unknown"}="GO:0005575";
	$terms_to_GOs{TERMS}{"cellular component unknown"}="GO:0005575";
	$terms_to_GOs{GOs}{"GO:0005575"}="cellular component unknown";
	$have_already_read_terms_file=1;
    }
    #print "kk : $term : $terms_to_GOs{TERMS}{$term} : $mode : \$terms_to_GOs{GOs}{$term}\n";
#    &debug("term : $term, id:$terms_to_GOs{TERMS}{$term}, id:$terms_to_GOs{TERMS}{$term} " );
#    print STDERR "term : $term, id:$terms_to_GOs{TERMS}{$term} \n";
    if($mode==0){
	die("1 Could not find $term ") unless defined($terms_to_GOs{TERMS}{$term});
	return($terms_to_GOs{TERMS}{$term}) ;
    }
    else{
	die("2 Could not find $term ($mode)") unless defined($terms_to_GOs{GOs}{$term});
	return($terms_to_GOs{GOs}{$term}) ;
    }
    
}

############################################################

sub print_results{
    my $cand=shift;
    my $m=shift;
    my ($case,$prec1,$prec2,$ccprob)=@_;
    $ccprob="" unless defined($ccprob);
    my $gcand=$cand;
    my @r;
    unless ($excel_style){
	defined($gold{$cand}) ?
	    ($gcand=color("$color").$cand.color("reset")) : 
		($gcand=$cand);
    }
    my $score=$fcands{$cand}{$m}{$case}{SCORE};
    if($score<=1e-7 && $score>0){
	$score=~/(.+)(e.+)/||die("No match : $score\n");
	$score= nearest(0.01,$1) . "$2";
    }  
    if($ccprob<=1e-7 && $ccprob>0){
	$ccprob=~/(.+)(e.+)/||die("No match : $ccprob\n");
	$ccprob= nearest(0.01,$1) . "$2";
    }  
    if($m eq 'P'){
	if (defined($found{$cand}{$m})){    
	    for (my $i=1; $i<=$fcands{$cand}{MODE}{$m}; $i++){
		if($excluded &&  defined($exclude{&terms_to_GOs($fcands{$cand}{$m}{$i}{go1},1,"P")})){ next;}
		if($excluded &&  defined($exclude{&terms_to_GOs($fcands{$cand}{$m}{$i}{go2},1,"P")})){ next;}
		push @r, "\tMODE:P\tInter1: $fcands{$cand}{$m}{$i}{INTER1}\tInter1_GO: \"" . &terms_to_GOs($fcands{$cand}{$m}{$i}{go1},1,"P")  ."\" ($fcands{$cand}{$m}{$i}{go1} $prec1) \tInter2: $fcands{$cand}{$m}{$i}{INTER2}\tInter2_GO: \"" . &terms_to_GOs($fcands{$cand}{$m}{$i}{go2},1,"P")  ."\" ($fcands{$cand}{$m}{$i}{go2} $prec2) \t$fcands{$cand}{$m}{$i}{SCORE}\t$ccprob\n";
	    }
	}
    }    
    else{
	$"=",";
	    my @class1_CC_annots=@{&get_class_CC($fcands{$cand}{$m}{$case}{CLASS1})} unless $m eq 'd';
	    my @class2_CC_annots=@{&get_class_CC($fcands{$cand}{$m}{$case}{CLASS2})} unless $m eq 'a' or $m eq 'd';

	if($m eq 'B' || $m eq 'b'){
	    unless($excluded &&  
		   defined($exclude{&terms_to_GOs($fcands{$cand}{$m}{$case}{go1},1,"P")}) ||
		   defined($exclude{&terms_to_GOs($fcands{$cand}{$m}{$case}{go2},1,"P")})){ 
		if ($very_verbose) {
		    push @r, "MODE:$m\tClass1: $fcands{$cand}{$m}{$case}{CLASS1}\t\"" . &terms_to_GOs($fcands{$cand}{$m}{$case}{go1},1,"P") . "\" ($fcands{$cand}{$m}{$case}{go1}, $prec1)\t@class1_CC_annots\tClass2: $fcands{$cand}{$m}{$case}{CLASS2}\t\"" . &terms_to_GOs($fcands{$cand}{$m}{$case}{go2},1,"P") . "\" ($fcands{$cand}{$m}{$case}{go2}, $prec2)\t@class2_CC_annots\t$ccprob\t$score\n";
		}
		else {
		    push @r, "\t$m: \"" . &terms_to_GOs($fcands{$cand}{$m}{$case}{go1},1,"P") . "\",@class1_CC_annots\t\"" . &terms_to_GOs($fcands{$cand}{$m}{$case}{go2},1,"P") . "\",@class2_CC_annots\t$ccprob\t$score\t($fcands{$cand}{$m}{$case}{go1}, $prec1)\t($fcands{$cand}{$m}{$case}{go2}, $prec2)\n";
		}
	    }
	}
	elsif($m eq 'i' || $m eq 'I'){
	    unless($excluded &&  
		   defined($exclude{&terms_to_GOs($fcands{$cand}{$m}{$case}{go1},1,"P")}) ||
		   defined($exclude{&terms_to_GOs($fcands{$cand}{$m}{$case}{go2},1,"P")})){ 
		if ($very_verbose) {
		    push @r, "MODE:$m\tClass1: $fcands{$cand}{$m}{$case}{CLASS1} ($fcands{$cand}{$m}{$case}{CARD1})\t\"" . &terms_to_GOs($fcands{$cand}{$m}{$case}{go1},1,"P") . "\" ($fcands{$cand}{$m}{$case}{go1}, $prec1)\t@class1_CC_annots\tClass2: $fcands{$cand}{$m}{$case}{CLASS2} ($fcands{$cand}{$m}{$case}{CARD2})\t\"" . &terms_to_GOs($fcands{$cand}{$m}{$case}{go2},1,"P") . "\" ($fcands{$cand}{$m}{$case}{go2}, $prec2)\t@class2_CC_annots\t$ccprob\t$score\n";
		}
		else {
		    push @r, "\t$m: \"" . &terms_to_GOs($fcands{$cand}{$m}{$case}{go1},1,"P") . "\",@class1_CC_annots\t\"" . &terms_to_GOs($fcands{$cand}{$m}{$case}{go2},1,"P") . "\",@class2_CC_annots\t$ccprob\t$score\t($fcands{$cand}{$m}{$case}{go1}, $prec1)\t($fcands{$cand}{$m}{$case}{go2}, $prec2)\n";		    
		}
	    }
	}
	elsif($m eq 'a'){
	    unless($excluded &&  
		   defined($exclude{&terms_to_GOs($fcands{$cand}{$m}{$case}{go1},1,"P")}) ||
		   defined($exclude{&terms_to_GOs($fcands{$cand}{$m}{$case}{go2},1,"P")})){ 
		push @r, "\tMODE:$m\t\"" . &terms_to_GOs($fcands{$cand}{$m}{$case}{go1},1,"P") . "\" ($fcands{$cand}{$m}{$case}{go1}, $prec1)\tClass: $fcands{$cand}{$m}{$case}{CLASS1}\t\"" . &terms_to_GOs($fcands{$cand}{$m}{$case}{go2},1,"P") . "\" ($fcands{$cand}{$m}{$case}{go2}, $prec2)\t\t@class1_CC_annots\t$score\t$ccprob\n";
	    }
	}
	elsif($m eq 'd'){#malaka
		unless($excluded &&  
		   defined($exclude{&terms_to_GOs($fcands{$cand}{$m}{$case}{go1},1,"P")}) ||
		   defined($exclude{&terms_to_GOs($fcands{$cand}{$m}{$case}{go2},2,"P")})){ 
	
		push @r, "\t MODE:$m\tCandGO: \"" . &terms_to_GOs($fcands{$cand}{$m}{$case}{go1},3,"P") . "\" ($fcands{$cand}{$m}{$case}{go1}, $prec1)\t$fcands{$cand}{$m}{$case}{INTER1}\t\"" . &terms_to_GOs($fcands{$cand}{$m}{$case}{go2},4,"P") . "\" ($fcands{$cand}{$m}{$case}{go2}, $prec2)\t$score\t$ccprob\n";
	    }
	}
	else{
	    unless($excluded &&  
		   defined($exclude{&terms_to_GOs($fcands{$cand}{$m}{$case}{go1},1,"P")}) ||
		   defined($exclude{&terms_to_GOs($fcands{$cand}{$m}{$case}{go2},1,"P")})){ 
		push @r, "\t MODE:$m\tCandGO: \"" . &terms_to_GOs($fcands{$cand}{$m}{$case}{go1},1,"P") . "\" ($fcands{$cand}{$m}{$case}{go1}, $prec1)\t@class1_CC_annots\tClass: $fcands{$cand}{$m}{$case}{CLASS1}\tClassGO: \"" . &terms_to_GOs($fcands{$cand}{$m}{$case}{go2},1,"P") . "\" ($fcands{$cand}{$m}{$case}{go2}, $prec2)\t@class2_CC_annots\t$score\t$ccprob\n";
	    }
	}
    }
    return(@r);
}

############################################################

sub check_gold{
    ## Load gold standard
    unless($gold{READ}){
	open(G,"$gold_file")||die("Could not open $gold_file : $!");
	while(<G>){
	    next if $.==1;
	    my @a=split(/\t/);
	    next if $a[2]==0; ## skip the ones we have removed from the GS
	    $gold{$a[0]}{ISO}=$a[1];
	    $gold{$a[0]}{GOOD}=$a[2];
	    $gold{$a[0]}{NUM}=$a[3];
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
    push@{$gold{$bait}{MISSED}{$mode}},"$class\t$class_annot\t$bgo\t$P";
    

}
#####################################################################
sub term_precision{
    my $mode=shift;
    if($mode==0){
	my $hash=shift;
	open(PP,">/tmp/$$")|| die("Could not open precision file /tmp/$$ for writing:$!\n");
	map{print PP "$_\n";}keys(%{$hash});
	close(PP);
	system("term_precision.pl -g $DATADIR/biological_process.genealogy  /tmp/$$ > /tmp/$$.b");
#	print STDERR "term_precision.pl $DATADIR/biological_process.genealogy  /tmp/$$ > /tmp/$$.b\n";
	open(PO,"/tmp/$$.b")|| die "Could not open percision outfile /tmp/$$.b:$!\n";
	while(<PO>){
	    chomp;
	    /(.+)\t(.+)/;
	    $precision{$1}=$2;
	}
	return();
    }
    else{
	my $go=shift;
	die("MISSING : $go ($mode)\n") unless defined($precision{$go});
	return(nearest(.0001,$precision{$go}));
    }
}
#####################################################################
sub get_class_CC{
    my $class=shift;
    if($class==0){ ## if this is the first run, read the file
	my $class_name;
	open(CC,"$cc_annots_file")||die "(sub get_class_CC) Could not open class annotation file $cc_annots_file\n";
	while(<CC>){
	    if(/^\[CLASS:\s*(\d+)/){$class_name=$1}
	    if(/^CA\s+(.+?)$/){
		my $kk=$1;
		$kk=~s/\[\'//;
		$kk=~s/\'\]//;
		$kk=~s/\', \'/ /g;
		my @GOs=();

		my @ko=split(/\s+/,$kk);
		foreach my $term(@ko){
#		    my $T=&terms_to_GOs($term,0,"C");
		    push @{$cc_annots{$class_name}},$term;
		}
	    }
	}
	close(CC);
    }
    else{return($cc_annots{$class});}
}
#####################################################################
sub gopair_prob{
    my $annot1=$_[0];
    my $annot2=$_[1];
    if ($annot1!~/^GO:\d+$/) {
	$annot1=&terms_to_GOs($_[0],0);
    }
    if ($annot2!~/^GO:\d+$/) {
	$annot2=&terms_to_GOs($_[1],0);
    }
    my $TESTS=$_[2]||1;
    my $pinter=$_[3]||undef;  ## call the sub with 4 arguments to return inter prob
    return(1) if $annot1 eq $annot2;
    my $gopair=join("_",sort {$b lt $a} ($annot1,$annot2));
    ############################
    # Declare database details #
    ############################
    my $host=$db_host;
    my $database="pogo";
    my $user="root";
    my $pw="yenapas";
    my $table="$species";
    $table="inter_$table" if $pinter;
    my $port=$db_port;
    ###############################################
    # If this is the first run, connect to the DB #
    ###############################################
    if (scalar keys (%cc_prob)==0) {
	$dsn = "DBI:mysql:database=$database;host=$host;port=$port;"; 
	##Select database 
	$dbh=DBI->connect($dsn,$user,$pw);
    }
    if (defined($cc_prob{$gopair})){
	return($cc_prob{$gopair}); 
    }
    else {
	my $query="SELECT P_low from $table where gopair='$gopair'";  
	debug($query);
	my $result=$dbh->prepare("$query");        
	## Run query
	$result->execute;                                          
	my $ref = $result->fetchrow_hashref();
	$cc_prob{$gopair}=$ref->{'P_low'} * $TESTS;
	$cc_prob{$gopair}=1 if $cc_prob{$gopair}>1;
    }
    $cc_prob{$gopair}=1 if $cc_prob{$gopair}>1;
    $cc_prob{$gopair}=1 unless defined($cc_prob{$gopair});
    return($cc_prob{$gopair});
    
}
#####################################################################
sub usage{
    $0=~/.+\/(.+)/;
    my $name = $1;
    my $us="[options] <ANNOTATIONS FILE> <NETWORK FILE>";
    my $desc="This script will parse the output of moonGO.pl. It expects files of the format foo.<MODE>.bar, where mode is one of the modes of moonGO (abBcimpPd). For example, human.b.moon.";
    my %opts=(
	      "usage" => $us,
	      "desc" => $desc,
	      "B" => "(B)ase directory. The data directory is assummed to be ./data and this is where the network files are found. Def: ~/research/new_moon/data",
	      "c" => "Minimum class (c)ardinality for mode i. Only candidates involving classes with  at least this many members will be returned.",
	      "C" => "(C)lass file. Def: <data dir>/<species>.clas",
	      "d" => "Debug, very very verbose",
	      "D" => "Data directory, this is where the script expects to find data files. Def:<BASE_DIR>/data.", 
	      "f" => "Return those candidates that were (f)ound by at least this many methods. So, '$name -f 3' will return those cands found by at least 3 different methods.",
	      "F" => "File containing a list of Gold Standard proteins (def: <DATA_DIR>dir>/gold_simple.txt).",
	      "g" => "Count (g)os. It will print a tab separated list of: \n" .
	      "\t\t1. each GO found in the file(s); \n" .
	      "\t\t2. the total number of protein candidates; \n" .
	      "\t\t3. the number of cases, ie the total number of cases for all candidates. \n" .
	      "\t\t   Some cands can be found by more then one method and more than once by each method; \n" .
	      "\t\t4. the total number of GOs found in the file(s);\n" .
	      "\t\t5. the number of times this GO was found and\n" .
	      "\t\t6. the percentage of the results (all cases) that involved this GO.",
	      "G" => "Only print results for (g)old standard proteins.",
	      "h" => "Print this (h)elp and exit.",
	      "l" => "Probability of association threshold for Ce(l)lular Component annotations. Def: no limit",
	      "L" => "Class file annotated for Cellular Component ontology \n\t\t\t[default: <DATA_DIR>/<species>.CC.clas].",
	      "m" => "The (m)odes I want (def = 'all'). If output files of other modes are given, they will be ignored.",
	      "M" => "For each candidate, print the (m)ethods it was found by.",
	      "n" => "Build (n)etworks and cytoscape scripts for each candidate. This will create the following files:\n" .
	      "\t\t<cand_name>.<mode>.<case_number>.sif : Network file.\n" .
	      "\t\t<cand_name>.<mode>.<case_number>.class.attr : Classes attribute file\n" .
	      "\t\t<cand_name>.<mode>.<case_number>.true_class.attr : True classes attribute file\n" .
	      "\t\t<cand_name>.<mode>.<case_number>.names.attr : Names attribute file\n" .
	      "\t\t<cand_name>.<mode>.<case_number>.cyto : Cytoscape script",            
	      "N" => "Network file. Def: <BASE_DIR>/<SPECIES>.gr",
	      "P" => "Print Count GO pairs (like -g above, but for pairs).",
	      "p" => "Maximum (p)robability threshold. Only cands with p<= this threshold will be printed.",
	      "t" => "Database port for connecting to PoGO (def: 3306).",
	      "H" => "Host name/ip for the PoGO database (def: 10.1.3.30).",
	      "r" => "Minimum p(r)ecision threshold. Only cands with AT LEAST 1 term with r>= this threshold will be printed.",
	      "R" => "Minimum p(r)ecision threshold. Only cands BOTH of whose terms have r>= this threshold will be printed.",
	      "s" => "The (s)um of the precisions of each go term. Only cases whose go precision values add up to more than the value given will be printed.",
	      "T" => "GO term id file. Def: <DATA_DIR>/GO.terms_alt_ids",
	      "v" => "Verbose, print progress reports",
	      "V" => "Verbose output, print more details for candidates",
	      "S" => "Species, def: guessed from input file name",

	      "w" => "List of wanted proteins. Only results involving these will be printed. It can be either a file with one name per line or a comma separated list.",
	      "x" => "List of proteins or terms to e(x)clude. Any results involving these will be ignored. It can be either a file with one name per line or, a comma separated list.",
	      "X" => "Produce .csv style format for import into Excel or libreoffice spreadsheets"
	     );
    print_help_exit(\%opts,0);
    exit();
}
