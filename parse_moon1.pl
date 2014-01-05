#!/usr/bin/perl -w
# use lib "/home/terdon/lib";
# use moon_cand;
use strict;
use Getopt::Std;
use Switch;
use Term::ANSIColor; 
use Math::Round;

my (%opts, %terms_to_GOs, %found,%gold,%found_gos,%precision);
my $have_already_read_terms_file=0;

getopts('MTgncnN:w:f:C:m:p:r:R:',\%opts);
my $color="bold blue";
my $only_gold=$opts{g}||undef;
my $wanted_modes=$opts{m}||'all';
my $count_methods=$opts{M}||undef;
my $min_prob=$opts{p}||1;
my $sort_by=$opts{s}||'2spec';
my $found_by=$opts{f}||undef;
my $count_gos=$opts{c}||undef;
my $DATADIR=$opts{D}||"/home/terdon/research/testing/new/data";
my $class_file=$opts{C}||"$DATADIR/human.clas";
my $network_file=$opts{N}||"$DATADIR/../HSHQ.gr";
my $go_terms_file=$opts{T}||"$DATADIR/GO.terms_alt_ids";
my @MODES;
my $least_likely=1 if $opts{p};
my $gold_file=$opts{G}||"$DATADIR/gold_simple.txt";
my $wanted=$opts{w}||undef;
my @files=@ARGV;
my %candidates;
my %fcands;
my $min_prec=0;
if($opts{r}){
    $min_prec=$opts{r};
}
if($opts{R}){
    $min_prec=$opts{R};
}
if($opts{r} && $opts{R}){die "-r means that EITHER go's precision must be below the threshold, -R that BOTH must be. The two options cannot be used together\n"}



&check_gold(); 
my %want;
if($wanted){
    if(-e $wanted){
	open(W,"$wanted")||die("Could not open file $wanted:$!\n");
	while(<W>){
	    chomp;
	    next if /^Name\t/;
	    /^(.+?)\s/;
	    $want{$1}++;
	}
    }
    else{
	my @s=split(/\s/,$wanted);
	map{$want{$_}++}@s;
    }
    die("Error, need a list of names, either a file or a space separated list with option -w\n") unless scalar(keys(%want))>0;
}

if($count_gos){&count_gos(\@ARGV)}
else{
    foreach my $file (@ARGV){
	$file=~/.+\.([BabimpPcd])\./;
	my $mode=$1;
	if($wanted_modes ne 'all'){
	    next unless $wanted_modes=~/$mode/;
	}
	push @MODES,$mode;
	open(A,"$file")||die("Cannot open $file : $!\n");
	while(<A>){
	    next if /^\s*$/;
	    chomp;
	    my @a=split(/\t/);
	    $found{$a[0]}{$mode}=1;
	    ## Check what type of file this is
	    switch($mode){
		case /[i]/{
		    $fcands{$a[0]}{NUMS}++;
		    $fcands{$a[0]}{MODE}{$mode}++;
		    $a[2]=~s/\s*$//;
		    $a[3]=~s/\s*$//;
		    $a[1]=~/(\d+),\s*(\d+)/;
		    ## I use this to deal with cases where the same cand
		    ## appears in a file >1 time
		    $fcands{$a[0]}{$mode}{$fcands{$a[0]}{MODE}{$mode}}{CLASS1}=$1;
		    $fcands{$a[0]}{$mode}{$fcands{$a[0]}{MODE}{$mode}}{CLASS2}=$2;
		    #$fcands{$a[0]}{$mode}{$fcands{$a[0]}{MODE}{$mode}}{go1}="\"" . $a[2] . "\"," . &terms_to_GOs($a[2],0)  ;
#		    $fcands{$a[0]}{$mode}{$fcands{$a[0]}{MODE}{$mode}}{go2}="\"" . $a[3] . "\"," . &terms_to_GOs($a[3],0)  ;
		    $fcands{$a[0]}{$mode}{$fcands{$a[0]}{MODE}{$mode}}{go1}= &terms_to_GOs($a[2],0) ;
		    $fcands{$a[0]}{$mode}{$fcands{$a[0]}{MODE}{$mode}}{go2}= &terms_to_GOs($a[3],0) ;
		    $found_gos{&terms_to_GOs($a[2],0)}++;
		    $found_gos{&terms_to_GOs($a[3],0)}++;
		    
		    $a[4]=~s/\s*:\s*//;
		    $fcands{$a[0]}{$mode}{$fcands{$a[0]}{MODE}{$mode}}{SCORE}=nearest(0.000001,$a[4]);
		    $a[4]=~/(\d+\.\d{0,3})/;
		    # my $o=$1;
		    # my $b="";
		    # if($a[4]=~/\d+(e\-\d+)$/){ $b=$1;}
		    # ## catch probs like 0.0000000002
		    # if($a[4]!~/e/ && $b eq ''){
		    # 	$fcands{$a[0]}{$mode}{$fcands{$a[0]}{MODE}{$mode}}{SCORE}=$a[4];
		    # }
		    # else{
		    # 	$fcands{$a[0]}{$mode}{$fcands{$a[0]}{MODE}{$mode}}{SCORE}=$o . $b;
		    # }
		}
		case /[pP]/{
		    $fcands{$a[0]}{NUMS}++;
		    $fcands{$a[0]}{MODE}{$mode}++;
		    $a[3]=~s/\s*:\s*//;
		    $a[1]=~/^(.+)\s\((GO:\d+)/;
		    $fcands{$a[0]}{$mode}{$fcands{$a[0]}{MODE}{$mode}}{INTER1}=$1;
		    $fcands{$a[0]}{$mode}{$fcands{$a[0]}{MODE}{$mode}}{go1}=$2 ;
		    $found_gos{$2}++;
		    $a[2]=~/^(.+)\s\((GO:\d+)/;
		    $fcands{$a[0]}{$mode}{$fcands{$a[0]}{MODE}{$mode}}{INTER2}=$1;
		    $fcands{$a[0]}{$mode}{$fcands{$a[0]}{MODE}{$mode}}{go2}=$2;
		    $found_gos{$2}++;

		    $a[3]=~s/\s*:\s*//;
		    $fcands{$a[0]}{$mode}{$fcands{$a[0]}{MODE}{$mode}}{SCORE}=nearest(0.000001,$a[3]);

		}
		case /[Bb]/{
		    #	for(my $i=0;$i<=$#a;$i++){print "$i:$a[$i]\n"}die();
		    $fcands{$a[0]}{NUMS}++;
		    $fcands{$a[0]}{MODE}{$mode}++;

		    $a[1]=~/(\d+)\s+(GO:\d+)/;
		    $fcands{$a[0]}{$mode}{$fcands{$a[0]}{MODE}{$mode}}{CLASS1}=$1;
		    $fcands{$a[0]}{$mode}{$fcands{$a[0]}{MODE}{$mode}}{go1}=$2;
		    $found_gos{$2}++;
		    $a[2]=~/^(\d+)\s+(GO:\d+)/;
		    $fcands{$a[0]}{$mode}{$fcands{$a[0]}{MODE}{$mode}}{CLASS2}=$1;
		    $fcands{$a[0]}{$mode}{$fcands{$a[0]}{MODE}{$mode}}{go2}=$2;
		    $found_gos{$2}++;
		    $a[3]=~s/\s*:\s*//;
		   #for (my $i=0; $i<=$#a;$i++){print "$i : $a[$i]\n";}
		    $fcands{$a[0]}{$mode}{$fcands{$a[0]}{MODE}{$mode}}{SCORE}=nearest(0.000001,$a[3]);
		}
		case /[cdm]/{
		    die("Not ready for mode $mode yet\n");
		    my @b=split(/\s+/,$a[0]);
		    $fcands{$b[0]}++;
		    @b=split(/\s+/,$a[2]);
		    $fcands{$b[0]}{NUMS}++;
		    $fcands{$b[0]}{MODE}{$mode}++;
		    $fcands{$b[0]}{$mode}{$fcands{$a[0]}{MODE}{$mode}}{SCORE}=$a[3];
		    $a[3]=~s/\s*:\s*//;
		    $fcands{$a[0]}{$mode}{$fcands{$a[0]}{MODE}{$mode}}{SCORE}=nearest(0.000001,$a[3]);
		}
		case /a/{
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
		    $a[4]=~s/\s*:\s*//;
		    $fcands{$a[0]}{$mode}{$fcands{$a[0]}{MODE}{$mode}}{SCORE}=nearest(0.000001,$a[4]);
		}
	    } ## end switch
	} ## end while(<A>)
   } ## end foreach my $file

    ## get precision values for all GOs found
    print STDERR "Precision...";
    &term_precision(0);
    print STDERR "Done\n";
    my $tot=0;


    ## If we just want those candidates that were found by
    ## at least $found_by methods
    if($found_by){
	cand:foreach my $cand (keys(%fcands)){
	    if($wanted && not defined($want{$cand})){ next cand;}
	    if($only_gold){
		next unless defined($gold{$cand}) ;
	    }
	    next unless scalar(keys(%{$found{$cand}}))>=$found_by;
	    $tot++;
	    my $gcand=color("$color").$cand.color("reset");
	    defined($gold{$cand}) ? (print "$gcand\t") : (print "$cand\t");
	    
	    my @modes=sort {uc($a) cmp uc($b)} keys(%{$fcands{$cand}{MODE}});
	    map{print "$_ " }@modes;
	    print "\n";
	    map{&print_results($cand,$_)}@modes;
	    if($opts{n}){
		foreach my $mode (@modes){
		    next unless $mode =~/[aBbi]/;
		    next unless defined($found{$cand}{$mode});
		    for (my $i=1; $i<=$fcands{$cand}{MODE}{$mode}; $i++){
			&make_networks($cand,$mode,$i,$fcands{$cand}{$mode}{$i}{CLASS1},$fcands{$cand}{$mode}{$i}{CLASS2});
		    }
		}
	    }
	}
	  print "TOTAL : $tot\n";
    }
    ## Just print what methods found each candidate 
    elsif($count_methods){
	foreach my $cand (keys(%fcands)){
	    if($only_gold){
		next unless defined($gold{$cand}) ;
	    }
	    if($wanted && not defined($want{$cand})){ next;}
	    my $gcand;
	    defined($gold{$cand}) ? ($gcand=color("$color").$cand.color("reset")) : ($gcand=$cand);
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
	my @modes;
	my %to_be_printed;
	if($wanted_modes ne 'all'){
	    @modes=split(//,$wanted_modes);
	}
	else{@modes=@MODES;}
	foreach my $cand (keys(%fcands)){
	    print "$cand\n";
	    if($only_gold){
		next unless defined($gold{$cand}) ;
	    }
	    if($wanted && not defined($want{$cand})){ next;}
	    my $gcand;
	    defined($gold{$cand}) ? ($gcand=color("$color").$cand.color("reset")) : ($gcand=$cand);
	    my $printed=0;
	    foreach my $m ((keys(%{$fcands{$cand}{MODE}}))){
		## for each instance of this cand in this file
		for (my $i=1; $i<=$fcands{$cand}{MODE}{$m}; $i++){
		    my $prec1=&term_precision(1,$fcands{$cand}{$m}{$i}{go1})||1;
		    my $prec2=&term_precision(1,$fcands{$cand}{$m}{$i}{go2})||1;
		    next unless $fcands{$cand}{$m}{$i}{SCORE}<=$min_prob;
		    if($opts{R}){			
			next unless $prec1>=$min_prec;
			next unless $prec2>=$min_prec;
		    }
		    elsif($opts{r}){
			unless( $prec1>=$min_prec ||
				$prec2>=$min_prec){
			    next;
			}
		    }
		    $tot++;
		    # if($printed==0){
		    # 	push @to_be_printed, "$gcand\n";
		    # 	$printed=1;
		    # }
		    $to_be_printed{$cand}{NAME}=$gcand;
		    map{$to_be_printed{$cand}{MODE}{$_}++}@modes;
		    if(($opts{n}) && ($m eq 'i' || $m eq 'B' || $m eq 'P' || $m eq 'b')){ 
			for (my $i=1; $i<=$fcands{$cand}{MODE}{$m}; $i++){
			    &make_networks($cand,$m,$i,$fcands{$cand}{$m}{$i}{CLASS1},$fcands{$cand}{$m}{$i}{CLASS2}) if defined($found{$cand}{$m});
			}
		    }
		    elsif(($opts{n}) && ($m eq 'a')){
			for (my $i=1; $i<=$fcands{$cand}{MODE}{$m}; $i++){
			    &make_networks($cand,$m,$i,$fcands{$cand}{$m}{$i}{CLASS1}) if defined($found{$cand}{$m});
			}
		    }
		} # end for $i
	    } ## end foreach $m
	} # end foreach cand
	 print "TOTAL : $tot\n";
    }
}

############################################################
sub make_networks{
    my $cand=shift;
    my $mode=shift;
    my $i=shift;
#    my $cand_attr=$cand . "." . $mode . ".cand.attr";
 #   open(A,">$cand_attr");
    if(($mode eq 'i')||($mode eq 'B') || ($mode eq 'b') || ($mode eq 'a')){
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
	    open(B,"moon_get_class_prots.pl $class2 $class_file |");
	    while(<B>){
		chomp;
		next if $cand eq $_;
		$names{$_}{$class2}++;
	    }
	    close(B);
	}
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
		print C1 "$NetName = $class1\n";
	    }
	    elsif(defined($names{$name}{$class2})){
		print C "$NetName = B\n";
		print C1 "$NetName = $class2\n";
	    }
	    else{die("Huh!!??\n$NetName,$class1,$class2,$class1,$class2\n");}
	}
	close(C);
	my $subnet_file=$cand . ".$i." . $mode . ".sif";
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
	
	my $cyto_script=$cand . ".$i." . $mode . ".cyto";
	open(C3,">$cyto_script")||die "Could not open $cyto_script:$!\n";
	my $cur_dir=`pwd`;
	chomp($cur_dir);
	print C3 "session open file=\"/home/terdon/research/testing/new/results/candidates.cys\"\n";
	print C3 "network import file=\"$cur_dir/$subnet_file\"\n";
	print C3 "node import attributes file=\"$cur_dir/$attr_file\"\n";
	print C3 "node import attributes file=\"$cur_dir/$attr_file1\"\n";
	print C3 "node import attributes file=\"$cur_dir/$attr_file2\"\n";

    } ## end if mode i || b
    # elsif($mode eq 'a'){
    # 	my $class1=shift;
    # 	my %names;
    # 	open(A,"moon_get_class_prots.pl $class1 $class_file |");
    # 	while(<A>){
    # 	    chomp;
    # 	    next if $cand eq $_;
    # 	    $names{$_}{$class1}++;
    # 	}
    # 	close(A);
    # }
}


############################################################
sub count_gos{
    my @files=@{$_[0]};
    my %gos;
    my $total_cands=0;
    foreach my $file (@files){
	$file=~/.+\.([BabimpPcd])\./;
	my $mode=$1;
	
	open(A,"$file")||die("Cannot open $file : $!\n");
	while(<A>){
	    next if /^\s*$/;
	    chomp;
	    $total_cands++;

	    if($mode eq 'i'){
		my @a=split(/\t/);
		$a[2]=~s/\s*$//;
		$a[3]=~s/\s*$//;
		$gos{&terms_to_GOs($a[2],0)}++;
		$gos{&terms_to_GOs($a[3],0)}++;
		
	    }
	    else{
		my @a=/(GO:\d+)/g;
		$gos{$a[0]}++;
		$gos{$a[1]}++;
	    }
	}
    }
    my @GOs=keys(%gos);
    my $tt=0;
    map{$tt+=$gos{$_}}@GOs;

    foreach my $go (keys(%gos)){
	print "$go\t$gos{$go}/$total_cands/$#GOs/$tt\t";
	printf("%.2f",($gos{$go}*100)/$tt);
	print "\t"  . &terms_to_GOs($go,1) . "\n";
	
    }


}

############################################################
    sub terms_to_GOs{
	my $term=shift;
	my $mode=shift; ## 0 will return GO:xxx, 1 will return term name
	$term=~s/_/ /g;
	my $subonto="P";
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

sub print_results{
    my $cand=shift;
    my $m=shift;
    my $gcand;
    defined($gold{$cand}) ? ($gcand=color("$color").$cand.color("reset")) : ($gcand=$cand);
    if($m eq 'P'){
	if (defined($found{$cand}{$m})){
	    for (my $i=1; $i<=$fcands{$cand}{MODE}{$m}; $i++){
		next unless $fcands{$cand}{$m}{$i}{SCORE}<=$min_prob;
		my $prec1=&term_precision(1,$fcands{$cand}{$m}{$i}{go1});
		my $prec2=&term_precision(1,$fcands{$cand}{$m}{$i}{go2});
		if($opts{R}){			
		    next unless $prec1>=$min_prec;
		    next unless $prec2>=$min_prec;
		}
		elsif($opts{r}){
		    unless( $prec1>=$min_prec ||
			    $prec2>=$min_prec){
			next;
		    }
		}
		
		print "\tMODE:P\tInter1: $fcands{$cand}{$m}{$i}{INTER1}\tInter1_GO: $fcands{$cand}{$m}{$i}{go1} ($prec1) \tInter2: $fcands{$cand}{$m}{$i}{INTER2}\tInter2_GO: $fcands{$cand}{$m}{$i}{go2} ($prec2) \t$fcands{$cand}{$m}{$i}{SCORE}\n";
	    }
	}
    }    
    else{
	if (defined($found{$cand}{$m})){
	    for (my $i=1; $i<=$fcands{$cand}{MODE}{$m}; $i++){
		next unless $fcands{$cand}{$m}{$i}{SCORE}<=$min_prob;
		my $prec1=&term_precision(1,$fcands{$cand}{$m}{$i}{go1})||1;
		my $prec2=&term_precision(1,$fcands{$cand}{$m}{$i}{go2})||1;
		if($opts{R}){			
		    next unless $prec1>=$min_prec;
		    next unless $prec2>=$min_prec;
		}
		elsif($opts{r}){
		    unless( $prec1>=$min_prec || $prec2>=$min_prec){
			next;
		    }
		}
		print "\tMODE:$m\tCandGO: \"" . &terms_to_GOs($fcands{$cand}{$m}{$i}{go1},1) . "\" ($fcands{$cand}{$m}{$i}{go1}, $prec1) Class: $fcands{$cand}{$m}{$i}{CLASS1}\tClassGO: \"" . &terms_to_GOs($fcands{$cand}{$m}{$i}{go2},1) . "\" ($fcands{$cand}{$m}{$i}{go2}, $prec2) \t$fcands{$cand}{$m}{$i}{SCORE}\n";
	    }
	}
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
#    print STDERR "1G : $bait,$bgo,$class_annot,$P,$class\n";

#    if(($mode eq 'a') || ($mode eq 'dis2class')){
    ($bgo,$class_annot,$P,$class)=@_;
    $class='' unless $class;
#    print STDERR "G : $bait,$bgo,$class_annot,$P,$class\n";
    push@{$gold{$bait}{MISSED}{$mode}},"$class\t$class_annot\t$bgo\t$P";
 #   }
    

}
#####################################################################
sub term_precision{
    my $mode=shift;
    if($mode==0){
	open(PP,">/tmp/$$")|| die("Could not open precision file /tmp/$$ for writing:$!\n");
	map{print PP "$_\n";}keys(%found_gos);
	close(PP);
	system("term_precision.pl $DATADIR/biological_process.genealogy  /tmp/$$ > /tmp/$$.b");
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
	return(nearest(.0001,$precision{$go}));
    }
}


#  malaka, fix mode P, only GO:, no term name
#malaka fix multiple modes printed  twice, eg RPB1_HUMAN
	# MODE:i	CandGO: "signal transduction" (GO:0007165, 0.4542) Class: 43	ClassGO: "transcription-coupled nucleotide-excision repair" (GO:0006283, 0.9011) 	0.001547
	# MODE:b	CandGO: "cell death" (GO:0008219, 0.4565) Class: 411	ClassGO: "histone H2A acetylation" (GO:0043968, 0.8821) 	0.411928
	# MODE:i	CandGO: "signal transduction" (GO:0007165, 0.4542) Class: 43	ClassGO: "transcription-coupled nucleotide-excision repair" (GO:0006283, 0.9011) 	0.001547
	# MODE:b	CandGO: "cell death" (GO:0008219, 0.4565) Class: 411	ClassGO: "histone H2A acetylation" (GO:0043968, 0.8821) 	0.411928
