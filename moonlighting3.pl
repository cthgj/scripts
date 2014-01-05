#!/usr/bin/perl -w
use strict;
use Getopt::Std;
use IO::File;
use Cwd; 
use Data::Dumper;

my (%known_function,%ancestors,%candidates,%prob,%skip_this,%missing,%singletons,%terms_to_GOs,%foundGOs,%proteins,%synonyms,%opts,%interactions,%found,%offspring, %interacting,%related);
#my (@interactions);
$synonyms{LOADED}=0;
my $DATADIR=$opts{D}||"./";
getopts('gdVvhmD:S:p:P:G:s:f:',\%opts) || do { print "Invalid option, try 'moonlighting.pl -h' for more information\n"; exit(1); };
my $very_verbose=$opts{V}||undef;
my $have_already_read_terms_file=0;
my $annotation_file_type;
my $gaf_annotations_file=$opts{f} || "gene_association.goa_human";
&usage() if $opts{h};
my @net_probs;
my $stats_dir="$DATADIR/gostats";
my $geneology_file=$opts{g}||"biological_process.genealogy";
my $cwd=cwd();
my $debug=$opts{d}||undef;
my $help=$opts{h}||undef;
my $synfile=$opts{s}||"synonyms_human";
my $annotations_file=$ARGV[0]|| die("Need an annotations file\n");
my $network_file=$ARGV[1]|| die("Need a network file\n");
my $low_prob=$opts{p}||0.00005;
my $high_prob=$opts{P}||0;
my $subonto=$opts{o}||"P";
my $species=$opts{S}|| "human";
my $verbose=$opts{v}||undef;
my $go_terms_file=$opts{G}||"GO.terms_ids_obs";
my $go_percentages=undef;
my $pairs=1;
if($opts{m}){$pairs=undef; $go_percentages=1};


if($subonto eq "p"){$subonto="P"}
elsif($subonto eq "c"){$subonto="c"}
elsif($subonto eq "f"){$subonto="f"}
elsif($subonto eq "F" || "P" || "C"){}
else{print STDERR "Subontology (-o) must be one of \"p,c,f\"\n"; exit(1)}
$verbose=1 if $debug;
my $date=localtime;
my $LOG = IO::File->new("> $$.log")|| die("Cannot open $network_file : $!\n");
print $LOG "********** $date ******************
Network file : $ARGV[1]
Class file: $ARGV[0]\n"; 
map{print $LOG "$_ : $opts{$_}\n"}keys(%opts);
print $LOG "\n\n"; 
## Read the network file
&load_network($network_file);
&load_ancestors();
## Load functional annotations
&parse_gaf_file();
&load_annotations($annotations_file);



## Do your thang!
&main();
system("cat $$.log >> moonlighting.log");
unlink("$$.log");

############################ SUBROUTINES ############################ 

sub main{
    print STDERR "Starting main...\n" if $verbose;

    if($pairs){
	my $koko=0;
	my $pp=scalar(keys(%interactions));
      bait:foreach my $bait (keys(%interactions)){
	  next bait if defined($skip_this{$bait});
	  my (@bait_gos,@target_gos,@target_names);
	  my(%bait_hash,%bait_go_count,%target_go_count)=();
	  $koko++;      
#	  printf(STDERR "$koko of $pp\r") if $verbose;
	  ##get all gos for bait protein
	hash:foreach my $h (@{$proteins{$bait}{$subonto}}){
	    if ($h eq "none"){
		push @bait_gos,"none";
	    }
	    else{
		%bait_hash=%{$h};
		&debug("baits go $bait ",$bait_hash{"goID"},$bait_hash{"Qualifier"} );
		$bait_go_count{$bait}{$bait_hash{"goID"}}++;
		next hash if $bait_hash{"Qualifier"} eq "NOT";
		## Some Gos are mentioned twice in the annotation file
		push @bait_gos,$bait_hash{"goID"} unless $bait_go_count{$bait}{$bait_hash{"goID"}}>1;
	    }
	} ## end foreach bait_hash
	  
	  ## Cycle through all the targets of this bait protein
	target:foreach my $target (@{$interactions{$bait}}){
	    next target if defined($skip_this{$target});
	    @target_gos=();
	  thash:foreach my $h (@{$proteins{$target}{$subonto}}){
	      if ($h eq "none"){
		  push @target_gos,"none";
		  push @target_names,$target;
	      }
	      else{	
		  my %target_hash=%{$h};
		  next thash if $target_hash{"Qualifier"} eq "NOT";

		  ## If the bait protein has even one of the target's GOs, 
		  ## skip target
		  if (exists($bait_go_count{$bait}{$target_hash{"goID"}})){
		      &debug("skipped $target because $bait shares $target_hash{goID}");
		      print $LOG "skipped $target because $bait shares $target_hash{goID}\n";
		      next target ;
		  }
		  foreach my $baitgo (@bait_gos){
		      ## If one of the bait protein's GOs is 
		      ## related to one of the targets, skip target
		      if (exists($offspring{$baitgo}{$target_hash{"goID"}}) || 
			  exists($offspring{$target_hash{"goID"}}{$baitgo})){
			  &debug("$baitgo is related to $target_hash{goID}");
			  print $LOG "aaa $bait $baitgo is related to $target $target_hash{goID}\n";
			  next target ;
		      }
		      ## Else, keep bait/target pair as a candidate
		      ## if the GO pair passes the thresholds  
		      else{
			  my $GOpair=$baitgo . "xxx" .$target_hash{"goID"};
			  ## Skip if either protein has no GO
			  next if $GOpair =~/none/;
			  my $P=&gopair_probability($baitgo,$target_hash{"goID"},0,$bait,$target);
			  if(($P>=$high_prob) && $P<=$low_prob){##fixme
			      $candidates{$bait}{$target}{STRING}=$GOpair . "xxx" . $bait_hash{"source"} . "xxx" . $target_hash{"source"};
			      $candidates{$bait}{$target}{P}=$P;
			  }
			  else{
			      print $LOG "Skipped $target, probability $GOpair = $P\n";
			      next target;
			  }
		      }
		  } ## end foreach $baitgo
	      }## end else
	  }## end foreach target_hash 	
	}## end forach $target
      }## end foreach bait


    }## end if $pairs
    elsif($go_percentages){
	my $koko=0;
	
	my %seen_probs=();
	my $pp=scalar(keys(%interactions));
      bait:foreach my $bait (keys(%interactions)) {
	  my (@bait_gos,@target_gos,@target_names);
	  my %bait_info;
	  next bait if defined($skip_this{$bait});
	  ## Skip unless this bait has enough interactors
	  ## to allow indentification of the odd one(s) out
	  next bait if scalar(@{$interactions{$bait}})<3;
#	  printf(STDERR "$koko of $pp\r") if $verbose;
	  $koko++;      
	  
	  &debug("\n==================BAIT : $bait===================\n");
	  my %bait_hash;
	  my %target_hash;
	bait_hash:foreach my $h (@{$proteins{$bait}{$subonto}}){
	    if ($h eq "none"){
		$bait_hash{"goID"}=$h;
	    }
	    else{
	    	%bait_hash=%{$h};
	    }
	    push @bait_gos, $bait_hash{"goID"};
	    push @{$bait_info{$bait}{GOs}}, $bait_hash{"goID"};
	} ## end foreach bait_hash
	target:foreach my $target (@{$interactions{$bait}}){
	    next target if defined($skip_this{$target});
	    push @{$bait_info{$bait}{TARGETS}},$target;
	  thash:foreach my $h (@{$proteins{$target}{$subonto}}){
	      if ($h eq "none"){
		  push @target_gos,"none";
		  push @target_names,$target;
		  $target_hash{NAME}=$target;
	      }
	      else{	
		   %target_hash=%{$h};
		   $target_hash{NAME}=$target;
		  next thash if $target_hash{"Qualifier"} eq "NOT";
		#  foreach my $baitgo (@bait_gos){	 
		      # ## If one of the bait protein's GOs is 
		      # ## related to one of the targets, skip target
		      # if (exists($offspring{$baitgo}{$target_hash{"goID"}}) || 
		      # 	  exists($offspring{$target_hash{"goID"}}{$baitgo})){
		      # 	  &debug("$baitgo is related to $target_hash{goID}");
		      # 	  print $LOG "aaa $bait $baitgo is related to $target $target_hash{goID}\n";
		      # 	  next target ;
		      # }
		      
	#	  }
		   push @target_gos, $target_hash{"goID"};
		   push @target_names,$target;
#		   push @{$bait_info{$bait}{$target}{GOs}},$target_hash{"goID"};
		   $bait_info{$bait}{$target}{GOs}{$target_hash{"goID"}}++;
	      }
	  } ## end thash
	   
	} ## end foreach target
	  #&calculate_percentages($bait,\@bait_gos,\@target_gos,\@target_names,\%bait_hash,\%target_hash,$koko,$pp) if scalar(@target_gos)>3;
	  my @cands=&calculate_percentages(\%bait_info,$koko,$pp);
	  print "cc $bait: @cands\n";
      } ## end foreach bait
    } ## end elsif $go_percentages
    print STDERR "\n" if $verbose;
    ## Now, print those bait/target pairs that pass the thresholds

    foreach my $bait (keys(%candidates)){print "aa $bait\n";
	foreach my $target (keys(%{$candidates{$bait}})){
	    if ($go_percentages){
		print "$candidates{$bait}{$target}{STRING}\n";
	    }
	    else{
		my @a=split(/xxx/,$candidates{$bait}{$target}{STRING});
		print "$bait $a[0]\t$a[2]\t$target\t$a[1]\t$a[3]\t" . $candidates{$bait}{$target}{P} . "\n";	
	    }	
	}
    }

    
}## end main
## Create probalities file for this network,so that we can run faster next time
unless (-e "$network_file.prob"){
    open(NET,">>$network_file.prob");
    map{print NET}@net_probs;
    system("cat $network_file.prob | sort|uniq > lililolo; mv lililolo $network_file.prob");
    close(NET);
}

############################################################
sub calculate_percentages{
    my %hash=%{(shift)};
    my $koko=shift;
    my $pp=shift;
    my ($mc,$mcgo,$perc)=0;
    my %gostats;
    my @keys=keys(%hash);
    my %good; ## this will hold the targets we are interested in     
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
#	    printf(STDERR "$bait<==>$target\n") if $verbose;
	    ## skip target if it is annotated to $mcgo
	    next target if defined($hash{$target}{GOs}{$mcgo});
    	    foreach my $tgo (keys(%{$hash{$target}{GOs}})){
		my $P=&gopair_probability($mcgo,$tgo,0,$bait,$target);
		## skip target if any of its gos are related to the mcgo
		if($P>=$low_prob){##fixme
		    $good{$target}=undef;
		    next target;
		}
		foreach my $bgo (@{$hash{GOs}}){
		    $P=&gopair_probability($bgo,$tgo,0,$bait,$target);
		    if($P>=$low_prob){
			## skip target if ANY of its gos 
			## are close to any bgo
			$good{$target}=undef;
			next target;
		    }
		    else{
			$good{$target}++;
		    }
		}
	    } ## end foreach my tgo
	} #end foreach target
	  
    } ## end if perc >=50
    return(keys(%good));
}
############################################################

sub load_network{
    print STDERR "Loading Network...\n" if $verbose;
     #open($A,"$network_file")|| die("Cannot open $_[0]:$!\n");
    my $A = IO::File->new("< $network_file")|| die("Cannot open $network_file : $!\n"); # even better!
    my $n=0;
    while(<$A>){
	next if /^\d+$/;
	my $line=$_; #$_ is empty after &get_name
	my ($bait,$target)=split(/\s+/,$_);
	$bait=&get_name($bait,"load_network");
	$target=&get_name($target,"load_network");
	$proteins{$bait}{F}[0]=$proteins{$bait}{P}[0]=$proteins{$bait}{C}[0]="none";
	$proteins{$target}{F}[0]=$proteins{$target}{P}[0]=$proteins{$target}{C}[0]="none";
	$found{$bait}=$found{$target}=1;
	#&get_annotations($bait,$target);
	my $a=$bait . "xxx" . $target;
	$interacting{$a}++;
	if(exists($interactions{$target})){
	    push @{$interactions{$target}},$bait;
	}
	else{
	    push @{$interactions{$bait}},$target;
	}
	$n++;
    }
    close(A);
## Make sure that no proteins exist as both baits and targets
## That will happen in a case like this:
## ProtB ProtA
## ProtA ProtC
## ProtA ProtD
    # my %pp;
    # foreach my $bait (keys(%interactions)){
    # 	foreach my $target (@{$interactions{$bait}}){
	    
    # 	    print "$bait : $target\n";
    # 	}
    # }die();
}

############################################################

sub load_annotations{
    print STDERR "Loading annotations...\n" if $verbose;
    my $prot;
    open(A,"$annotations_file")|| die("Cannot open $_[0]:$!\n");
    my $singleton=0;
    my @GOs;
    while(<A>){
	if($.==1){
	    /^\[CLASS/ ? ($annotation_file_type=1) : ($annotation_file_type=0);
	}
	s/\t\t/\tnull\t/g ; ## Deal with empty fields
	s/\t\n/\tnull\n/g ; ## Deal with empty fields
        ## If this is a GO GAF file
	if($annotation_file_type==0){
	    next if /^!/;
	    next if /^\s*$/;
	    chomp;
	    my @tmparray=split(/\t/);
#	    if (scalar(@tmparray)!=18){die("Baaad format annotation file\n@tmparray\n" . scalar(@tmparray) . "\n");};
	    next unless $tmparray[8] eq $subonto;
	    if(defined($found{&get_name($tmparray[1],"found?")})){
		#&debug("Adding(1): @tmparray\n");
		&add_protein(@tmparray,"GAF");
	    }
	    
	}
	## If this is a class file
	else{
	    if(/^CA/){
		@GOs=();
		$singleton=0;
		## if the file has GO ids		
		if(/^CA.+:\d+\)/){die("The file has go ids\n")}
		else{ ## if the file has only go terms
		    /^CA\s+(.+)$/;
		    my $kk=$1;
		    my @ko=split(/\s+/,$kk);
		    foreach my $term(@ko){
			my $T=&terms_to_IDs($term);
			push @GOs,$T;
			$foundGOs{$T}++;
			#&debug( "ff : \$foundGOs{$T}: $foundGOs{$T}");
		    }
		}
	    }
	    if(/^P\#\s+(\d+)/){
		$singleton=1 if $1 == 1;
	    }
	    if(/^PN\s*(.+)/){
		my $kk=$1;
		if($singleton==1){
		    $singletons{$kk}=1;
		}
		my @prots=split(/,\s+/,$kk);
		my %hash;
		foreach my $protein(@prots){
		    for (my $n=0; $n<scalar(@GOs); $n++){
			$hash{"goID"}=$GOs[$n];
#			&debug("$protein also has $GOs[$n]");
			$hash{"Qualifier"}="IS";
			$hash{"source"}="inferred";
			
			push @{$proteins{$protein}{$subonto}},\%hash;
		    }
		}		
	    }		
	}
    }
    close(A);
    my @keys=keys(%found);
    if($annotation_file_type==1){
	for (my $n=0; $n<scalar(@keys); $n++){
	    my $name=$keys[$n];
# Find those proteins that have no class annotation
	    # &debug("111 \$proteins{$name}{$subonto}[0] :  $proteins{$name}{$subonto}[0]" );
	    if($proteins{$name}{$subonto}[0] eq "none"){
		for(my $n=0;$n<scalar(@{$proteins{$name}{$subonto}});$n++){
		  #  &debug("MISSING $n : $name : $proteins{$name}{$subonto}[$n]");
		}
		$missing{$name}++;
	    }
	}
    }
}
############################################################
sub parse_gaf_file{
    print STDERR "Parsing GAF file...\n" if $verbose;
    open(A,"$gaf_annotations_file")|| die("Cannot open $gaf_annotations_file:$!\n");
    while(<A>){
	my $sillyvar="!";
	next if /^$sillyvar/; ## silly thing cause emacs screws up the syntax highlihting when /!/;
	chomp;

	my @tmparray=split(/\t/,$_);
	next unless $tmparray[8] eq $subonto;
	my $nn=&get_name($tmparray[2],"parse_gaf_file");
	
#	map{print "x : $_\n"}keys(%proteins); die();
	if(exists($proteins{$nn})){
	    &add_protein(@tmparray,"GAF");
#	    print "a1a $nn ${$proteins{$nn}{$subonto}}[0]\n"; 
	    $known_function{$nn}{$tmparray[4]}++;
	    map{$known_function{$nn}{$_}++;}@{$ancestors{$tmparray[4]}};
	} 
	
	# foreach my $h (@{$proteins{CEP63_HUMAN}{P}}){
	#     next if $h eq "none";
	#     my %hash=%{$h};
	#     print "a : " . $hash{"goID"} . "\n" if $nn =~ /CEP63/;
	# }
    }
    close(A);
  
}

############################################################
sub load_ancestors{
    open (ANC,"$geneology_file")|| die("cannot open $geneology_file:$!\n");
    while(<ANC>){
	chomp;
	my @terms=split(/\t/);
	my $child=shift(@terms);
	@{$ancestors{$child}}=@terms;
	map{$offspring{$_}{$child}++}@terms;
    }
    close(ANC);
}


####################################################################### 


sub add_protein{
#    &debug("\nAdding : @_\n");
    
    my ($annotation_source)=$_[$#_];
    my($db) = shift;
    my($dbID) = shift;
    my $name=&get_name($dbID,"add_protein for");
        die("$annotation_source\n") if $name=~ "Q7Z6U2_HUMAN"; 
        die("$annotation_source\n") if $dbID=~ "Q7Z6U2_HUMAN"; 
    my($dbSymbol) = shift;
    my($Qualifier) = shift;
    $Qualifier = "IS" if $Qualifier eq "";
    my $goID=shift;
    $foundGOs{$goID}++;
#	push my (@goIDs), $goID;
    my($dbRef) = shift;
    my($evid) = shift;
    my($with) = shift;
    my($aspect) = shift;
    my($dbObjName) = shift;
    my($synonym) = shift;
    my($type) = shift;	
    my($taxon) = shift;
    my($Date) = shift;
    my($AssBy) = shift;
    my($AnnotExt) = shift;
    my($GPFormId) = shift;
    # &debug("===============================================");
    # &debug("name             \t$name ");
    # &debug("db              \t$db");
    # &debug("dbID             \t$dbID");
    # &debug("dbSymbol             \t$dbSymbol");
    # &debug("Qualifier             \t$Qualifier");
    # &debug("goID             \t$goID");
    # &debug("dbRef             \t$dbRef");
    # &debug("evid             \t$evid");
    # &debug("with             \t$with");
    # &debug("aspect             \t$aspect");
    # &debug("dbObjName             \t$dbObjName");
    # &debug("Synonym             \t$synonym");
    # &debug("type             \t$type");
    # &debug("taxon             \t$taxon");
    # &debug("Date             \t$Date");
    # &debug("AssBy             \t$AssBy");
    # &debug("AnnotExt             \t$AnnotExt");
    # &debug("GPFormId             \t$GPFormId") if defined($GPFormId);
    # &debug("source             \t$annotation_source ");
    # &debug("===============================================\n");
#    map{$synonyms{$_}=$dbSymbol;}split(/\|/,$synonym);
#    $synonyms{$name}=$dbSymbol;
    my %hash= (
	"name"          => $name || "null",
	"db"		=> $db||"null",
	"dbID"		=> $dbID||"null",
	"dbSymbol"	=> $dbSymbol||"null",
	"Qualifier"	=> $Qualifier||"null",
	"goID"		=> $goID||"null",
	"dbRef"		=> $dbRef||"null",
	"evid"		=> $evid||"null",
	"with"		=> $with||"null",
	"aspect"	=> $aspect||$subonto,
	"dbObjName"	=> $dbObjName||"null",
	"Synonym"	=> $synonym||"null",
	"type"   	=> $type||"null",
	"taxon"		=> $taxon||"null",
	"Date"		=> $Date||"null",
	"AssBy"		=> $AssBy||"null",
	"AnnotExt"	=> $AnnotExt||"null",
	"GPFormId"	=> $GPFormId||"null",
	"source"        => $annotation_source || "null"
        );
    die("xaxa no $name : $aspect : $subonto " .  $hash{"aspect"} . "\n")  unless defined($proteins{$name}{$aspect});
    if(exists($proteins{$name}{$subonto}[0])){
	if ($proteins{$name}{$subonto}[0] eq "none"){
	    shift(@{$proteins{$name}{$subonto}}); 
	}
    }
    push @{$proteins{$name}{$aspect}},\%hash;
#    map{&debug("$_ : -$hash{$_}-")}keys(%hash);
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
		## network_name\tcurrent_name\tAccession
		if(/^(.+)\t(.+)\t(.+)/){
		   
		    $synonyms{ACC}{$3}=$3;
		    $synonyms{ACC}{$2}=$3;
		    $synonyms{$1}=$2;
		    $synonyms{$3}=$2;
		}	
                ## If this is a type 2 synonyms file, ie, 
		## gene name\tlist of synonyms
		   elsif(/^(.+)\t(.+)/){
		       my $nn=$1;
		       ## Discard protein if there are naming problems
		       if (defined($synonyms{$nn})){
			   $skip_this{$nn}++;
		       }
		       my @syns=split(/\s/,$2);
		       $synonyms{GOOD}{$nn}++;
		       for (my $n=0;$n<scalar(@syns); $n++){
			   $synonyms{$syns[$n]}=$nn;
		       }
		   }
		   else{
		       die("Bad format synonyms file : $synfile\n");
		   }
	    }
	    close(S);
	}
	else{die();
	    print STDERR "No synonyms file found, parsing annotations... \n NEED TO MODIFY THIS BECAUSE OF THE CHANGE IN SYNONYMS FORMAT" if $verbose;
	    open(A,"$ARGV[0]")|| die("Cannot open synonyms $ARGV[0]:$!\n");
	    while(<A>){
		my $silivar="!";
		next if /^$silivar/; ## silly thing cause emacs screws up the syntax highlihting when /!/;
		my @tmparray=split(/\t/);
		map{$synonyms{$_}=$tmparray[2]}split(/\|/,$tmparray[10]);
	    }
	    close(A);
	}	
 	$synonyms{LOADED}=1;
    }
    my $hname;
    $name =~ /_$species/i ? ($hname=$name) : ($hname=$name . "_" . uc($species));
    if(defined($synonyms{GOOD}{$hname})){
	$name=$hname;
	#die("shiiit1 : $name=$hname\n") if $synonyms{GOOD}{$hname} >1;
    }
    elsif(defined($synonyms{$name})){
	$name=$synonyms{$name} ;
    }
    elsif ($name=~/(.+)_$species/i && defined($synonyms{$1})){
	$name=$synonyms{$1};
    }
    else{$synonyms{$name}=$name;}
#    &debug("NNNNNN : $old_name => $name");
    if ($very_verbose){&debug("Called by : $called_by for $old_name ($synfile), returned $name"); }
    return $name;
    
}
############################################################
sub gopair_probability{
    #  $gopair e.g. : GO:0000001xxxGO:0050876
    my $go1=shift;
    my $go2=shift;
    my $tail=shift;
    my ($high,$low);
    my $gopair=$go1 . "xxx" . $go2;
    my $gopair2=$go2 . "xxx" . $go1;
    
    $go1=~/GO:((.....).+)$/;
    my ($a1,$a2)=($1,$2);
    $go2=~/GO:((.....).+)$/;
    my ($b1,$b2)=($1,$2);
    unless(defined($prob{$gopair}{"high"})){
	if (-e "$network_file.prob"){
	    open(A,"$network_file.prob");
	    while(<A>){
		chomp;
		if(( /$b1/) &&( /$a1/)){
		    my @tt=split(/\t/);
		    $prob{$gopair2}{"high"}=$prob{$gopair}{"high"}=$tt[2];
		    $prob{$gopair2}{"low"}=$prob{$gopair}{"low"}=$tt[3];
		}
	    }
	    close(A);
	}
    }
    unless(defined($prob{$gopair}{"high"})){
	## I have calculated the probabilities of ALL possible
	## GO pairs and split them into files according to the GO
	## number. So, GO:0019952 will be in the file 00199
	my $file;
	$a1>$b1 ? ($file=$b2 . ".prob") : ($file=$a2 . ".prob")  ;
# #    If we have already seen this gopair, skip
	unless (defined($prob{$gopair}{"high"}) || defined($prob{$gopair2}{"high"})){
	    if(-e "$stats_dir/$file.gz") {
		open(A,"zcat $stats_dir/$file |")|| die("cannot open $stats_dir/$file : $!\n$a1,$a2:$b1,$b2,$go1,$go2,$gopair\n");
	    }
	    elsif(-e "$stats_dir/$file") {
		open(A,"$stats_dir/$file")|| die("cannot open $stats_dir/$file : $!\n$a1,$a2:$b1,$b2,$go1,$go2,$gopair\n");
	    }
	    else{die("bad stats filname: $stats_dir/$file:$1\n")}
	    while(<A>){
		chomp;
		if(( /$b1/) &&( /$a1/)){
		    my @tt=split(/\t/);
		    $prob{$gopair2}{"high"}=$prob{$gopair}{"high"}=$tt[2];
		    $prob{$gopair2}{"low"}=$prob{$gopair}{"low"}=$tt[3];
		    last;
		}
	    }
	    close(A);
	}
## Create probalities file for this network,so that we can run faster next time
## it is printed at the end of main  0051301

	
    }
    ## if the second parameter is 1 return the high tail probability
    ## else,  return the low tail probability
    unless( ((defined($prob{$gopair}{"high"})) && defined($prob{$gopair}{"low"})) ||
	   (defined($prob{$gopair2}{"high"}) && defined($prob{$gopair2}{"low"})) ) {
		die("crap : $gopair/$gopair2: " . $prob{$gopair}{"high"} . ":". $prob{$gopair2}{"high"} ."\n") ;
		print STDERR "Probability problem, missing $gopair\n";
	$prob{$gopair}{"high"}=$prob{$gopair}{"low"}=-1;
	$prob{$gopair2}{"high"}=$prob{$gopair2}{"low"}=-1;
#	print "$gopair\n"; die();
    }

    push @net_probs, "$go1\t$go2\t" . $prob{$gopair}{"high"} . "\t" . $prob{$gopair}{"low"} . "\n";
    #print "NN : @net_probs\n";

    

    $tail == 1 ? return($prob{$gopair}{"high"}) : return($prob{$gopair}{"low"});
}

############################################################

sub usage{   
    open(HELP, "| more") ;
    print HELP <<EndOfHelp;

USAGE:  
	moonlighting.pl [options] <ANNOTATIONS FILE> <NETWORK FILE>

moonlighting.pl will take a class annotation file and a network and
look for interacting proteins whose GO annotations are very disimilar.
It will return a list of proteins with interaction probabilities 
below a given threshold. It also requires certain files that are
described at the end of this help.

COMMAND-LINE OPTIONS:

    -d : Debugging output, VERY verbose
    -D : Data directory, default is current directory.
    -f : Gene Ontology GAF format annotation file (def ./gene_association.goa_human)
    -g : geneology file, output of OntoGenealogy.pl (def ./biological_process.genealogy)
    -h : Print this help and exit
    -G : GO terms file (def ./GO.terms_ids_obs)
    -o : Gene Ontology:
           p : biological process (default)
	   c : cellular compartment
	   f : biological function
    -p : minimum probability threshold (def= 0.0005)
    -P : maximum probability threshold (def= 0.5)
    -s : Synonyms file (def ./synonyms_human)
    -S : Species, (def: HUMAN) this is useful to avoid counting names
         such as CBL and CBL_HUMAN twice. 
    -v : Verbose output
     
FILES:

There are a number of files required to run moonlihting.pl. Most can 
be specified by command line options. The GO association probabilities 
cannot. The program expects to find a subdirectory ./gostats containing
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
	  P#	2
	  PN    CCDC6_HUMAN, DAPK1_HUMAN


    Gene Ontology GAF format annotation file:
         This is the standard format GAF file as downloaded from 
	 geneontology.org. For an explanation of the expected format, 
	 see http://www.geneontology.org/GO.format.gaf-2_0.shtml

    Geneology File:
	 This is the output of OntoGenealogy.pl, each line is divided 
	 into tab-separated columns. The first column has the parent 
	 GO term and each of the other column its children:
	 
	 parent          child1          child2
	 GO:0071555	GO:0071554	GO:0008150

    GO terms file:
	 This is a tab separated list of GO term names, their
	 corresponding GO numbers and associated ontology:
	 downloaded from http://www.geneontology.org/doc/GO.terms_and_ids

	 GO:0000001	mitochondrion inheritance	P	

    Network File:
	 The network to be analysed. One interacting pair per line 
         and, optionally, the number of interactions in the first line.
         eg:
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

############################################################
sub terms_to_IDs{
    my $term=shift;
    $term=~s/_/ /g;
    if($have_already_read_terms_file==0){
	open(T,"$go_terms_file")|| die("Cannot open terms file : $!\n");
	while(my $line=<T>){
	    next if $line=~/^\!/; 
	    chomp;
	    my @a=split(/\t/, $line);
	    $terms_to_GOs{$a[1]}=$a[0];
	}
	close(T);
	$have_already_read_terms_file=1;
    }
    die("shit $term\n") unless defined($terms_to_GOs{$term});
#    &debug("term : $term, id:$terms_to_GOs{$term}" );
    return($terms_to_GOs{$term});
    
}
############################################################
sub debug{
    if ($debug)
    {
	print STDERR "@_\n";
    }
}


#################################################################################
############################# PROGRAM END #######################################
#################################################################################


############################ Obsolete Subroutines ############################


############################################################
# sub share_annotations{
#     my $a=shift;
#     print "a : $a\n"; die();
#     print "$gos{$_}"; die();

# }



############################################################
# sub is_related{
#     my $is_related=undef;
#     if(exists($offspring{$_[0]}{$_[1]}) || exists($parents{$_[1]}{$_[0]})
#        || exists($offspring{$_[1]}{$_[0]}) || exists($parents{$_[0]}{$_[1]})
# 	){
# 	$is_related=1;
#     }
#     return($is_related);
# }



############################################################
# sub load_ancestors1(){
#     print STDERR "Loading ancestors...\n" if $verbose;
#     if( -e $ancestor_file){
# 	open(A,"$ancestor_file")||die();;
#     }
#     else{	
# 	open(A,"echo \"select * from ancestors\" | psql -U cchapple -h 10.1.1.53  -q -d obo |")||die("db2 problem\n");
#     }
#     my $a=0;
#     while(<A>){
# 	$a++;
# 	next unless $a>2;
# 	next unless /.+\|.+\|/;
# 	next if /\d+\s*rows/;
# 	next if /^\s*$/;
# 	my @a=split;
# 	#print "cc $a[0] : $a[2] :: $foundGOs{$a[0]} $foundGOs{$a[2]}\n";
# #	push @{$parents{$a[0]}},$a[2];
# #	print "$a[2]:$foundGOs{$a[2]}\n" if $a[2]=~/GO:0050794/;	
# # if (($a[2]=~/GO:0050794/) && ($a[0]=~/GO:0051436/)){
# # 		print "ccc $a[0] $a[2] $children{$a[0]}{$a[2]} $parents{$a[2]}{$a[0]}\n";

# # 	    }
# 	if ($a[2]=~/GO:0050794/){
# 	    print "ccc $a[0] $a[2] ::  $offspring{$a[0]}{$a[2]} :  $parents{$a[2]}{$a[0]}\n";
# 	}
# 	if ($a[2]=~/GO:0050794/){
# 	    print "ccc $a[0] $a[2] ::  $offspring{$a[0]}{$a[2]} :  $parents{$a[2]}{$a[0]}\n";
# 	}
# 	next if $a[0] eq $a[2];
# 	if(exists($foundGOs{$a[0]}) && exists($foundGOs{$a[2]})) {
# 	    $offspring{$a[0]}{$a[2]}=1;
# 	    $parents{$a[2]}{$a[0]}=1;
# 	}
#     }
#     close(A);
# }





############################################################

# sub probability()
# { 
#     print STDERR "Probability...\n" if $verbose;
#     my %GOpairs;
#     my $GOpair;
#     if( -e $stats_file){
# 	if($stats_file =~ /\.gz$/) {
# 	    open(A,"zcat $stats_file |")||die("Could not open $stats_file :$!\n");
# 	}
# 	else{
# 	    open(A," $stats_file")||die("Could not open $stats_file :$!\n");
# 	}
#     }
    
#     while(<A>){
# 	my @tt=();
# 	print if /GO:0048609/ && /GO:0000738/;
# 	print STDERR "." if $. %10000 == 0 && $verbose;
# 	print "[$.]\n"  if $. %500000 == 0 && $verbose;
# ## Is this the raw database dump of temp_asso
# ## or is it the output of my read_temp_asso.pl?
# 	if($.==1 && /gene/){
# 	    close(A);
# 	    die("please run read_temp_asso.pl on the $stats_file file\n");
# 	}
# 	next if $.==1;
# 	chomp;
# 	@tt=split(/\t/);
# 	next unless defined($foundGOs{$tt[0]});
# 	next unless defined($foundGOs{$tt[1]});
# 	$prob{$tt[0] . "xxx" . $tt[1]}{HIGH}=$prob{$tt[1] . "xxx" . $tt[0]}{HIGH}=$tt[1];
# 	$prob{$tt[0] . "xxx" . $tt[1]}{LOW}=$prob{$tt[1] . "xxx" . $tt[0]}{LOW}=$tt[2];
#     }
# print "[$.]\n"  if $verbose;
#     # load gene names from temp_asso
#     # 	foreach pair of gos in temp_asso
#     # 	      nb_asso is a list of all GOs and the number of genes annotated to each of them (implicit and explicit)
#     close(A);

# }



############################################################
# sub load_associations(){
#     print STDERR "Loading associations...\n" if $verbose;

#     my @a=keys(%interacting);
#     map{
# 	my ($b,$t)=split(/xxx/);
#     }@a;

#     my $table="associations_" . $species . "_go";
#     if(-e $associations_file){
# 	if ($associations_file =~ /\.gz$/) {open(A,"zcat $associations_file |")||die("Cannot open associations, $associations_file: $!\n");}
# 	else {open(A,"$associations_file")||die("Cannot open associations, $associations_file: $!\n");}
#     }
#     else{
# 	open(A,"echo \"select * from $table\" | nice -10 psql -U cchapple -h 10.1.1.53  -q -d obo |")||die("Could not open $table: $!\n");
#     }
#     my $a=0;
#     while(<A>){
# 	s/^\s+//g;
# 	s/\s+\|\s+/\|/g;
# 	my @fields=split(/\|/);
# 	$associations{GENE_NUM}=$fields[3] unless defined($associations{GENE_NUM});
# #	print "ff :$fields[0] $fields[1] $fields[2] $fields[3] $fields[4] $fields[5]\n";
# 	$a++;
# 	$verbose && print STDERR "." if $a % 10000 == 0; 
# 	$verbose && print STDERR "\t[$a]\n" if $a % 500000 ==0; 


# 	next unless $a>2;
# 	next if /\d+\s*rows/;
# 	next if /^\s*$/;
# 	my @a=split;
# #	print "$a[0] : $a[2]\n";
# #	push @{$parents{$a[0]}},$a[2];
# #	my $n=$fields[0] . "xxx" . $fields[1];
# 	$associations{$fields[0]}{NUM}=$fields[4]; # nmbr of genes with go1
# 	$associations{$fields[1]}{NUM}=$fields[5]; # nmbr of genes with go2
# 	$associations{$fields[0]}{$fields[1]}{NUM}=$fields[6]; #nmbr of genes with both gos
# 	$associations{$fields[0]}{$fields[1]}{CMNANC}=$fields[6]; # common ancestry measure

# 	## Inherit implied associations
	
# 	die("$fields[2] $fields[1] defined\n") if defined($associations{$fields[2]}{$fields[1]});
#     }

#     close(A);
# }




############################################################
# sub get_statistics(){
#     print STDERR "Loading Stats stats...\n" if $verbose;
#     my @prot_names=keys(%known_function);
#     foreach my $prot (@prot_names){
# 	my @GOs=keys(%{$known_function{$prot}});
# 	for (my $n=0;$n<scalar(@GOs);$n++)
# 	{
# 	    for (my $k=$n; $k<scalar(@GOs);$k++){
# 		next if $GOs[$n] eq $GOs[$k];
# 		my $go1=$GOs[$n];
# 		my $go2=$GOs[$k];
# 		my $GOpair=$go1 . "xxx" . $go2;
# 		if(defined($stats{$go1 . "xxx" . $go2})){
# 		    $GOpair=$go1 . "xxx" . $go2;
# 		}
# 		elsif(defined($stats{$go2 . "xxx" . $go1})){
# 		    $GOpair=$go2 . "xxx" . $go1;
# 		}
# 		else{
# 		    $stats{$go1 . "xxx" . $go2}=$stats{$go2 . "xxx" . $go1}=0;
# 		    $GOpair=$go1 . "xxx" . $go2;	
# 		}
# 		$stats{$GOpair}++;
# 	    }
# 	}
#     }
# }


# ############################################################
# sub load_offspring_R(){
#     print STDERR "Loading offspring...\n" if $verbose;
#     if($fetch_offspring){
# 	my ($c1,$count)=0;
# 	my $tot=scalar(keys(%foundGOs));
# 	my $rscript="./.rscript";
# 	unlink($rscript);
# 	open(KK,">$$.gos");
# 	open(NN,">$$.names");
# 	map{ 
# 	    print KK "$_\n"; 
# 	    s/:/_/; 
# 	    print NN "$cwd/offspring_files/$_\n";
# 	}keys(%foundGOs);
# 	close(KK);
# 	close(NN);
# 	open (R,">>$rscript");
# 	print R "gos=scan(\"$cwd/$$.gos\",what=\"character\")\n";
# 	print R "outnames=scan(\"$cwd/$$.names\",what=\"character\")\n";
# 	print R "library(\"GO.db\")\n";
# 	open(RR,">.rscript_loop");
# 	print RR "offspring <- get(gos[i], GOBPOFFSPRING)\nout <- c(gos[i],offspring)\nwrite(out,file=outnames[i],ncolumns=length(out),append=FALSE,sep=\"\t\")\n";
# 	close(RR);
# 	print R "for (i in 1:length(gos)){try(source(\"$cwd/.rscript_loop\"))}\n";
# 	close(R);
	
# 	system("R CMD BATCH $rscript");
#     }
#     unless(-e $offspring_file){system("cat ./offspring_files/* > $offspring_file");}
#     open(OFFSPRING,"$offspring_file")||die("Cannot open $offspring_file : $!\n");
#     while(<OFFSPRING>){
# 	my @a=split(/\t/);
# 	my $dad=shift(@a);
# 	map{
# 	    $offspring{$dad}{$_}++;
# 	    $parents{$_}{$dad}++;
# 	}@a;
#     }
    
# }



########### TODO
# 1. Use precision
# 2. change style, forach possible prot pair, go through ALL their GOs and discard any whose association does not pass threshold
## Include frequency of each gopair in the dataset in prog output 



####### FILES
#read_temp_asso.pl -v temp_asso.gz> stats


###### Ancestor stuff R
#library(FactoMineR)
#library(GO.db)
# get("GO:0050794", GOBPOFFSPRING)
# xx <-as.list(GOBPOFFSPRING)
# xx <- xx[!is.na(xx)]
# goids <- xx[1:length(xx)]
# for(i in 1:length(xx)){
# cc <- c(goids[i][1])
# write.infile(cc,file="aa",sep=";", append = TRUE)
# }

