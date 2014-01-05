#!/usr/bin/perl -w
use strict;
use Getopt::Std;
use IO::File;
use Cwd; 

my (%stats,%known_function,%go_annotations,%ancestors,%candidates,%prob,%skip_this,%missing,%singletons,%terms_to_GOs,%foundGOs,%proteins,%synonyms,%opts,%gos,%interactions,%found,%parents,%offspring, %associations,%interacting);
#my (@interactions);
$synonyms{LOADED}=0;
getopts('gdvhSup:P:G::o:s:a:A:t:f:',\%opts) || do { print "Invalid option, try 'moonlighting.pl -h' for more information\n"; exit(1); };
my $have_already_read_terms_file=0;
my $annotation_file_type;
my $gaf_annotations_file=$opts{f} || "gene_association.goa_human";
&usage() if $opts{h};

my $stats_dir="/home/terdon/research/testing/gostats";
my $geneology_file="biological_process.genealogy";
my $cwd=cwd();
my $print_gos=$opts{g}|| undef;
my $fetch_offspring=$opts{c}||undef;
my $debug=$opts{d}||undef;
my $help=$opts{h}||undef;
my $synfile=$opts{s}||"synonyms_human";
my $stats_file=$opts{t}||"hyper";
my $ancestor_file=$opts{a}||"ancestors";
my $associations_file=$opts{A}||"associations.gz";
my $annotations_file=$ARGV[0]|| die("Need a annotations file\n");
my $network_file=$ARGV[1]|| die("Need a network file\n");
my $offspring_file=$opts{o}||$network_file . ".gokids";
my $low_prob=$opts{p}||0.0005;
my $high_prob=$opts{P}||0.5;
my $subonto=$opts{O}||"P";
my $print_unique_interactions=$opts{u}|| undef;
my $silent=$opts{S}|| undef;
my $species=$opts{s}|| "human";
my $verbose=$opts{v}||undef;
my $candidates=0;
my $go_terms_file=$opts{G}||"GO.terms_ids_obs";
if($subonto eq "p"){$subonto="P"}
elsif($subonto eq "c"){$subonto="P"}
elsif($subonto eq "f"){$subonto="P"}
elsif($subonto eq "F" || "P" || "C"){}
else{&usage("Subontology (-o) must be one of \"p,c,f\"");}
$verbose=1 if $debug;
print STDERR "#################################################################\n" if $verbose;
print STDERR "#################################################################\n" if $verbose;
print STDERR "#################################################################\n" if $verbose;
print STDERR "#################################################################\n" if $verbose;

open(E,">er"); ## logfile

## Read the network file
&load_network($network_file);
## Load functional annotations
&load_annotations($annotations_file);
&load_ancestors();
&parse_gaf_file();

## Do your thang!
&main();

close(E);
############################ SUBROUTINES ############################ 

sub main{
    print STDERR "Starting main...\n" if $verbose;
    my $koko=0;
    my $pp=scalar(keys(%interactions));
  bait:foreach my $bait (keys(%interactions)){
      next bait if defined($skip_this{$bait});
      my (@bait_gos,@target_gos,@target_names);
      my(%bait_hash,%bait_go_count,%target_go_count)=();
      $koko++;      
      printf(STDERR "$koko of $pp\r");
      ##get all gos for bait protein
      hash:foreach my $h (@{$proteins{$bait}{$subonto}}){
	  if ($h eq "none"){
	      push @bait_gos,"none";
	  }
	  else{
	      %bait_hash=%{$h};
	      &debug("baits go $bait",$bait_hash{"goID"});
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
      hash:foreach my $h (@{$proteins{$target}{$subonto}}){
	  if ($h eq "none"){
	      push @target_gos,"none";
	      push @target_names,$target;
	  }
	  else{	
	      my %target_hash=%{$h};
	      ## If the bait protein has even one of the target's GOs, 
	      ## skip target
	      if (exists($bait_go_count{$bait}{$target_hash{"goID"}})){
		  &debug("skipped $target because $bait shares $target_hash{goID}");
		  next target ;
		  }
	      foreach my $baitgo (@bait_gos){	 
	      ## If one of the bait protein's GOs is 
	      ## related to one of the targets, skip target
		  if (exists($offspring{$baitgo}{$target_hash{"goID"}}) || 
		      exists($offspring{$target_hash{"goID"}}{$baitgo})){
		      &debug("$baitgo is related to $target_hash{goID}");
		      next target ;
		  }
	      ## Else, keep bait/target pair as a candidate
	      ## if the GO pair passes the thresholds  
		  else{
		      my $GOpair=$baitgo . "xxx" .$target_hash{"goID"};
		      ## Skip if either protein has no GO
		      next if $GOpair =~/none/;
		      my $P=&gopair_probability($GOpair,0)|| print STDERR "Probability problem, missing $GOpair\n";
		      if(($P<$high_prob) && $P>$low_prob){
			 $candidates{$bait}{$target}=$GOpair;
		      }
		      else{
			  next target;
		      }
		  }
	      }
	  }## end else
      }## end foreach target_hash 	
    }## end forach $target
  }## end foreach bait

    ## Now, print those bait target pairs that pass the thresholds
    foreach my $bait (keys(%candidates)){
	foreach my $target (keys(%{$candidates{$bait}})){
	    my @a=split(/xxx/,$candidates{$bait}{$target});
	    print "$bait ($a[0])\t$target($a[1])\n";
	}	
    }
    
}## end main
############################################################

sub load_network{
    print STDERR "Loading Network...\n" if $verbose;
     #open($A,"$network_file")|| die("Cannot open $_[0]:$!\n");
    my $A = IO::File->new("< $network_file")|| die("Cannot open $network_file : $!\n"); # even better!
    my $n=0;
    while(<$A>){
	next if /^\d+$/;
	&debug("============nwk line==============\n $_");
	my $line=$_; #$_ is empty after &get_name
	my ($bait,$target)=split(/\s+/,$_);
	
	$proteins{$bait}{F}[0]=$proteins{$bait}{P}[0]=$proteins{$bait}{C}[0]="none";
	$proteins{$target}{F}[0]=$proteins{$target}{P}[0]=$proteins{$target}{C}[0]="none";
	&debug("axa \$proteins{$bait}{F}[0] : $proteins{$bait}{F}[0]=$proteins{$bait}{P}[0]=$proteins{$bait}{C}[0]\n\$proteins{$target}{F}[0] : $proteins{$target}{F}[0]=$proteins{$target}{P}[0]=$proteins{$target}{C}[0]");
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
	    my $sillyvar="!";
	    ## silly thing because emacs screws up the syntax 
	    ## highlihting when /!/;
	    next if /^$sillyvar/;
	    chomp;
	    my @tmparray=split(/\t/);
	    if (scalar(@tmparray)!=18){die("Bad format annotation file\n");};
	    die("\n\nBadly formatted GAF line:\n$_\n\nGAF file must be edited with the desired protein names, ie those in the network, as the first field on each line\ne.g.:\nSMRC2_HUMAN	UniProtKB	Q8TAQ2	SMARCC2	NOT	GO:0005730	PMID:18029348	IDA		C	SWI/SNF complex subunit SMARCC2	SMARCC2|BAF170|IPI00216047|IPI00150057|SMRC2_HUMAN|Q92923|Q96E12|Q96GY4	protein	taxon:9606	20090616	HPA\n\nTry again...\n") unless scalar(@tmparray)==18;
	    next unless $tmparray[9] eq $subonto;
	    if(defined($found{$tmparray[0]})){
		&debug("Adding(1): @tmparray\n");
		&add_protein(@tmparray);
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
			&debug( "ff : \$foundGOs{$T}: $foundGOs{$T}");
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
			&debug("$protein also has $GOs[$n]");
			$hash{"Qualifier"}="is";
			${$proteins{$protein}{$subonto}}[$n]=\%hash;
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
	     &debug("111 \$proteins{$name}{$subonto}[0] :  $proteins{$name}{$subonto}[0]" );
	    if($proteins{$name}{$subonto}[0] eq "none"){
		for(my $n=0;$n<scalar(@{$proteins{$name}{$subonto}});$n++){
		    &debug("MISSING $n : $name : $proteins{$name}{$subonto}[$n]");
		}
		$missing{$name}++;
	    }
	}
    }
}
############################################################
sub parse_gaf_file{
    print STDERR "Getting missing annotations...\n" if $verbose;
    open(A,"$gaf_annotations_file")|| die("Cannot open $gaf_annotations_file:$!\n");
    while(<A>){
	my $sillyvar="!";
	next if /^$sillyvar/; ## silly thing cause emacs screws up the syntax highlihting when /!/;
	chomp;
	my @tmparray=split(/\t/);
	next unless $tmparray[8] eq $subonto;
	my $nn=&get_name($tmparray[2],"missing1");
	die("nn : $nn\n") if $nn =~ /2DMB_HUMAN/;
	if(defined($missing{$nn})){
	    &debug("Adding(2) $nn : @tmparray\n");
	    &add_protein($nn,@tmparray);
	}
	$known_function{$nn}{$tmparray[4]}++;
	map{$known_function{$nn}{$_}++}@{$ancestors{$tmparray[4]}};
	
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
    }
    close(ANC);
}


####################################################################### 


sub add_protein{
#    &debug("\nAdding : @_\n");
    my $name=shift;
    $name=&get_name($name,"add_protein for");
    my($db) = shift;
    my($dbID) = shift;
    my($dbSymbol) = shift;
    my($Qualifier) = shift;
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
	"aspect"	=> $aspect||"null",
	"dbObjName"	=> $dbObjName||"null",
	"Synonym"	=> $synonym||"null",
	"type"   	=> $type||"null",
	"taxon"		=> $taxon||"null",
	"Date"		=> $Date||"null",
	"AssBy"		=> $AssBy||"null",
	"AnnotExt"	=> $AnnotExt||"null",
	"GPFormId"	=> $GPFormId||"null"
        );
    die("xaxa no $name\n")  unless exists($proteins{$name}{$aspect}[0]);
    if ($proteins{$name}{$aspect}[0] eq "none"){
	shift(@{$proteins{$name}{$aspect}}); 
    }
    push @{$proteins{$name}{$aspect}},\%hash;
    map{&debug("$_ : -$hash{$_}-")}keys(%hash);
}
############################################################
sub get_name{
    my $name=$_[0];
    my $old_name=$name;
    my $called_by=$_[1]||"null";
    my $reverse=$_[2]||undef;
    &debug("Called by : $called_by for $name ($synfile)");
    if ($synonyms{LOADED}==0){
	if(-e $synfile){
	    open(S,$synfile);
	    while(<S>){
		/^(.+)\t(.+)/ || die("Each line of the synonyms file should contain a gene name (left) and its desired synonyms (right), separated by a tab.\n");
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
	    close(S);
	}
	else{
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
	die("shiiit1\n") if $synonyms{GOOD}{$hname} >1;
    }
    elsif(defined($synonyms{$name})){
	$name=$synonyms{$name} 
    }
    elsif ($name=~/(.+)_$species/i && defined($synonyms{$1})){
	$name=$synonyms{$1};
    }
    else{$synonyms{$name}=$name;}
    &debug("NNNNNN : $old_name => $name");
    return $name;
    
}
############################################################
sub gopair_probability{
    #  $gopair e.g. : GO:0000001xxxGO:0050876
    my $gopair=shift;
    my $tail=shift;
    my ($high,$low);
    my ($go1,$go2)=split(/xxx/,$gopair);
    
    ## I have calculated the probabilities of ALL possible
    ## GO pairs and split them into files according to the GO
    ## number. So, GO:0019952 will be in the file 00199
    $go1=~/GO:((.....).+)$/;
    my ($a1,$a2)=($1,$2);
    $go2=~/GO:((.....).+)$/;
    my ($b1,$b2)=($1,$2);
    my $file;
    $a1>$b1 ? ($file=$b2 . ".prob.uniq") : ($file=$a2 . ".prob.uniq")  ;
#    If we have already seen this gopair, skip
    unless (defined($prob{$gopair}{"high"})){
	open(A,"$stats_dir/$file")|| die("cannot open $stats_dir/$file : $!\n$a1,$a2:$b1,$b2,$go1,$go2,$gopair\n");
	while(<A>){
	    if(( /$b1/) &&( /$a1/)){
		my @tt=split(/\t/);
		$prob{$gopair}{"high"}=$tt[1];
		$prob{$gopair}{"low"}=$tt[2];
		last;
	    }
	}
    }
    ## if the second parameter is 1 return the high tail probability
    ## else,  return the low tail probability
    $tail == 1 ? return($prob{$gopair}{"high"}) : return($prob{$gopair}{"low"});
}

############################################################

sub usage{
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
    &debug("term : $term, id:$terms_to_GOs{$term}" );
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

