#!/usr/bin/perl -w
use strict;
use Getopt::Std;
use IO::File;
use Cwd; 

my (%candidates,%prob,%skip_this,%missing,%singletons,%terms_to_GOs,%foundGOs,%stats,%proteins,%synonyms,%opts,%gos,%interactions,%found,%parents,%offspring, %associations,%interacting);
#my (@interactions);
$synonyms{LOADED}=0;
getopts('gdvhSuG::o:s:a:A:t:f:',\%opts) || do { print "Invalid option, try 'moonlighting.pl -h' for more information\n"; exit(1); };
my $have_already_read_terms_file=0;
my $annotation_file_type;
my $gaf_annotations_file=$opts{f} || "gene_association.goa_human";
&usage() if $opts{h};
my $cwd=cwd();
my $print_gos=$opts{g}|| undef;
my $fetch_offspring=$opts{c}||undef;
my $debug=$opts{d}||undef;
my $help=$opts{h}||undef;
my $synfile=$opts{s}||"synonyms_human";
my $temp_asso=$opts{t}||"stats";
my $ancestor_file=$opts{a}||"ancestors";
my $associations_file=$opts{A}||"associations.gz";
my $annotations_file=$ARGV[0]|| die("Need a annotations file\n");
my $network_file=$ARGV[1]|| die("Need a network file\n");
my $offspring_file=$opts{o}||$network_file . ".gokids";
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
&load_network($network_file);
&load_annotations($annotations_file);
# &load_ancestors();
&load_offspring();
#&load_associations();
print STDERR "Starting main...\n" if $verbose;
&main();
close(E);
############################ SUBROUTINES ############################ 

sub main{
#for(my $n=0;$n<scalar(@interactions);$n++){
  bait:foreach my $bait (keys(%interactions)){
      next bait if defined($skip_this{$bait});
	my @target_names;
	my (@bgos,@tgos);
	my (%b,%t);
	foreach my $h (@{$proteins{$bait}{$subonto}}){
	    &debug("H1: \@{\$proteins{$bait}{$subonto}} : @{$proteins{$bait}{$subonto}}\nH2: $h");
	    if ($h eq "none"){
		push @bgos,"none";
	    }
	    else{
		my %hash=%{$h};
		&debug("baits go $bait",$hash{"goID"});
		$b{$hash{"goID"}}++;
		next if $hash{"Qualifier"} eq "NOT";
		## Some Gos are mentioned twice in the annotation file
		push @bgos,$hash{"goID"} unless $b{$hash{"goID"}}>1;
	    }
	}
  #   bgo:foreach my $bgo (@bgos){
# 	  next bgo if $bgo=~/none/;
# 	  foreach my $h (@{$proteins{$bait}{$subonto}}){
# 	      my %hash=%{$h};
# #	      print STDOUT "xxx $bgo, " . $hash{"goID"} . "\n";
# 	      exists($parents{$bgo}{$hash{"goID"}}) ? print STDERR "- $bgo is a parent of " .  $hash{"goID"} . "\n"  && shift(@bgos): print "";
# 	      exists($offspring{$hash{"goID"}}{$bgo}) ? print STDERR "::- $bgo is a child of " .  $hash{"goID"} . "\n"  : print "";
# 	  }
#       }
      
      target:foreach my $target (@{$interactions{$bait}}){
	  next target if defined($skip_this{$target});
	  @tgos=();
	  foreach my $h (@{$proteins{$target}{$subonto}}){
		 	   
	      if ($h eq "none"){
		  push @tgos,"none";
		  push @target_names,$target;
	      }
	      else{
		  my %hash=%{$h};		 
		  ## If the bait protein has even one of the target's GOs, 
		  ## skip target
		  next target if exists($b{$hash{"goID"}});
		  ## If one of the bait protein's GOs is 
		  ## related to one of the targets, skip target
		  map{	    
#		      print "xxx2 $_ : " . $hash{"goID"} ."\n";
		      if (exists($offspring{$_}{$hash{"goID"}}) || 
			  exists($offspring{$hash{"goID"}}{$_})){
#			  print STDERR "Skipped $_, " . $hash{"goID"} . "\n" if $verbose;
			  next target ;
		      }
		  }@bgos;

		  $t{$hash{"goID"}}++;
		  unless ($t{$hash{"goID"}}>1){
		      push @tgos,$hash{"goID"};
		      push @target_names,$target;
		  }
	      }# else, if $h != none
	  }## end foreach %h
	  map{
	      for (my $n=0;$n<scalar(@tgos);$n++){
		  #$a : eg GO:0045184xxxGO:0045762
		  my $a=$_."xxx".$tgos[$n];
		  $a=$_."xxx".$tgos[$n];
		  #$bb : eg FLNAxxxCALCR
		  my $bb=$bait . "xxx" . $target;
		  $gos{$bb}{PROT}=$bb;
		  $gos{$bb}{NUM}++;
		  push @{$gos{$bb}{GOs}},$a;;
	      }
	  }@bgos;
	  
	  if (($target=~/ATX2_HUMAN/ || $bait=~/ATX2_HUMAN/) && ($target=~/ATX1_HUMAN/ || $bait=~/ATX1_HUMAN/)){
	      my @kk=keys(%b);
	      print "target : $target, gos : @kk\n";
	      print "bait : $bait, gos : @tgos\n";
	      die("whiiiit\n");
	  }
      }## end foreach target
    } ## end foreach $bait
## Assess if specific interactions are statistically significant
    &probability();
    # ## Print go terms that interact in just ONE targat/bait pair
    unless($silent){	## if this is a class file
	if($annotation_file_type==1)
	{
	    foreach my $singleton(keys(%singletons)){
		#	print "$singleton is a singleton\n";
	    }
	}
      protpair:foreach my $protpair (keys(%gos)){
	  foreach my $gopair (@{$gos{$protpair}{GOs}}){
	      unless($gopair=~/none/){
		  if(defined($prob{$gopair}{LOW})){
		      if($prob{$gopair}{LOW}<0.0005){
			  my @bb=split(/xxx/,$gos{$protpair}{PROT});
			  my @gg=split(/xxx/,$gopair);

			  if ($bb[0]=~/ATX1_HUMAN/ && $bb[1]=~/ATX2_HUMAN/){
			      die("crap")
			      }

			  unless(&is_related($gg[0],$gg[1])){
			      print "$bb[0]($gg[0])\t\t$bb[1]($gg[1]), $prob{$gopair}{LOW}\n" unless $candidates{$bb[0]}{$bb[1]};
			     # print "@bb\n@gg\n"; 
			      $candidates{$bb[0]}{$bb[1]}++;
			      $candidates{$bb[1]}{$bb[0]}++;
			  }
		      }
		  }##defined($prob{$_}{LOW}))
		  else{
		      print E "$gopair";
		  }
	      }## unless none
	  }
      }
    } ## silent
    print STDERR scalar(keys(%candidates)) . " candidates found\n" if $verbose;
    print STDERR "####################################\n";

    foreach my $cand (keys(%candidates)){
#	print "cc : $cand\n"; 
	
    }




}## main
############################################################

sub load_network{
    print STDERR "Loading Network...\n" if $verbose;
     #open($A,"$network_file")|| die("Cannot open $_[0]:$!\n");
    my $A = IO::File->new("< $network_file")|| die("Cannot open $network_file : $!\n"); # even better!
#    my $fh=open($network_file,"<")||die("crap=\n");
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
	    my $silivar="!";
	    next if /^$silivar/; ## silly thing cause emacs screws up the syntax highlihting when /!/;
	    chomp;
	    my @tmparray=split(/\t/);
	    if (scalar(@tmparray)!=18){print; die("Bad format annotation file\n");};
	    die("\n\nBadly formatted GAF line:\n$_\n\nGAF file must be edited with the desired protein names, ie those in the network, as the first field on each line\ne.g.:\nSMRC2_HUMAN	UniProtKB	Q8TAQ2	SMARCC2	NOT	GO:0005730	PMID:18029348	IDA		C	SWI/SNF complex subunit SMARCC2	SMARCC2|BAF170|IPI00216047|IPI00150057|SMRC2_HUMAN|Q92923|Q96E12|Q96GY4	protein	taxon:9606	20090616	HPA\n\nTry again...\n") unless scalar(@tmparray)==18;
	    #if(defined($proteins{&get_name($tmparray[2],"one")})){
	    next unless $tmparray[9] eq $subonto;
	    if(defined($found{$tmparray[0]})){
		&debug("Adding(1): @tmparray\n");
		&add_protein(@tmparray);
#		print  "a :  \@{\$proteins{&get_name($tmparray[2])}{$subonto}} @{$proteins{&get_name($tmparray[2])}{$subonto}} \n-" . &get_name($tmparray[2]) . "-\n";##
	    }
	}
	## If this is a class file
	else{
	    
	    if(/^CA/){
		@GOs=();
		$singleton=0;
		if(/^CA.+:\d+\)/){die("The file has go ids\n")} ## if the file has GO ids
		
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
#		foreach my $goterm(@GOs){
		foreach my $protein(@prots){
		    for (my $n=0; $n<scalar(@GOs); $n++){
			$hash{"goID"}=$GOs[$n];
			&debug("$protein also has $GOs[$n]");
			$hash{"Qualifier"}="is";
#			the problem is here
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
# FInd those proteins that have no class annotation
	     &debug("111 \$proteins{$name}{$subonto}[0] :  $proteins{$name}{$subonto}[0]" );

	    if($proteins{$name}{$subonto}[0] eq "none"){
		for(my $n=0;$n<scalar(@{$proteins{$name}{$subonto}});$n++){
#		&debug("1x1 \@{\$proteins{$name}{$subonto}}: @{$proteins{$name}{$subonto}}");
		    &debug("MISSING $n : $name : $proteins{$name}{$subonto}[$n]");
		}
		$missing{$name}++;
	    }
	}
    }
    &get_missing_annotations();
}
############################################################

sub get_missing_annotations{
    print STDERR "Getting missing annotations...\n" if $verbose;
    open(A,"$gaf_annotations_file")|| die("Cannot open $gaf_annotations_file:$!\n");
    while(<A>){
	my $silivar="!";
	next if /^$silivar/; ## silly thing cause emacs screws up the syntax highlihting when /!/;
	chomp;
	my @tmparray=split(/\t/);
	next unless $tmparray[8] eq $subonto;
	my $nn=&get_name($tmparray[2],"missing1");
	if(defined($missing{$nn})){
	    &debug("Adding(2) $nn : @tmparray\n");
	    &add_protein($nn,@tmparray);
	}
    }
    close(A);
}



############################################################
sub load_offspring(){
    print STDERR "Loading offspring...\n" if $verbose;
    if($fetch_offspring){
	my ($c1,$count)=0;
	my $tot=scalar(keys(%foundGOs));
	my $rscript="./.rscript";
	unlink($rscript);
	open(KK,">$$.gos");
	open(NN,">$$.names");
	map{ 
	    print KK "$_\n"; 
	    s/:/_/; 
	    print NN "$cwd/offspring_files/$_\n";
	}keys(%foundGOs);
	close(KK);
	close(NN);
	open (R,">>$rscript");
	print R "gos=scan(\"$cwd/$$.gos\",what=\"character\")\n";
	print R "outnames=scan(\"$cwd/$$.names\",what=\"character\")\n";
	print R "library(\"GO.db\")\n";
	open(RR,">.rscript_loop");
	print RR "offspring <- get(gos[i], GOBPOFFSPRING)\nout <- c(gos[i],offspring)\nwrite(out,file=outnames[i],ncolumns=length(out),append=FALSE,sep=\"\t\")\n";
	close(RR);
	print R "for (i in 1:length(gos)){try(source(\"$cwd/.rscript_loop\"))}\n";
	close(R);
	
	system("R CMD BATCH $rscript");
    }
    unless(-e $offspring_file){system("cat ./offspring_files/* > $offspring_file");}
    open(OFFSPRING,"$offspring_file")||die("Cannot open $offspring_file : $!\n");
    while(<OFFSPRING>){
	my @a=split(/\t/);
	my $dad=shift(@a);
	map{
	    $offspring{$dad}{$_}++;
	    $parents{$_}{$dad}++;
	}@a;
    }
    $verbose && print STDERR "\n";
    
}
sub load_ancestors1(){
    print STDERR "Loading ancestors...\n" if $verbose;
    if( -e $ancestor_file){
	open(A,"$ancestor_file")||die();;
    }
    else{	
	open(A,"echo \"select * from ancestors\" | psql -U cchapple -h 10.1.1.53  -q -d obo |")||die("db2 problem\n");
    }
    my $a=0;
    while(<A>){
	$a++;
	next unless $a>2;
	next unless /.+\|.+\|/;
	next if /\d+\s*rows/;
	next if /^\s*$/;
	my @a=split;
	#print "cc $a[0] : $a[2] :: $foundGOs{$a[0]} $foundGOs{$a[2]}\n";
#	push @{$parents{$a[0]}},$a[2];
#	print "$a[2]:$foundGOs{$a[2]}\n" if $a[2]=~/GO:0050794/;	
# if (($a[2]=~/GO:0050794/) && ($a[0]=~/GO:0051436/)){
# 		print "ccc $a[0] $a[2] $children{$a[0]}{$a[2]} $parents{$a[2]}{$a[0]}\n";

# 	    }

	if ($a[2]=~/GO:0050794/){
	    print "ccc $a[0] $a[2] ::  $offspring{$a[0]}{$a[2]} :  $parents{$a[2]}{$a[0]}\n";
	    }
	if ($a[2]=~/GO:0050794/){
	    print "ccc $a[0] $a[2] ::  $offspring{$a[0]}{$a[2]} :  $parents{$a[2]}{$a[0]}\n";
	    }
	next if $a[0] eq $a[2];
	if(exists($foundGOs{$a[0]}) && exists($foundGOs{$a[2]})) {
	    $offspring{$a[0]}{$a[2]}=1;
	    $parents{$a[2]}{$a[0]}=1;
	    
	 
	}
    }
    close(A);
}

############################################################
sub load_associations(){
    print STDERR "Loading associations...\n" if $verbose;

    my @a=keys(%interacting);
    map{
	my ($b,$t)=split(/xxx/);
    }@a;

    my $table="associations_" . $species . "_go";
    if(-e $associations_file){
	if ($associations_file =~ /\.gz$/) {open(A,"zcat $associations_file |")||die("Cannot open associations, $associations_file: $!\n");}
	else {open(A,"$associations_file")||die("Cannot open associations, $associations_file: $!\n");}
    }
    else{
	open(A,"echo \"select * from $table\" | nice -10 psql -U cchapple -h 10.1.1.53  -q -d obo |")||die("Could not open $table: $!\n");
    }
    my $a=0;
    while(<A>){
	s/^\s+//g;
	s/\s+\|\s+/\|/g;
	my @fields=split(/\|/);
	$associations{GENE_NUM}=$fields[3] unless defined($associations{GENE_NUM});
#	print "ff :$fields[0] $fields[1] $fields[2] $fields[3] $fields[4] $fields[5]\n";
	$a++;
	$verbose && print STDERR "." if $a % 10000 == 0; 
	$verbose && print STDERR "\t[$a]\n" if $a % 500000 ==0; 


	next unless $a>2;
	next if /\d+\s*rows/;
	next if /^\s*$/;
	my @a=split;
#	print "$a[0] : $a[2]\n";
#	push @{$parents{$a[0]}},$a[2];
#	my $n=$fields[0] . "xxx" . $fields[1];
	$associations{$fields[0]}{NUM}=$fields[4]; # nmbr of genes with go1
	$associations{$fields[1]}{NUM}=$fields[5]; # nmbr of genes with go2
	$associations{$fields[0]}{$fields[1]}{NUM}=$fields[6]; #nmbr of genes with both gos
	$associations{$fields[0]}{$fields[1]}{CMNANC}=$fields[6]; # common ancestry measure

	## Inherit implied associations
	
	die("$fields[2] $fields[1] defined\n") if defined($associations{$fields[2]}{$fields[1]});
    }

    close(A);
}


############################################################

sub kk {

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
#    &debug("1hhaa :$dbSymbol :: $proteins{$dbSymbol}{$aspect}[0]:: $hash{'goID'}\n");
  #  die("name : -$dbSymbol-\n") unless exists($proteins{$dbSymbol}{$aspect}[0] );
#    &debug("xaxa no $name") unless exists($proteins{$name}{$aspect}[0]);
    die("xaxa no $name\n")  unless exists($proteins{$name}{$aspect}[0]);
    if ($proteins{$name}{$aspect}[0] eq "none"){
	#&debug("\@{\$proteins{$name}{$aspect}} : @{$proteins{$name}{$aspect}}"); 
	shift(@{$proteins{$name}{$aspect}}); 
#	&debug("\@{\$proteins{$name}{$aspect}} : @{$proteins{$name}{$aspect}}");
    }
    push @{$proteins{$name}{$aspect}},\%hash;
 #   &debug( "3hhaa : @{$proteins{$name}{$aspect}}::");
    map{debug("$_ : -$hash{$_}-")}keys(%hash)
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
#		${$synonyms{$nn}}[0]=$2;
		$synonyms{GOOD}{$nn}++;
		for (my $n=0;$n<scalar(@syns); $n++){
		    #${$synonyms{$nn}}[$n]=$syns[$n];
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

sub probability()
{
    my %GOpairs;
    my $GOpair;
    map{
	map{$GOpairs{$_}++}@{$gos{$_}{GOs}}
    }keys(%gos);
    print STDERR "Probability...\n" if $verbose;
    if( -e $temp_asso){
	if($temp_asso =~ /\.gz$/) {
	    open(A,"zcat $temp_asso |")||die("Could not open $temp_asso :$!\n");
	}
	else{
	    open(A," $temp_asso")||die("Could not open $temp_asso :$!\n");
	}
    }
    else{print STDERR "Querying database for temp_asso...\n" if $verbose;
	open(AA,"echo \"select * from temp_asso\" | psql -U cchapple -h 10.1.1.53  -q -d obo |")||die("temp_asso problem\n");
	while(<AA>){
	    print;
	}
	 close(AA);
    }
    my $database_output=0;
    while(<A>){
## Is this the raw database dump of temp_asso
## or is it the output of my read_temp_asso.pl?
	if($.==1 && /gene/){$database_output=1}
	if($database_output==1){
	    next if $.<3;
	    next unless /\|/;
	    print STDERR "." if $. %100000 == 0 && $verbose;
	    print "[$.]\n"  if $. %10000000 == 0 && $verbose;
	    s/^\s*//;
	    my @tt=split(/\s+\|\s*/);
	    my $GOpair=$tt[1] . "xxx" . $tt[2];
	    die("crap, defined $tt[2]" . "xxx" ."$tt[1] \n") if defined($stats{$tt[2] . "xxx" .$tt[1] });
	    die("crap, defined2 $tt[2]" . "xxx" ."$tt[1] \n") if defined($gos{$tt[2] . "xxx" .$tt[1] });
	    die("crap, defined3 $tt[2]" . "xxx" ."$tt[1] \n") if defined($gos{$GOpair});
	    
	    $stats{$GOpair}++;
	}
	else{
	    next if $.==1;
	    chomp;
	    my @tt=split(/\t/);
	    if(defined($GOpairs{$tt[0] . "xxx" . $tt[1]})){
		$GOpair=$tt[0] . "xxx" . $tt[1];
	    }
	    elsif(defined($GOpairs{$tt[1] . "xxx" . $tt[0]})){
		$GOpair=$tt[1] . "xxx" . $tt[0];
	    }
	    else{next}
#		print "$tt[0]-$tt[1]-\n";
	    print STDERR "." if $. %1000 == 0 && $verbose;
	    print "[$.]\n"  if $. %10000 == 0 && $verbose;
	    $prob{$GOpair}{HIGH}=$tt[2];
	    $prob{$GOpair}{LOW}=$tt[3];
	    $stats{$GOpair}++;
	    
	}
    }
    # load gene names from temp_asso
    # 	foreach pair of gos in temp_asso
    # 	      nb_asso is a list of all GOs and the number of genes annotated to each of them (implicit and explicit)
    close(A);

}
############################################################

sub is_related{
    my $is_related=undef;
    if(exists($offspring{$_[0]}{$_[1]}) || exists($parents{$_[1]}{$_[0]})
       || exists($offspring{$_[1]}{$_[0]}) || exists($parents{$_[0]}{$_[1]})
	){
	$is_related=1;
    }
    return($is_related);
}
############################################################
sub share_annotations{
    my $a=shift;
    print "a : $a\n"; die();
    print "$gos{$_}"; die();

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
sub debug
{
    if ($debug)
    {
	print STDERR "@_\n";
    }
}

########### TODO
# 1. Use precision
# 2. change style, forach possible prot pair, go through ALL their GOs and discard any whose association does not pass threshold



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

