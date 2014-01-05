#!/usr/bin/perl -w
use strict;
use Getopt::Std;

my (%singletons,%terms_to_GOs,%foundGOs,%stats,%proteins,%synonyms,%opts,%gos,%interactions,%found,%parents,%children, %associations,%interacting);
#my (@interactions);
$synonyms{LOADED}=0;
getopts('gdvhSuG::o:s:a:A:t:',\%opts) || do { print "Invalid option, try 'moonlighting.pl -h' for more information\n"; exit(1); };
my $have_already_read_terms_file=0;
my $annotation_file_type;
&usage() if $opts{h};
my $print_gos=$opts{g}|| undef;
my $debug=$opts{d}||undef;
my $help=$opts{h}||undef;
my $synfile=$opts{s}||"synonyms";
my $temp_asso=$opts{t}||"temp_asso.gz";
my $ancestor_file=$opts{a}||"ancestors";
my $associations_file=$opts{A}||"associations.gz";
my $annotations_file=$ARGV[0]|| die("Need a annotations file\n");
my $network_file=$ARGV[1]|| die("Need a network file\n");
my $subonto=$opts{o}||"P";
my $print_unique_interactions=$opts{u}|| undef;
my $silent=$opts{S}|| undef;
my $species=$opts{s}|| "human";
my $verbose=$opts{v}||undef;
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

 &load_network($network_file);
 &load_annotations($annotations_file);
 &load_ancestors();
#&load_associations();
print STDERR "Starting main...\n" if $verbose;
&main();

############################ SUBROUTINES ############################ 

sub main{
#for(my $n=0;$n<scalar(@interactions);$n++){
    foreach my $bait (keys(%interactions)){
	my @aa; ## @aa holds target names
	my (@bgos,@tgos);
	my (%b,%t);
	foreach my $h (@{$proteins{$bait}{$subonto}}){
	    &debug("H1: \@{\$proteins{$bait}{$subonto}} : @{$proteins{$bait}{$subonto}}\nH2: $h");
	    if ($h eq "none"){
		push @bgos,"none";
	    }
	    else{
		my %hash=%{$h};
		$b{$hash{"goID"}}++;
		next if $hash{"Qualifier"} eq "NOT";
		## Only keep most specific GO term for each lineage
		push @bgos,$hash{"goID"} unless $b{$hash{"goID"}}>1;
		# foreach my $bgo (@bgos){
		#     print STDOUT "xxx $bgo, " . $hash{"goID"} . "\n";
		#     exists($parents{$bgo}{$hash{"goID"}}) ? print "$bgo is a parent of " .  $hash{"goID"} . "\n"  && shift(@bgos): print "";
		#     exists($children{$bgo}{$hash{"goID"}}) ? print "$bgo is a child of " .  $hash{"goID"} . "\n"  : print "";
		# }
	    }
	}
	
      target:foreach my $target (@{$interactions{$bait}}){
	  @tgos=();
	  foreach my $h (@{$proteins{$target}{$subonto}}){
	      if ($h eq "none"){
		  push @tgos,"none";
		  push @aa,$target;
	      }
	      else{
		  my %hash=%{$h};
		  ## If the bait protein has even one of the target's GOs, 
		  ## skip target
		  next target if exists($b{$hash{"goID"}});

		  ## If one of the bait protein's GOs is 
		  ## related to one of the targets, skip target
		  map{next target if &is_related($_,$hash{"goID"})}@bgos;

		  $t{$hash{"goID"}}++;
		  unless ($t{$hash{"goID"}}>1){
		      push @tgos,$hash{"goID"};
		      push @aa,$target;
		  }
	      }
	  }
	  map{
	      for (my $n=0;$n<scalar(@tgos);$n++){
		  #$a : eg GO:0045184xxxGO:0045762
		  my $a=$_."xxx".$tgos[$n];
		  $gos{$a}{NUM}++;
		  #$bb : eg FLNAxxxCALCR
		  my $bb=$bait . "xxx" . $aa[$n];
		  $gos{$a}{PROT}=$bb;
	      }
	  }@bgos;
      }
    }

## Assess if specific interactions are statistically significant
    #&probability();

    # ## Print go terms that interact in just ONE targat/bait pair
    unless($silent){

	## if this is a class file
	if($annotation_file_type==1)
	{
	    foreach my $singleton(keys(%singletons)){
	#	print "$singleton is a singleton\n";
	    }
	    
	}
	map{
	    unless(/none/){
		if($gos{$_}{NUM}<2){
		    my @bb=split(/xxx/,$gos{$_}{PROT});
		    my @gg=split(/xxx/);
#		unless($parents{$bb[0]}{$bb[1]} || $children{$bb[0]}{$bb[1]}){
		    unless(&is_related($gg[0],$gg[1])){
			#   unless(&share_annotations($bb[0],$bb[1])){
			print "$bb[0]($gg[0])\t\t$bb[1]($gg[1])\n";
#		    print "$bb[0]\n$bb[1]\n";
			#}
			
		    }
		}
	    }
	     
	}keys(%gos);
    }
}
############################################################

sub load_network{
    print STDERR "Loading Network...\n" if $verbose;
     open(A,"$network_file")|| die("Cannot open $_[0]:$!\n");
     my $n=0;
    while(<A>){
	next if /^\d+$/;
	&debug("============nwk line==============\n $_");
	my $line=$_; #$_ is empty after &get_name
	my ($b,$t)=split(/\s+/,$_);
	my $bait=(&get_name($b,"bait " . $line));
	my $target=(&get_name($t,"target " . $line));
	
	$proteins{$bait}{F}[0]=$proteins{$bait}{P}[0]=$proteins{$bait}{C}[0]="none";
	$proteins{$target}{F}[0]=$proteins{$target}{P}[0]=$proteins{$target}{C}[0]="none";
		 &debug("bb : $bait : $target");
	&debug("axa $proteins{$bait}{F}[0]=$proteins{$bait}{P}[0]=$proteins{$bait}{C}[0]") if $b =~ "/PROM/";
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
	   
	## If this is a GO GAF file
	if($annotation_file_type==0){
	    my $silivar="!";
	    next if /^$silivar/; ## silly thing cause emacs screws up the syntax highlihting when /!/;
	    chomp;
	    my @tmparray=split(/\t/);
	    #if(defined($proteins{&get_name($tmparray[2],"one")})){
	    next unless $tmparray[8] eq $subonto;
	    if(defined($found{&get_name($tmparray[2],"one")})){
		&add_protein(@tmparray);
#		print  "a :  \@{\$proteins{&get_name($tmparray[2])}{$subonto}} @{$proteins{&get_name($tmparray[2])}{$subonto}} \n-" . &get_name($tmparray[2]) . "-\n";##malaka
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
			push @GOs,&terms_to_IDs($term);
		    }
		}
	    }
	    if(/^P\#\s+(\d+)/){
		$singleton=1 if $1 == 1;
	    }
	    if(/^PN\s*(.+)/){
		my $kk=$1;
		if($singleton==1){
		    $singletons{&get_name($kk,"singleton")}=1;
		}
		my @prots=split(/,\s+/,$kk);
		my %hash;
		foreach my $goterm(@GOs){
		    map{
			$hash{"goID"}=$goterm;
			$hash{"Qualifier"}="is";
			push @{$proteins{$_}{$subonto}},\%hash;
		    }@prots;
		}		
	    }		
	}
    }
    close(A); 
}
############################################################


sub load_ancestors(){
    print STDERR "Loading ancestors...\n " if $verbose;
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
#	print "$a[0] : $a[2]\n";
#	push @{$parents{$a[0]}},$a[2];
	next unless defined($foundGOs{$a[0]}) || defined($foundGOs{$a[2]});
	next if $a[0] eq $a[2];
	$children{$a[0]}{$a[2]}=1;
	$parents{$a[2]}{$a[0]}=1;

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
    &debug("\nAdding : @_\n");
    my($db) = shift;
    my($dbID) = shift;
    my($old_dbSymbol) = shift;
    my $dbSymbol=&get_name($old_dbSymbol,"three");
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
    map{$synonyms{$_}=$dbSymbol;}split(/\|/,$synonym);
    $synonyms{$dbID}=$dbSymbol;
    my %hash= (
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
    &debug("xaxa no $dbSymbol") unless exists($proteins{$dbSymbol}{$aspect}[0]);
    if ($proteins{$dbSymbol}{$aspect}[0] eq "none"){&debug("\@{\$proteins{$dbSymbol}{$aspect}} : @{$proteins{$dbSymbol}{$aspect}}"); shift(@{$proteins{$dbSymbol}{$aspect}}); &debug("\@{\$proteins{$dbSymbol}{$aspect}} : @{$proteins{$dbSymbol}{$aspect}}");}
    push @{$proteins{$dbSymbol}{$aspect}},\%hash;
    &debug( "3hhaa : @{$proteins{$dbSymbol}{$aspect}}::");
    map{debug("$_ : -$hash{$_}-")}keys(%hash)
}

############################################################
sub get_name{
    my $name=$_[0];
    my $old_name=$name;
    my $called_by=$_[1]||"null";
    &debug("Called by :$called_by");
    if ($synonyms{LOADED}==0){
	if(-e $synfile){
	    open(S,$synfile);
	    while(<S>){
		/^(.+)\t(.+)/ || die("Each line of the synonyms file should contain a gene name (left) and its desired synonym (right), separated by a tab.\n");
		$synonyms{$1}=$2;
	    }
	    close(S);
	}
	else{
	    print STDERR "No synonyms file found, parsing annotations... \n" if $verbose;
	    open(A,"$ARGV[0]")|| die("Cannot open $ARGV[0]:$!\n");
	    while(<A>){
		my $silivar="!";
		next if /^$silivar/; ## silly thing cause emacs scres up the syntax highlihting when /!/;
		my @tmparray=split(/\t/);
		map{$synonyms{$_}=$tmparray[2]}split(/\|/,$tmparray[10]);
	    }
	    close(A);
	}	
 	$synonyms{LOADED}=1;
    }
    if(defined($synonyms{$name})){
	$name=$synonyms{$name} 
    }
    elsif ($name=~/(.+)_$species/i && defined($synonyms{$1})){
	$name=$synonyms{$1};
    }
    else{$synonyms{$name}=$name;}
       &debug("NNNNNN : $old_name => $name");#malakadie("shit1: $synonyms{$name}") if $name =~ /PROM1/;
    return $name;
}
############################################################

sub probability()
{
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
	open(A,"echo \"select * from temp_asso\" | psql -U cchapple -h 10.1.1.53  -q -d obo |")||die("temp_asso problem\n");
	while(<A>){
	    print
	}
    }
    while(<A>){
	
	next if $.<3;
	next unless /\|/;
	print STDERR "." if $. %100000 == 0 && $verbose;
	print "[$.]\n"  if $. %10000000 == 0 && $verbose;
	s/^\s*//;
	my @tt=split(/\s+\|\s*/);
	my $GOpair=$tt[1] . "xxx" . $tt[2];
	die("crap, defined $tt[2]" . "xxx" ."$tt[1] batman\n") if defined($stats{$tt[2] . "xxx" .$tt[1] });
	die("crap, defined2 $tt[2]" . "xxx" ."$tt[1] batman\n") if defined($gos{$tt[2] . "xxx" .$tt[1] });
	die("crap, defined3 $tt[2]" . "xxx" ."$tt[1] batman\n") if defined($gos{$GOpair});

	$stats{$GOpair}++;
    }
    print STDERR "[$.]\n" if $verbose;
    
    
    # load gene names from temp_asso
    # 	foreach pair of gos in temp_asso
    # 	      nb_asso is a list of all GOs and the number of genes annotated to each of them (implicit and explicit)



    close(A);

}
############################################################

sub is_related{
    my $is_related=0;
    if(defined($children{$_[0]}{$_[1]}) || defined($parents{$_[1]}{$_[0]})){
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
	open(T,"$go_terms_file")|| die("Cannot open  : $!\n");
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
    &debug("term : $term, id :$terms_to_GOs{$term}" );
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
# 2. Use precision
# 2. Find a way to load both class annotations and GAF file for those with no class
