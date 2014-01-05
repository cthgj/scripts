#!/usr/bin/perl -w
use strict;
use Getopt::Std;

my (%proteins,%synonyms,%opts,%gos,%interactions,%found,%parents,%children, %associations,%interacting);
#my (@interactions);
$synonyms{LOADED}=0;
getopts('gdvhSu:o:s:a:A:',\%opts) || do { print "Invalid option, try 'moonlighting.pl -h' for more information\n"; exit(1); };

&usage() if $opts{h};
my $print_gos=$opts{g}|| undef;
my $debug=$opts{d}||undef;
my $help=$opts{h}||undef;
my $synfile=$opts{s}||"synonyms";
my $ancestor_file=$opts{a}||undef;
my $associations_file=$opts{A}||undef;
my $annotations_file=$ARGV[0]|| die("Need a annotations file\n");
my $network_file=$ARGV[1]|| die("Need a network file\n");
my $subonto=$opts{o}||"P";
my $print_unique_interactions=$opts{u}|| undef;
my $silent=$opts{S}|| undef;
my $species=$opts{s}|| "human";
my $verbose=$opts{v}||undef;
if($subonto eq "p"){$subonto="P"}
elsif($subonto eq "c"){$subonto="P"}
elsif($subonto eq "f"){$subonto="P"}
elsif($subonto eq "F" || "P" || "C"){}
else{&usage("Subontology (-o) must be one of \"p,c,f\"");}
$verbose=1 if $debug;
print STDOUT "#################################################################\n";
print STDOUT "#################################################################\n";
print STDOUT "#################################################################\n";
print STDOUT "#################################################################\n";

&load_network($network_file);die();
&load_annotations($annotations_file);
#&load_ancestors();
#&load_associations();
print STDERR "Starting main...\n" if $verbose;
&main();

############################ SUBROUTINES ############################ 

sub main{
#for(my $n=0;$n<scalar(@interactions);$n++){
    foreach my $bait (keys(%interactions)){
	&debug("Bait: -$bait-");
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
		push @bgos,$hash{"goID"} unless $b{$hash{"goID"}}>1;
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
		  ## If the bait protein has even one of the target's GOs, skip target
		  next target if exists($b{$hash{"goID"}}); 
		  
		  $t{$hash{"goID"}}++;
		  unless ($t{$hash{"goID"}}>1){
		      push @tgos,$hash{"goID"};
		      push @aa,$target;
		  }
	      }
	      
	  }
	  map{
	      for (my $n=0;$n<scalar(@tgos);$n++){
		  my $a=$_."xxx".$tgos[$n];
		  $gos{$a}{NUM}++;
		  my $bb=$bait . "xxx" . $aa[$n];
		  $gos{$a}{PROT}=$bb;
	      }
	  }@bgos;
      }
    }

## Print go terms that interact in just ONE targat/bait pair
    if($print_unique_interactions){
	map{
	    if($gos{$_}{NUM}<2){
		my @bb=split(/xxx/,$gos{$_}{PROT});
		print "$bb[0] - $bb[1]\n";
		@bb=split(/xxx/);
		print "$bb[0] - $bb[1]\n";
	    }
	    
	    
	} keys(%gos);
    }
    unless($silent){
	map{
	    if($gos{$_}{NUM}<2){
		my @bb=split(/xxx/,$gos{$_}{PROT});
		unless($parents{$bb[0]}{$bb[1]} || $children{$bb[0]}{$bb[1]}){
		    print "$bb[0]\t\t$bb[1]\n";
		    @bb=split(/xxx/);
		    print "$bb[0]\t$bb[1]\n";
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
	
	print STDERR " \$_ is 3now $_\n";
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

}

############################################################

sub load_annotations{
    
    print STDERR "Loading annotations...\n" if $verbose;
    my $prot;
    open(A,"$annotations_file")|| die("Cannot open $_[0]:$!\n");
    my $hh=0;
    my $annotation_file_type;
    while(<A>){
	$hh++;
	if($hh==1){
	    /^\[CLASS/ ? $annotation_file_type=1 : $annotation_file_type=0;
	}
	## If this is a GO GAF file
	if($annotation_file_type==0){
	    my $silivar="!";
	    next if /^$silivar/; ## silly thing cause emacs screws up the syntax highlihting when /!/;
	    chomp;
	    my @tmparray=split(/\t/);
	    #if(defined($proteins{&get_name($tmparray[2],"one")})){
	    if(defined($found{&get_name($tmparray[2],"one")})){
		&add_protein(@tmparray);
	    }
	    
	    
#	else{print STDOUT " $tmparray[2] NOT FOUND \n" }
#	push @{$prots{%{$prot}->{'dbID'}}},$prot;
	}
	else{
	    my @GOs=();
	    if(/^CA/){
		my (@GOs)=(/\((.+?)\)/g);
	    }
	    elsif(/^PN\s*(.+)/){
		my $kk=$1;
		my @PROTS=split(/, /,$kk);
		my @prots;
		map{push @prots,&get_name($_),"two"}@PROTS;
		print "pp : @prots\n";
		
		my %hash;
		foreach my $goterm(@GOs){
		    map{
			$hash{"goID"}=$goterm;
			print "xx : -$_- : -$goterm-\n";
			push @{$proteins{$_}{$subonto}},\%hash;
		    }@prots;

		}		
	    }
	    else{}
	    chomp;
	}
    }
    close(A); 

}
############################################################


sub load_ancestors(){
    print STDERR "Loading ancestors...\n " if $verbose;
    if($ancestor_file){
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
	$parents{$a[0]}{$a[2]}=1;
	$children{$a[2]}{$a[0]}=1;
    }

}

############################################################
sub load_associations(){
    print STDERR "Loading associations...\n" if $verbose;

    my @a=keys(%interacting);
    map{
	my ($b,$t)=split(/xxx/);
	print " aa $_ $t : $b :  \@{$proteins{$b}{$subonto}}\n$synonyms{$t} : $synonyms{$b} : \n";
	
    }@a;

    my $table="associations_" . $species . "_go";
    if($associations_file){
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
	
	die("$fields[2] $fields[1] defined\n") if defined($associations{$fields[2]}{$fields[1]});
    }

    close(A);
}


############################################################

sub kk {

}



####################################################################### 


sub add_protein{
    my($db) = shift;
    my($dbID) = shift;
    my($old_dbSymbol) = shift;
    my $dbSymbol=&get_name($old_dbSymbol,"three");
    &debug("\nAdding -$dbSymbol- : @_\n");
    my($Qualifier) = shift;
    my $goID=shift;
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
	"db"		=> $db,
	"dbID"		=> $dbID,
	"dbSymbol"	=> $dbSymbol,
	"Qualifier"	=> $Qualifier,
	"goID"		=> $goID,
	"dbRef"		=> $dbRef,
	"evid"		=> $evid,
	"with"		=> $with,
	"aspect"	=> $aspect,
	"dbObjName"	=> $dbObjName,
	"Synonym"	=> $synonym,
	"type"   	=> $type,
	"taxon"		=> $taxon,
	"Date"		=> $Date,
	"AssBy"		=> $AssBy,
	"AnnotExt"	=> $AnnotExt,
	"GPFormId"	=> $GPFormId
        );
#    &debug("1hhaa :$dbSymbol :: $proteins{$dbSymbol}{$aspect}[0]:: $hash{'goID'}\n");
  #  die("name : -$dbSymbol-\n") unless exists($proteins{$dbSymbol}{$aspect}[0] );
    &debug("xaxa no $dbSymbol") unless exists($proteins{$dbSymbol}{$aspect}[0]);
    if ($proteins{$dbSymbol}{$aspect}[0] eq "none"){ shift(@{$proteins{$dbSymbol}{$aspect}})}
    push @{$proteins{$dbSymbol}{$aspect}},\%hash;
    &debug( "3hhaa : @{$proteins{$dbSymbol}{$aspect}}:: $hash{'goID'}\n");
}

############################################################
sub get_name{
    my $name=$_[0];
    my $old_name=$name;
    my $called_by=$_[1];
    &debug("Called by :$called_by");
    if ($synonyms{LOADED}==0){
	if(-e $synfile){
	    open(S,$synfile);
	    while(<S>){
		/^(.+)\t(.+)/ || die("Each line of the synonyms file should contain a gene name (left) and its desired synonym (right), separated by a tab.\n");
		$synonyms{$1}=$2;
	    }		
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

sub probability();
{
    load gene names from temp_asso
	foreach pair of gos in temp_asso
	      nb_asso is a list of all GOs and the number of genes annotated to each of them (implicit and explicit)



}








############################################################

sub usage{
    print STDERR "\n***** @_ *****\n\n" if @_;
    exit(1);

}


sub debug
{
    if ($debug)
    {
	print STDERR "@_\n";
    }
}

########### TODO
# 2. Use precision













