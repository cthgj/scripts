#!/usr/bin/perl -w
use strict;
use Getopt::Std;
use IO::File;

my (%not_moonlighting,%seen,%moonlighting,%ancestors,%prob,%offspring,%opts,%proteins,%synonyms,%classes);
getopts('gdVvhma:D:S:p:P:G:s:f:',\%opts) || do { print "Invalid option, try 'moonlighting.pl -h' for more information\n"; exit(1); };
my (@net_probs);
my $very_verbose=$opts{V}||undef;
my $prob_file_read=0; ## is 1 if the prob file for this network has been read, gopair_probability
my $gaf_annotations_file=$opts{f} || "gene_association.goa_human";
&usage() if $opts{h};
my $DATADIR=$opts{D}||".";
my $stats_dir="$DATADIR/gostats";
my $geneology_file=$opts{g}||"$DATADIR/biological_process.geneology";

my $debug=$opts{d}||undef;
my $help=$opts{h}||undef;
my $synfile=$opts{s}||"$DATADIR/synonyms_human";
my $annotations_file=$opts{a}||undef;
my $network_file=$ARGV[0]|| die("Need a network file\n");
my $low_prob=$opts{p}||0.00005;
my $subonto=$opts{o}||"P";
my $species=$opts{S}|| "human";
my $verbose=$opts{v}||undef;
my $go_terms_file=$opts{G}||"$DATADIR/GO.terms_ids_obs";
my $go_percentages=undef;
my $pairs=1;
if($opts{m}){$pairs=undef; $go_percentages=1};
$synonyms{LOADED}=0;
if($opts{m}){$pairs=undef; $go_percentages=1};
if($subonto eq "p"){$subonto="P"}
elsif($subonto eq "c"){$subonto="c"}
elsif($subonto eq "f"){$subonto="f"}
elsif($subonto eq "F" || "P" || "C"){}
else{print STDERR "Subontology (-o) must be one of \"p,c,f\"\n"; exit(1)}
$verbose=1 if $debug;
my $date=localtime;
my $LOG = IO::File->new(">> $DATADIR/moonlighting.log")|| die("Cannot open moonlighting.log : $!\n");
print $LOG "********** $date ******************
Network file : $ARGV[0]\n"; 
map{print $LOG "$_ : $opts{$_}\n"}keys(%opts);
print $LOG "\n\n"; 


my $multi_only=1;
my $use_classes=0;
############################### Main program ###############################
&load_network($network_file);
&parse_gaf_file();
&load_ancestors();
&load_annotations($annotations_file) if $annotations_file;


my $protnum=scalar(keys(%proteins));
my $c=0; ## counter
## look for interactions between different classes
if ($use_classes){
    bait:foreach my $bait (keys(%proteins)){
	$c++;
	target:foreach my $target (keys(%{$proteins{$bait}{INTERACTORS}})){
	    next target if defined($not_moonlighting{$bait}{$target});
	    ##skip target if it shares a class with bait;
	    foreach my $class (keys(%{$proteins{$bait}{CLASS}})){
		next target if defined($proteins{$target}{CLASS}{$class});
	    }
	    foreach my $bgo (keys(%{$proteins{$bait}{GOs}})){
		foreach my $tgo (keys(%{$proteins{$target}{GOs}})){
		    printf STDERR ("$c of $protnum\r");
		    my $P=&gopair_probability($bgo,$tgo);
		    if ($P<=$low_prob){
			## Check ancestor GOs
			&check_ancestor_gos($bgo,$tgo);
			
			$moonlighting{$bait}{$target}="$bait ($bgo)  \t:\t$target ($tgo)  \t: $P";
			print STDERR "$bait ($bgo)  \t:\t$target ($tgo)  \t: $P\n"
		    }
		    else{
			$not_moonlighting{$bait}{$target}=1;
			$moonlighting{$bait}{$target}=undef;
			next target;
		    }
		}
	    }
	} ## end foreach target
    } ## end foreach bait
} ## end if $use_classes
elsif($multi_only){
  bait:foreach my $bait (keys(%proteins)){
      $c++;
      next bait unless keys(%{$proteins{$bait}{CLASS}}) > 1;
    target:foreach my $target (keys(%{$proteins{$bait}{INTERACTORS}})){
	next target if defined($moonlighting{$bait}{$target});
	next target if defined($moonlighting{$target}{$bait});
	&debug("$bait ::: $target");
	printf STDERR ("$c of $protnum\r");
	next target if defined($not_moonlighting{$bait}{$target});
	foreach my $bgo (keys(%{$proteins{$bait}{GOs}})){
	    next if $bgo eq 'GO:0008150';
	    foreach my $tgo (keys(%{$proteins{$target}{GOs}})){
		next if $tgo eq 'GO:0008150';
		my $P=&gopair_probability($bgo,$tgo,"main");
		if ($P<=$low_prob){
		    if (&check_ancestor_gos($bgo,$tgo,$bait,$target) == 1){
			$moonlighting{$bait}{$target}="$bait ($bgo)  \t:\t$target ($tgo)  \t: $P";
			#print STDERR "got $bait $target $bgo $tgo\n" if $verbose;
			next bait;
		    }
		}
		else{
		    $not_moonlighting{$bait}{$target}=1;
		    $moonlighting{$bait}{$target}=0;
		    next target;
		}
	    }
	}
    } ## end foreach target
  } ## end foreach bait
} ## end multi_only 

print "\n";
foreach my $bait (keys(%moonlighting)){
#    print "bb : $bait\n";
    foreach my $target (keys(%{$moonlighting{$bait}})){
	next if defined($not_moonlighting{$bait}{$target});
	print "$moonlighting{$bait}{$target}\n";

    }
}

## Create network specific probability file
#unless (-e "$network_file.prob"){
    open(NET,">$network_file.prob");
    map{print NET}@net_probs;
    system("cat $network_file.prob | sort|uniq > lililolo; mv lililolo $network_file.prob");
    close(NET);
#}


############################### SUBROUTINES ################################
sub load_network {
    print STDERR "Loading Network...\n" if $verbose;
    my $A = IO::File->new("< $network_file")|| die("Cannot open $network_file : $!\n");
    while(<$A>){
	next if /^\d+$/;
	my ($bait,$target)=split(/\s+/,$_);
	$bait=&get_name($bait,"load_network");
	$target=&get_name($target,"load_network");
	## %proteins will hold all the info on all the 
	## network's proteins
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
	else{die("Could not open $synfile\n" );
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
    }
    elsif(defined($synonyms{$name})){
	$name=$synonyms{$name} ;
    }
    elsif ($name=~/(.+)_$species/i && defined($synonyms{$1})){
	$name=$synonyms{$1};
    }
    else{$synonyms{$name}=$name;}
    if ($very_verbose){&debug("Called by : $called_by for $old_name ($synfile), returned $name"); }
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
	next if $tmparray[3] =~/NOT/;
	my $name=&get_name($tmparray[1],"parse_gaf_file");
	## If this protein is in the graph
	if(exists($proteins{$name})){
	    $proteins{$name}{GOs}{$tmparray[4]}++;
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
	my $child=shift(@terms);
	$ancestors{$child}=[@terms];
	map{$offspring{$_}{$child}++}@terms;
    }
    close(ANC);
}
############################################################

sub load_annotations{
    print STDERR "Loading annotations...\n" if $verbose;
    my $prot;
    open(A,"$annotations_file")|| die("Cannot open $_[0]:$!\n");
    my $k=0;
    while(<A>){
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
	map{$proteins{$_}{CLASS}{$k}++}@b;
	$k++;
    }
}
############################################################
sub gopair_probability{

    my $go1=shift;
    my $go2=shift;
    my $source=shift;
    my $low;
    my $gopair=$go1 . "xxx" . $go2;
    my $gopair2=$go2 . "xxx" . $go1;
#    print "$gopair\n";
    if ($go1 eq $go2){
	$prob{$gopair2}=$prob{$gopair}=1;
    }
    unless(defined($prob{$gopair})){
	die("crap ::: $gopair ::: $prob{$gopair} ::: $prob{$gopair2}\n") if defined($prob{$gopair2});
	if ((-e "$network_file.prob") && ($prob_file_read==0)){
	    open(A,"$network_file.prob");
	    print STDERR "\nReading $network_file.prob ($gopair)..." if $verbose;
	    while(<A>){
		chomp;
		my @tt=split(/\t/);
		my $tempgo1=$tt[0] . "xxx" . $tt[1];
		my $tempgo2=$tt[1] . "xxx" . $tt[0];
		$prob{$tempgo2}=$prob{$tempgo1}=$tt[2];
	    }
	    close(A);
	    $prob_file_read=1;
	    print STDERR "\n Done\n" if $verbose;
	}
	$go1=~/GO:((.....).+)$/||die("cannot match go1 : $go1\n");
	my ($a1,$a2)=($1,$2);
	$go2=~/GO:((.....).+)$/||die("cannot match go2 : $go2\n");
	my ($b1,$b2)=($1,$2);
	## if this prob was not in net prob file
	unless(defined($prob{$gopair})){
	    ## I have calculated the probabilities of ALL possible
	    ## GO pairs and split them into files according to the GO
	    ## number. So, GO:0019952 will be in the file 00199
	    my $file;
	    $a1>$b1 ? ($file=$b2 . ".prob") : ($file=$a2 . ".prob")  ;
	    if(-e "$stats_dir/$file.gz") {
		open(A,"zcat $stats_dir/$file |")|| die("cannot open $stats_dir/$file : $!\n$a1,$a2:$b1,$b2,$go1,$go2,$gopair\n");
	    }
	    elsif(-e "$stats_dir/$file") {
		open(A,"$stats_dir/$file")|| die("cannot open $stats_dir/$file : $!\n$a1,$a2:$b1,$b2,$go1,$go2,$gopair\n");
	    }
	    else{die("bad stats filename: $stats_dir/$file:$1\n")}
	    while(<A>){
		chomp;
		if(( /$b1/) &&( /$a1/)){
		    my @tt=split(/\t/);
		    $prob{$gopair2}=$prob{$gopair}=$tt[2];
		    last;
		}
	    }
	    close(A);
	} ## end unless(defined($prob{$gopair})) 2	
    } ## end unless(defined($prob{$gopair})) 1
    
    unless( (defined($prob{$gopair})) || (defined($prob{$gopair2})) ) {
	die("crap2 : $gopair/$gopair2: " . $prob{$gopair} . ":". $prob{$gopair2} ."\n") ;
	print STDERR "Probability problem, missing $gopair\n";
	$prob{$gopair}=$prob{$gopair2}=-1;
    }
    ## @net_probs will be used to create a prob file for this network
    push @net_probs, "$go1\t$go2\t" . $prob{$gopair} . "\n";
    return($prob{$gopair});
}
############################################################

sub check_ancestor_gos{
    my $go1=shift;
    my $go2=shift;
    my $bait=shift;
    my $target=shift;
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
	   #   print STDERR "$n/$k $#go1s/$#go2s $bait :: $target\n";
#	      printf STDERR ("$c of $protnum $k/$#go2s\r");
	      if(&gopair_probability($go1s[$n],$go2s[$k],"anc")>$low_prob){
		  $seen{$go1}{$go2}=$seen{$go2}{$go1}=0;
		  print $LOG "Skipped $bait,$target ($go1s[$n],$go2s[$k])\n";
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
