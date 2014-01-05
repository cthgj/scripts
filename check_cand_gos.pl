#!/usr/bin/perl -w
use Getopt::Std;
my %opts;

require "MY_SUBS.pl";
getopts('s:f:o:eh',\%opts);
usage() if $opts{h};
usage() unless $ARGV[0];

my $species=$opts{s}||die("Need a species -s\n");
my ($gaf_annotations_file,$synfile);
my $cands_file=$ARGV[0];
my $DATADIR=$opts{d}||"$ENV{HOME}/research/new_moon/data/";
my $genealogy_file="$DATADIR/all_three.genealogy";
my $onto=$opts{o}||undef;
my $exp=$opts{e}||undef;
my (%proteins,%ancestors,%offspring,%synonyms);
$synonyms{LOADED}=0;

## get files
if   (($species eq "human") || ($species eq "hum") || ($species eq "hs")){
    $species='human';
    $synfile="$DATADIR/human.map";
    $gaf_annotations_file=$opts{f} || "$DATADIR/gene_association.human";
}
elsif(($species eq "yeast") || ($species eq "sgd") || ($species eq "scc")){
    $species='yeast';
    $synfile="$DATADIR/yeast.map";
    $gaf_annotations_file=$opts{f} || "$DATADIR/gene_association.yeast";
}
## Collect cand names
open(A, "$cands_file")|| die("Could not open $cands_file: $!\n");
while (<A>) {
    chomp;
    my @a=split(/\t/);
    my $pair=join("_", sort  {$b lt $a} ($a[1],$a[2]));
    $proteins{$a[0]}{PAIRS}{$pair}++;
}

&load_ancestors();
&parse_gaf_file();

foreach my $cand (keys(%proteins)) {
    foreach my $pair (keys(%{$proteins{$cand}{PAIRS}})) {
	my ($go1,$go2)=split(/_/,$pair);
	my $score=0;
	$score++ if defined($proteins{$cand}{GOs}{$go1});
	$score++ if defined($proteins{$cand}{GOs}{$go2});
	print "$cand\t $pair\t$score\n";
    }
}

############################################################
sub parse_gaf_file{    
    open(A,"$gaf_annotations_file")|| die("Cannot open $gaf_annotations_file:$!\n");
    while(<A>){
	next if /^!/;
	chomp;
	my @tmparray=split(/\t/,$_);
	next if $tmparray[3] !~/^$/;
	if ($onto) {
	    next unless $tmparray[8] eq $onto;
	} 
	my $name=&get_name($tmparray[1],"parse_gaf_file");
	## If this protein is in the graph
	if(exists($proteins{$name})){
	    $proteins{$name}{GOs}{$tmparray[4]}="DIR";
	    ## Inherit ancestor annotations
	    map{$proteins{$name}{GOs}{$_}="ANC"}@{$ancestors{$tmparray[4]}};
	}

    }
    close(A);
}
############################################################
sub load_ancestors{
    open (ANC,"$genealogy_file")|| die("cannot open $genealogy_file:$!\n");
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
sub get_name{
    my $name=$_[0];
    my $old_name=$name;
    my $called_by=$_[1]||"null";
    my $reverse=$_[2]||undef;
    if ($synonyms{LOADED}==0){
	if(-e $synfile){
	    open(S,$synfile);
	    while(<S>){
		/^([^\t]+?)\t([^\t]+)\n/;
		$synonyms{NAME}{$2}=$1; 
		$synonyms{ACC}{$2}=$2; 
		$synonyms{NAME}{$1}=$1; 
		$synonyms{ACC}{$1}=$2; 
		# print "\$synonyms{NAME}{$2} : $synonyms{NAME}{$2}\n\$synonyms{ACC}{$2} : $synonyms{ACC}{$2}";
	    }
	    close(S);
	}
	else{die("Could not open $synfile\n" );	}	
 	$synonyms{LOADED}=1;
    }
    if(defined($synonyms{NAME}{$name})){
	$name=$synonyms{NAME}{$name} ;
    }
    return $name;
}
############################################################
sub usage{
   my $us="[options] <candgos file>";
   my $desc="This script will parse a .candgos file and, for each protein in that file, will print how many of the GOs assigned to that protein in the candgos file the protein was already annotated to. The format of the candgos file is (tab separated):\n\n\nPROT_A   GO1   GO2\nPROT_B   GO1   GO2\n\nIt assumes the presence of certain files in ~/research/new_moon/data";

    my %opts=(
	      's' => "Species, def:human",
	      'f' => "GAF annotations file, def:<DATADIR>/gene_association.<SPECIES>",
	      'o' => "Ontology, def:all",
	      'd' => "Data directory, def:~/research/new_moon/data/",
	      'e' => "Only count GOs from experimental data.",
	      'h' => "Print this help and exit."	      
	     );

print_help_exit(\%opts,0);

}
