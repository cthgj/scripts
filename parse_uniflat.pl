#!/usr/bin/perl
use strict;
use warnings;
use Data::Dumper;
use Getopt::Std;
use SWISS::Entry;
my %opts;
getopts('vhIFef:n:t:i:',\%opts) || do { print "Invalid option\n"; exit(1); };
&usage() if $opts{h};
&print_ID_types if $opts{I};
my $toID=$opts{t}||"AC";
my $fromID=$opts{f}||'ID';
my $flat=$ARGV[0]||die("Need a flat file\n");
my $names = $opts{n}||undef;
my $IDlist;
my $verbose=$opts{v}||undef;
my %entries;
my %ids;
my $return_fasta=$opts{F}||undef;
my $looking_for;
my $return_entry=$opts{e}|| undef;

## get desired ids
if($names){
    map{
	$ids{$_}++;
	$looking_for++;
    }split(/\s+/,$names);
}
else{
    $IDlist=$opts{i}||die("Need a list of desired IDs. Either as a text file (-i) or as a space separated list (-n)\n");
## Read desired IDs. Expecting a text file, one ID per line
    open(ID,$IDlist)||die("Cannot open ID list $IDlist :$!\n");
    while(<ID>){
	chomp;
	$ids{$_}++;
	$looking_for++;
    }
}


# Change the line termination string so we read an entire entry at a time
local $/ = "\n//\n";

if($flat=~/\.gz$/){
    open(A,"zcat $flat |")|| die("cannot open flat file $flat : $!\n");
}
else{
    open(A,"$flat")|| die("cannot open flat file $flat : $!\n");
}

while (<A>) {
    exit(0) if $looking_for==0;
    printf STDERR "($looking_for,$.)\r" if $verbose;
    # Read in all the entries and fill %entries 
    my $entry = SWISS::Entry->fromText($_);
    my $id;
    my $found=0;
    foreach my $key ($entry->IDs->elements){$id=$key}
    &get_all_ids($entry,$id) unless $fromID eq 'ID';
    if(($return_entry) && (defined($ids{$id}))){
	print STDOUT $entry->toText;
	$looking_for--;
	$ids{FOUND}{$id}++;
	$found=1;
    }
    elsif($return_fasta){
	if(defined($ids{$id})){
	    my $l = length($entry->SQs->seq);
	    my $i = 0;
	    print ">$ids{$id} $id \n";
	    while($i<$l)   ## convert to FASTA
	    {
		print substr($entry->SQs->seq,$i,60) . "\n"; 
		$i=$i+60;
	    }
	    $looking_for--;
	    $ids{FOUND}{$id}++;
	    $found=1;
	}
    }
    else{
	if(defined($ids{$id})){
	    $looking_for--;
	    $ids{FOUND}{$id}++;
	    $found=1;
	    if($toID eq "AC"){
		my $array=$entry->ACs->list;
		printf STDERR "($looking_for)\r" if $verbose;
		print "$id\t${$array}[0]\n";
		$ids{FOUND}{${$array}[0]}++;
	    }
	    else{
		
		foreach my $list ($entry->DRs->list){
		    foreach my $array (@{$list}){
			foreach my $name (@{$array}){
			    print "nn : $name\n";
			}
			if(${$array}[0] eq $toID){
			    print STDOUT "$id\t$ids{$id}\t";
			    if($toID eq "WormBase"){
				map{print STDOUT "$id\t$_\n" if /^WBG/}@{$array};
			    }
			    elsif($toID eq "HGNC"){
				for(my $n=1; $n<$#{$array}; $n++){
				    print STDOUT "$$array[$n]\t";
				}
				print STDOUT @{$array}[$#{$array}], "\n";
			    }
			    else{
				print STDOUT "${$array}[1]\n";
			    }
			}
		    }
		}
	    }
	}
	## If we have been given a gene name
	elsif((defined($entry->GNs->text)) &&  defined( $ids{$entry->GNs->text})){
	    print STDOUT $entry->GNs->text,"\t$id\n"; 
	    printf STDERR "($looking_for)\r" if $verbose;
	    $ids{FOUND}{$entry->GNs->text}++;
	    $looking_for--;
	    $found=1;
	}
	else{
	    ## foreach DR line
	    foreach my $array ($entry->DRs->elements){
		## foreach DR element in DR line
		if(${$array}[0] eq $toID){
		    for (my $n=0;$n<scalar(@{$array});$n++){
			if(defined($ids{$$array[$n]})){
			    $looking_for--;
			    $ids{FOUND}{$id}++;
			    print STDOUT "${$array}[1]\t${$array}[2]\t$id\t", ${$entry->ACs->list}[0],"\n";
			    printf STDERR "($looking_for)\r" if $verbose;
			    $found=1;
			}
		    }
		}
	    }
	    # print "xx " , $entry->GNs->text , "\n";
	    # foreach my $key ($entry->GNs->elements){
	    # 	my $a=$key->get('Names');
	    # 	print "aa : $a\n";
	    # 	print "kk : " , $key->list ,"\n";
	    # 	print "ll : ",  @{$key->Names} , "\n";
	    # 	map{print "oo $_ : $$key{$_}\n"}keys(%{$key});
	    # }
	    # my $aa=$entry->GNs->elements ;
	    # print STDOUT "aa @{$aa} bb \n"; die();
	}
	

	## If we haven't found anything interesting in this
	## entry, try gene synonyms
	unless($found==1){
	    if(defined($entry->GNs->text)){
		my @names=split(/\s*OR\s*/,$entry->GNs->text);
		map{
		    if(defined($ids{$_})){
			print "$_\t$id\n";
			$looking_for--;
			$found=1;
		    }
		}@names;
	    }
	}
    }
}

unless($looking_for==0){
    print STDERR "$looking_for IDs were not found:\n";
    map{
	unless(/FOUND/){
	    print STDERR "$_ " unless defined($ids{FOUND}{$_})
	}
    }keys(%ids),
    
}
print STDERR "\n";
exit();
sub get_all_ids{
    my $entry=shift;
    my $ID=shift;
    ## foreach DR line
    foreach my $array ($entry->DRs->elements){
	## foreach DR element in DR line
	my $toID=${$array}[0];
	if(defined($ids{$$array[1]})){
	    $ids{$ID}=$$array[1];
	}
    }
    foreach my $array ($entry->ACs->list){
	if(defined($ids{$$array[0]})){
	    $ids{$ID}=${$array}[0];
	}

    }
}
sub usage{ 
    $0=~s/.+?([^\/]+)$/$1/;
    open(HELP, "| more") ;
    print HELP <<EndOfHelp;

$0 parses UniProt Flat Files an returns either entire entries, or an entry\'s 
FASTA sequence or translates between IDs.

USAGE:  
        $0 [options] <FLAT FILE>

EXAMPLES:
  
    $0 -n CATA_HUMAN  -t AC human.flat
    $0 -n P04040  -t AC human.flat
    $0 -n HS_CAT -t AC human.flat
    $0 -fn HS_CAT human.flat


OPTIONS:
 	-e : Return entire entry in Flat File format.
	-f : What ID type the input IDs are.
	-F : Return FASTA sequence.
	-h : Print this help and exit.
	-i : List of IDs to look for.
	-I : Print the (long) list of IDs $0 can cope with.
	-n : What follows is a QUOTED, space separated list of the names desired.
	-t : What ID type the input IDs should be converted to.
	

EndOfHelp
close(HELP);
    exit();
}

sub print_ID_types{
    $0=~s/.+?([^\/]+)$/$1/;
    print  <<EndOfHelp;
These are the various database identifiers that $0 can deal with. 
By default, if no other options are given, it will tra translate 
between UniProt accessions and ids.


EXAMPLES:
  
    $0 -n CATA_HUMAN  -t AC human.flat
    $0 -n P04040  -t AC human.flat
    $0 -n HS_CAT -t AC human.flat


    Database Name              Example Entry
Aarhus/Ghent-2DPAGE	:	1524
ArrayExpress		:	P04040
Bgee			:	P04040
BioCyc			:	MetaCyc:MONOMER66-341
BRENDA			:	1.11.1.6
CleanEx			:	HS_CAT
CTD			:	847
DrugBank		:	DB01213
eggNOG			:	prNOG05863
EMBL			:	X04085
Ensembl			:	ENST00000241052
Gene3D			:	G3DSA:2.40.180.10
GeneCards		:	GC11P034460
GeneID			:	847
Genevestigator		:	P04040
GermOnline		:	ENSG00000121691
GO			:	GO:0005782
HGNC			:	HGNC:1516
H-InvDB			:	HIX0009550
HOVERGEN		:	HBG003986
HPA			:	CAB001515
InParanoid		:	P04040
IntAct			:	P04040
InterPro		:	IPR002226
IPI			:	IPI00465436
KEGG			:	hsa:847
MIM			:	115500
MINT			:	MINT-1210583
NextBio			:	3550
neXtProt		:	NX_P04040
NMPDR			:	fig|9606.3.peg.5453
OGP			:	P04040
OMA			:	SHLAARE
Orphanet		:	926
OrthoDB			:	EOG9KPWX2
PANTHER			:	PTHR11465
Pathway_Interaction_DB	:	foxopathway
PDB			:	1DGB
PDBsum			:	1DGB
PeptideAtlas		:	P04040
PeroxiBase		:	5282
Pfam			:	PF00199
PharmGKB		:	PA26099
PhosphoSite		:	P04040
PhylomeDB		:	P04040
PIR			:	A23646
PRIDE			:	P04040
PRINTS			:	PR00067
PROSITE			:	PS00437
ProteinModelPortal	:	P04040
Reactome		:	REACT_1698
RefSeq			:	NP_001743.1
REPRODUCTION-2DPAGE	:	IPI00465436
SMR			:	P04040
STRING			:	P04040
SUPFAM			:	SSF56634
SWISS-2DPAGE		:	P04040
UCSC			:	uc001mvm.1
UniGene			:	Hs.502302

EndOfHelp
exit();
}

