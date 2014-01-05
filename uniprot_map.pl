#!/usr/bin/perl
use strict;
use warnings;
use LWP::UserAgent;
use Getopt::Std;
my %opts;
getopts('Lsvuhf:F:t:l', \%opts);
my $base = 'http://www.uniprot.org';
my $verbose=$opts{v}||undef;
my $tool = 'mapping';
my (@NAMES,@names);
my $print_only_first=$opts{l}||undef;
my $allow_trembl=$opts{u}||undef;
list_ids() if $opts{L};
&usage() if $opts{h};

&usage() unless $ARGV[0];

####################################
# Collect the IDs we want to match #
####################################
foreach my $ar(@ARGV){
    if (-e $ar){
	open(A, $ar);
	while(<A>){
	    chomp;
	    push @names,$_;
	}
    }
    else{push @names,$ar};
}

######################################
# What are we matching FROM and TO ? #
######################################
my $from=$opts{f}||'ID';
my $to = $opts{t}||'ACC';

if($opts{t}){
    $to=translate_terms($to);
}
if($opts{f}){
    $from=translate_terms($from);
}
#################
# Output format #
#################
my $format = $opts{F}||'tab';
#my @query=


my $agent = LWP::UserAgent->new;
push @{$agent->requests_redirectable}, 'POST';

my $params= {
    'from' => $from,
    'to' => $to,
    'format' => $format,
    'query' => "@names"
};

print STDERR "Requesting...\t" if $verbose;
my $try_again=1;
my $response;
while($try_again==1){
    $response = get_response($params);
    if(defined($response->is_success)){$try_again=0}
    elsif($verbose){print STDERR "Failed, retrying...\n"; }
}

if( $print_only_first) 
{
    if($response->is_success){
	my %k;
	my @a=split(/\n/,$response->content);
	foreach my $line(@a){
	    next if $line=~/From/;
	    my @b=split(/\t/,$line);
	    next if(defined($k{$b[0]}));
	    $k{$b[0]}=$b[1];
	}
	map{print STDOUT "$_\t$k{$_}\n"}keys(%k);
    }
    else{
	die 'Failed, got ' . $response->status_line . ' for ' . $response->request->uri . "\n";
    }
}
else{
    if($response->is_success){ 
	my %map;
	my @a=split(/\n/,$response->content);
	#############################
        # remove "From   To"   line #
        #############################
	shift(@a) if $a[0] =~/^From/i;
	my %kk;
	map{
	    ########################################################
            # UniProt sometimes returns the SAME line twice	   #
            ########################################################
	    $kk{$_}++;
	    #print STDOUT "$_\n" unless $kk{$_}>1; 
	    next if $kk{$_}>1; 
	    my @b=split(/\t/);
	    push @{$map{$b[0]}},$b[1];
	}@a;

	foreach my $query (keys(%map)) {
	    my %sw;
	    foreach my $target (@{$map{$query}}) {
		#############################################
                # Does this prot have SWISSPROT IDs?	    #
                #############################################
		if ($to eq 'ID')
		{
		    my @b=split(/_/,$target);
		    $sw{$target}++ if length($b[0])!=6;
		}
	    }
	    foreach my $target (@{$map{$query}}) {
		################################################
		# If a protein has a least one SWISSPROT ID    #
		# skip printing its trEMBL ones unless -u      #
		################################################
		if ($to eq 'ID' and not $allow_trembl and scalar(keys(%sw)) > 0) {
		    print "$query\t$target\n" if defined($sw{$target});
		}
		else {
		    print "$query\t$target\n";
		}
	    }
	}
	
	## Count how many were actually found
	my %found;
	foreach my $resp (@a){
	    #print STDERR "R:$resp\t";
	    foreach my $id (@names){
		#print STDERR "I:$id\n";
		$found{$id}++ if $resp=~/^$id\t/;
	    }
	    
	}
	print STDERR scalar(keys(%found)) . " of " . scalar @names . "\n";

    }
	else{
	    die 'Failed, got ' . $response->status_line . 
		' for ' . $response->request->uri . "\n";
	}
}
################################################
# This will translate $to and $from values     #
# between the names used in uniprot flat files #
# and those used by the uniprot database       #
################################################
sub translate_terms{
    my $id=shift();
    my %hash=(
	"AC" => "ACC", "ID" => "ID", "UniParc" => "UPARC", "UniRef50" => "NF50", "UniRef90" => "NF90", "UniRef100" => "NF100", "EMBL/GenBank/DDBJ" => "EMBL_ID", "EMBL/GenBank/DDBJ CDS" => "EMBLD", "PIR" => "PIR", "UniGene" => "UNIGENE_ID", "GI number*" => "P_GI", "IPI" => "P_IPI", "RefSeq" => "P_REFSEQ_AC", "PDB" => "PDB_ID", "DisProt" => "DISPROT_ID", "HSSP" => "HSSP_ID", "DIP" => "DIP_ID", "MINT" => "MINT_ID", "MEROPS" => "MEROPS_ID", "PeroxiBase" => "PEROXIBASE_ID", "PptaseDB" => "PPTASEDB_ID", "REBASE" => "REBASE_ID", "TCDB" => "TCDB_ID", "Aarhus/Ghent-2DPAGE" => "AARHUS_GHENT_2DPAGE_ID", "ECO2DBASE" => "ECO2DBASE_ID", "World-2DPAGE" => "WORLD_2DPAGE_ID", "Ensembl" => "ENSEMBL_ID", "Ensembl	Protein" => "ENSEMBL_PRO_ID", "Ensembl	Trscrpt" => "ENSEMBL_TRS_ID", "Ensembl	Genomes" => "ENSEMBLGENOME_ID", "Ensembl	Genomes" => "Protein", "Ensembl	Genomes" => "Transcript", "GeneID" => "P_ENTREZGENEID", "GenomeReviews" => "GENOMEREVIEWS_ID", "KEGG" => "KEGG_ID", "TIGR" => "TIGR_ID", "UCSC" => "UCSC_ID", "VectorBase" => "VECTORBASE_ID", "AGD" => "AGD_ID", "ArachnoServer" => "ARACHNOSERVER_ID", "CGD" => "CGD", "ConoServer" => "CONOSERVER_ID", "CYGD" => "CYGD_ID", "dictyBase" => "DICTYBASE_ID", "EchoBASE" => "ECHOBASE_ID", "EcoGene" => "ECOGENE_ID", "euHCVdb" => "EUHCVDB_ID", "EuPathDB" => "EUPATHDB_ID", "FlyBase" => "FLYBASE_ID", "GeneCards" => "GENECARDS_ID", "GeneDB_Spombe" => "GENEDB_SPOMBE_ID", "GeneFarm" => "GENEFARM_ID", "GenoList" => "GENOLIST_ID", "H-InvDB" => "H_INVDB_ID", "HGNC" => "HGNC_ID", "HPA" => "HPA_ID", "LegioList" => "LEGIOLIST_ID", "Leproma" => "LEPROMA_ID", "MaizeGDB" => "MAIZEGDB_ID", "MIM" => "MIM_ID", "MGI" => "MGI_ID", "NMPDR" => "NMPDR", "Orphanet" => "ORPHANET_ID", "PharmGKB" => "PHARMGKB_ID", "PseudoCAP" => "PSEUDOCAP_ID", "RGD" => "RGD_ID", "SGD" => "SGD_ID", "TAIR" => "TAIR_ID", "TubercuList" => "TUBERCULIST_ID", "WormBase" => "WORMBASE_ID", "WormBase Trscrpt" => "WORMBASE_TRS_ID", "WormBase Protein" => "WORMBASE_PRO_ID", "Xenbase" => "XENBASE_ID", "ZFIN" => "ZFIN_ID", "eggNOG" => "EGGNOG_ID", "HOGENOM" => "HOGENOM_ID", "HOVERGEN" => "HOVERGEN_ID", "OMA" => "OMA_ID", "OrthoDB" => "ORTHODB_ID", "ProtClustDB" => "PROTCLUSTDB_ID", "BioCyc" => "BIOCYC_ID", "Reactome" => "REACTOME_ID", "CleanEx" => "CLEANEX_ID", "GermOnline" => "GERMONLINE_ID", "DrugBank" => "DRUGBANK_ID", "NextBio" => "NEXTBIO_ID"
	);
    defined($hash{$id}) ? 
	return($hash{$id}) : 
	return($id) ;
}

sub get_response{
    my $params=shift;
    my $response = $agent->post("$base/$tool/", $params);
    while (my $wait = $response->header('Retry-After')) {
	print STDERR "Waiting ($wait)...\n";
	sleep $wait;
	$response = $agent->get($response->base);
    }
    return($response);

}





sub usage{
    $0=~s/.+?([^\/]+)$/$1/;


    print STDERR <<EndOfHelp;

$0 will connect to UniProt and map IDs from the input file to ehatever ID was requested.

USAGE:
       $0 [options] <ID_LIST>

OPTIONS:
      -h : Print this help and exit.
      -f : From ID type. This is the type of ID found in the input file.
      -F : Format. See http://www.uniprot.org/faq/28#batch_retrieval_of_entries for details. 

      -l : When a requested ID maps to many output IDs, print only the first (probably the best) match.
      -L : List all recognised ID types.
      -u : Print unreviewed (trEMBL) prots as well. By default, if an ID maps to SWISSPROT and trEMBL IDs 
           only the SWISSPROTs will be printed, use this flag to print trEMBL as well.
      -t : To ID. This is the ID type that the input IDs will be mapped to.
      -v : Verbose. Print progress messages to STDERR.
EndOfHelp

    exit(0);

}
sub list_ids{
    open(HELP, "| less") ;
    print HELP <<EndOfHelp;

Database            Abbreviation
===================== UniProt =============================
UniProtKB AC/ID	:	ACC+ID
UniProtKB AC	:	ACC
UniProtKB ID	:	ID
UniParc	        :	UPARC
UniRef50	:	NF50
UniRef90	:	NF90
UniRef100	:	NF100

=============== Other sequence databases ==================
EMBL/GenBank/DDBJ:	EMBL_ID
EMBL/GenBank/DDBJ CDS:	EMBLD
PIR		:	PIR
UniGene		:	UNIGENE_ID
GI number*	:	P_GI
IPI		:	P_IPI
RefSeq		:	P_REFSEQ_AC

================== 3D structure databases =================
PDB		:	PDB_ID
DisProt 	:	DISPROT_ID
HSSP		:	HSSP_ID

============ Protein-protein interaction databases ========
DIP		:	DIP_ID
MINT		:	MINT_ID

============ Protein family/group databases ===============
MEROPS		:	MEROPS_ID
PeroxiBase	:	PEROXIBASE_ID
PptaseDB	:	PPTASEDB_ID
REBASE		:	REBASE_ID
TCDB		:	TCDB_ID

================== 2D gel databases =======================
Aarhus/Ghent-2DPAGE:  AARHUS_GHENT_2DPAGE_ID
ECO2DBASE	:	ECO2DBASE_ID
World-2DPAGE	:	WORLD_2DPAGE_ID

============= Genome annotation databases =================
Ensembl	        :	ENSEMBL_ID
Ensembl	Protein	:	ENSEMBL_PRO_ID
Ensembl	Trscrpt	:	ENSEMBL_TRS_ID
Ensembl	Genomes	:	ENSEMBLGENOME_ID
Ensembl	Genomes	:	Protein
Ensembl	Genomes	:	Transcript
GeneID	        :	P_ENTREZGENEID
GenomeReviews	:	GENOMEREVIEWS_ID
KEGG		:	KEGG_ID	
TIGR		:	TIGR_ID	
UCSC		:	UCSC_ID	
VectorBase	:	VECTORBASE_ID

============ Organism-specific gene databases =============
AGD		:	AGD_ID	
ArachnoServer	:	ARACHNOSERVER_ID
CGD		:	CGD	
ConoServer	:	CONOSERVER_ID
CYGD		:	CYGD_ID	
dictyBase	:	DICTYBASE_ID
EchoBASE	:	ECHOBASE_ID
EcoGene		:	ECOGENE_ID
euHCVdb		:	EUHCVDB_ID
EuPathDB	:	EUPATHDB_ID
FlyBase		:	FLYBASE_ID
GeneCards	:	GENECARDS_ID
GeneDB_Spombe	:	GENEDB_SPOMBE_ID
GeneFarm	:	GENEFARM_ID
GenoList	:	GENOLIST_ID
H-InvDB		:	H_INVDB_ID
HGNC		:	HGNC_ID	
HPA		:	HPA_ID	
LegioList	:	LEGIOLIST_ID
Leproma		:	LEPROMA_ID
MaizeGDB	:	MAIZEGDB_ID
MIM		:	MIM_ID	
MGI		:	MGI_ID	
NMPDR		:	NMPDR	
Orphanet	:	ORPHANET_ID
PharmGKB	:	PHARMGKB_ID
PseudoCAP	:	PSEUDOCAP_ID
RGD		:	RGD_ID	
SGD		:	SGD_ID	
TAIR		:	TAIR_ID	
TubercuList	:	TUBERCULIST_ID
WormBase	:	WORMBASE_ID
WormBase Trscrpt:	WORMBASE_TRS_ID
WormBase Protein:	WORMBASE_PRO_ID
Xenbase		:	XENBASE_ID
ZFIN		:	ZFIN_ID

================= Phylogenomic databases ==================
eggNOG		:	EGGNOG_ID
HOGENOM		:	HOGENOM_ID
HOVERGEN	:	HOVERGEN_ID
OMA		:	OMA_ID
OrthoDB		:	ORTHODB_ID
ProtClustDB	:	PROTCLUSTDB_ID

=================== Enzyme and pathway ====================
BioCyc		:	BIOCYC_ID
Reactome	:	REACTOME_ID

================= Gene expression databases ===============
CleanEx		:	CLEANEX_ID
GermOnline	:	GERMONLINE_ID

========================== Other ==========================
DrugBank	:	DRUGBANK_ID
NextBio		:	NEXTBIO_ID

EndOfHelp
close(HELP);
    exit(0);
}
