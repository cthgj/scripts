#!/usr/bin/perl -w
use strict;
use Getopt::Std;
use feature qw(switch say);
use Data::Dumper;
use Carp;
my %opts;
require "MY_SUBS.pl";
getopts('pvdbhs:m:f:t:T:S:',\%opts);
usage() if $opts{h};
usage() unless $ARGV[0];

our $verbose=$opts{v}||undef;
our $debug=$opts{d}||undef;
my $pair_map=$opts{p}||undef;
my $missing_file=$opts{m}||undef;
my $flat_file=$opts{f}||die "Need a flat file (-f)!\n";
my $binary=$opts{b}||undef;
my $species=$opts{s}||undef;
my $to=$opts{t}||"ID";
my $tmp_dir=$opts{T}||"/tmp";
my $spacer=$opts{S}||"\t";
#my $data_dir="./uniprot";
my (%ids, %printed,%missed_prots);
my (@tbe_mapped, @pairs);

#############################################
# Create the tmp directory unless it exists #
#############################################
mkdir $tmp_dir unless -d $tmp_dir;

###########################
# Get the species' name   #
###########################
if($species){
    $species=get_species($species,1);
}
else{
    $species=guess_species($flat_file,1);
}
if ($species == -1){die("Please choose a species (-s)\n");}
#############################################
# Parse the input psicquic file and collect #
# the ids that need to be mapped. 	    #
#############################################
my $psi_file=$ARGV[0] || die "Need an input .psi file!\n";
my ($m, $s, $p)=collect_ids($psi_file);
my %tmap=%{$m};
my %seen_ids=%{$s};
my %psi_pairs=%{$p};




##########################################
# Map the missing IDs from the flat file #
# using parse_uniprot.pl                 #
##########################################
my %mapped;
##I am passing %mapped as an empty has for historical reasons.
my ($aa,$bb)=map_ids(\%tmap, \%mapped); 
%mapped=%{$aa};
my %tried=%{$bb};




##########################################
# Try and map ids directly from UniProt  #
##########################################
#($m, $s)=map_from_UniProt(\%tmap);
($m, $s)=map_from_UniProt(\%tmap, \%mapped, $to, $species, $tmp_dir, $spacer);
%mapped=%{$m};
my %not_mapped=%{$s};

##########################################################
# UniProt's mapping service does not find obsolete IDs	 #
# but UniProt's sequence retrieval service does(!). So,  #
# run uniprot_fix_ids.pl which will attempt to retrieve  #
# FASTAs of each ID and then parse the FASTA header to	 #
# get the updated ID.					 #
##########################################################
if ($to eq 'ID' || $to eq 'AC') {
    ($aa,$bb)=run_uniprot_fix(\%tried, \%mapped, "$tmp_dir/$species.names.$$.tofix", $spacer);
    %mapped=%{$aa};
    %not_mapped=%{$bb};
}

######################
# Print the map file #
######################
foreach my $id (keys(%mapped)){
    print "$id";
    #################################################
    # There are sometimes one-to-many relationships #
    #################################################
    foreach my $mapped_id (keys(%{$mapped{$id}})){
	print "\t$mapped_id";
    }
    print "\n";
}

###########################
# Print the pair map file #
###########################
my %hq=();
if($binary){%hq=%{get_hq_MIs("b")}}

if ($pair_map) {
    $psi_file=~s/.psi.gz//;
    $psi_file=~s/.psi//;
#    $psi_file=~s/[\.\/]*//;
    $psi_file=~/([^\/]*)$/;
    my $mapf=$1 . ".pairmap";
    my $fh=check_file($mapf,'w' );
    foreach my $pair (keys(%psi_pairs)) {
	my @a=split(/\t/, $pair);
	next unless $mapped{$a[0]};
	next unless $mapped{$a[1]};
	# Psi_ID1
	print $fh "$a[0]\t";
	############################################
	# Print the UniProt name(s) for each prot  #
	############################################
	local $"=",";		## print comma separated arrays 
	my @kk=keys(%{$mapped{$a[0]}});
	#Net_ID 1
	print $fh "@kk\t";
	#Psi_ID2
	print $fh "$a[1]\t";
	#Net_ID2
	@kk=keys(%{$mapped{$a[1]}});
	print $fh "@kk\t";
	#########################################
        # Print this pair's detection method(s)	#
        #########################################
	@kk=keys(%{$psi_pairs{$pair}{DET}});
	map{$_.="($hq{$_})"}@kk if $binary;
	print $fh "@kk\t";
	#####################################
	# Print this pair's supporting pubs #
	#####################################
	@kk=keys(%{$psi_pairs{$pair}{PUB}});
	#map{s/pubmed://g}@kk;
	print $fh "@kk\t";
	#################################
	# Print this pair's source DBs  #
	#################################
	@kk=keys(%{$psi_pairs{$pair}{SOURCE}});
	map{s/psi.+?\((.+?)\)/$1/g}@kk;
	print $fh "@kk\n";
	

    }
}

#########################
# Print the missing ids #
#########################
if (scalar keys(%not_mapped) >0) {
    if($missing_file){
	open(my $fh,">", "$missing_file") or die "cannot open > $missing_file: $!\n";
	v_say("\tPrinting missed...", $spacer);
	foreach my $id (keys(%not_mapped)){
	    print $fh "$id\n";
	}
	close($fh);
    }
}






## Collect the ids that we want to try and map
sub collect_ids{
    my $file=shift;
    my $fh=check_file($file,'r' );
    my (%to_map, %seen, %mapped, %pairs);
    ###############################################
    # If we are filtering by MI detection method, #
    # read the hash of accepted MIs		  #
    ###############################################
    my %hq=();
    if($binary){%hq=%{get_hq_MIs("b")}}
    ############################
    # Read/parse the .psi file #
    ############################
    v_say("Reading PSI file $file...", $spacer, 1);
    while(<$fh>){
	chomp;
#	v_say(".", 1) if $. % 10000==0;
#	v_say(". [$.]\n", 1) if $. % 100000==0;
	v_say(".",1) if $. % 100000==0;
	next if /^Total:/;
	next unless /\w/;
	####################################################
	# If this is a list of IDs taken from a psicquic   #
	# file, read the names into %seen		   #
	####################################################
	if(!/\t/){
	    $seen{$_}++;
	}
	else{
	    my ($idA, $idB, $detMethod, $orgA, $orgB, $intType, $pub, $sourceDb, $conf, $kk) = split(/\t/);
#	    print "aa $orgA : $orgB : $species\n";
	    $detMethod=~s/.+(MI:\d+).+?$/$1/;
	   
	    ################################################
            # Skip if this is not the right species	   #
            ################################################
	    $orgA=~/taxid:([-\d]+)/||do{next};
	    next unless  $1 == $species;
	    $orgB=~/taxid:([-\d]+)/||do{next};
	    next unless  $1 == $species;

	    ##############################################################
            # If we will only be dealing with binary interactions,	 #
	    # skip all lines that do not have a HQ MI id.		 #
            ##############################################################
	    $binary && do {
		next unless defined($hq{$detMethod});
	    };
	    ########################################################################
	    # Some intact interactions (maybe others) have information	           #
	    # on the isoform. Since most databases do not, ignore isoform details  #
	    ########################################################################
	    $idA=~s/-\d$//;
	    $idB=~s/-\d$//;
	    $idA=~s/\.\d$//;
	    $idB=~s/\.\d$//;
	    $seen{$idA}=$seen{$idB}=1;
	    my $pair=join("\t", sort($idA,$idB)); 
	    # Some interactions have multiple publications
	    my @ko=split(/\|/,uc($pub));
	    map{$pairs{$pair}{PUB}{$_}++}@ko;
	    $sourceDb=~s/psi.+?\((.+?)\)/$1/g;
	    ## Combinatory DBs, like MINT, often have mult. sources
	    @ko=split(/\|/,uc($sourceDb));
	    map{$pairs{$pair}{SOURCE}{$_}++}@ko;
	    $pairs{$pair}{DET}{$detMethod}++;
	}
    }
    v_say("Done", $spacer);
    close($fh);

    foreach my $id (keys(%seen)){
	my $id_extracted=0;
	#########################################
	# If this IS a uniprot accession	#
	#########################################
	if (length($id)==6 && $id=~/^[A-Z\d]+$/){
	    $to_map{AC}{$id}=$id;
	    #############################################################
	    # The %{$from{REV}} hash is a reverse map, mapping each     #
	    # extracted id (flybase, refseq, embl, whatever) back to    #
	    # the psicquic id field it was extracted from. It is used   #
	    # to avoid searching for the same prot multiple times under #
	    # different ids.     					#
	    #############################################################
	    $to_map{REV}{$id}=$id;
	    $id_extracted++;
	}
	##########################################
	# If a uniprot accession has been given	 #
	##########################################
	if ($id=~/Swiss-Prot:([\w\d]+)/i){
	    $to_map{AC}{$id}=$1;
	    $to_map{REV}{$1}=$id;
	    $id_extracted++;
	}
	if ($id=~/uniprot:([\w\d]+)/i){
	    $to_map{AC}{$id}=$1;
	    $to_map{REV}{$1}=$id;
	    $id_extracted++;
	}
	###################################
	# If a uniprot ID has been given  #
	###################################
	if ($id=~/([A-Z\d]{0,8}_[A-Z]{5})/){
	    $to_map{ID}{$id}=$1;
	    $to_map{REV}{$1}=$id;
	    $id_extracted++;
	}
	if ($id=~/flybase:([\w\d]+)/){
	    $to_map{FlyBase}{$id}=$1;
	    $to_map{REV}{$1}=$id;
	    $id_extracted++;
	}
	if ($id=~/(CG\d{4,6})/){
	    $to_map{FlyBase}{$id}=$1;
	    $to_map{REV}{$1}=$id;
	    $id_extracted++;
	}
	if ($id=~/wormbase:([\w\d]+)/){
	    $to_map{WormBase}{$id}=$1;
	    $to_map{REV}{$1}=$id;
	    $id_extracted++;
	}
	if ($id=~/wormbase:([\w\d]+\.\d+)/){
	    $to_map{WormBase}{$id}=$1;
	    $to_map{REV}{$1}=$id;
	    $id_extracted++;
	}
	if ($id=~/wormbase:\"\w+?:([\w\d]+)/){
	    $to_map{WormBase}{$id}=$1;
	    $to_map{REV}{$1}=$id;
	    $id_extracted++;
	}
	if ($id=~/locuslink:(\d+)/){
	    $to_map{GeneID}{$id}=$1;
	    $to_map{REV}{$1}=$id;
	    $id_extracted++;
	}
	if ($id=~/genbank[^:]+:(\d+)/){
	    $to_map{GeneID}{$id}=$1;
	    $to_map{REV}{$1}=$id;
	    $id_extracted++;
	}
	if ($id=~/refseq:([\w\d_]+)/){
	    $to_map{RefSeq}{$id}=$1;
	    $to_map{REV}{$1}=$id;
	    $id_extracted++;
	}
	if ($id=~/refseq:([\w\d_]+)/){
	    $to_map{RefSeq}{$id}=$1;
	    $id_extracted++;
	    $to_map{REV}{$1}=$id;
	}
	if ($id=~/emb:([\w\d]+)/i){
	    $to_map{EMBL}{$id}=$1;
	    $to_map{REV}{$1}=$id;
	    $id_extracted++;
	}
	if ($id=~/gb:([\w\d]+)/){
	    $to_map{EMBL}{$id}=$1;
	    $to_map{REV}{$1}=$id;
	    $id_extracted++;
	}
	if ($id=~/\W(ENS\w{0,5}.\d+)/){
	    $to_map{Ensembl}{$id}=$1;
	    $to_map{REV}{$1}=$id;
	    $id_extracted++;
	}
	if ($id=~/(DIP-[\w\d]+)/){
	    $to_map{DIP}{$id}=$1;
	    $to_map{REV}{$1}=$id;
	    $id_extracted++;
	}
	if ($id=~/pdb:([a-zA-Z\d]+)/i){
	    $to_map{PDB}{$id}=$1;
	    $to_map{REV}{$1}=$id;
	    $id_extracted++;
	}
	#####################################################
	# If none of the known ID types listed above	    #
	# were found, split $id into its constituent IDs    #
	# and try with each				    #
	#####################################################
	if($id_extracted==0){
	    my @ids=($id=~/:(.+?)[\|\s]/g);
	    map{
		$to_map{UNK}{$id}=$_;
		$to_map{REV}{$_}=$id;
		$id_extracted++;
	    }@ids;
	}
    }
    return(\%to_map, \%seen, \%pairs);
}



## For each type of missing ID found (EMBL, RefSeq etc), 
## print the missing id into a tmp file and then map
sub map_ids{
    my %to_map=%{$_[0]};
    my %mapped_ids=%{$_[1]};
    my %ids_seen;
    ####################################################
    # Since some psicquic id lines have >1 type of id  #
    # sort %to_map according to the number of ids of   #
    # each $from. That way, we minimise the chances of #
    # searching for the same id twice.		       #
    ####################################################
    my @froms=sort{
	if($a eq 'UNK'){return 1}
	elsif($b eq 'UNK'){return -1}
	else{
	    return(scalar(keys(%{$to_map{$b}})) <=> scalar(keys(%{$to_map{$a}})));
	}
    } keys(%to_map);
    ## For each type of ID.
from:foreach my $from (@froms){
	########################################################################
        # When the IDs to be mapped are already of the desired ID type,	       #
	# add them to the map anyway so that we can track their origin.	       #
        ########################################################################
	if ($from eq $to) {	    

	    map{
		$ids_seen{$_}++;
		$mapped_ids{$_}{$to_map{$from}{$_}}++;
	    }keys(%{$to_map{$from}});
	} 
	else {	
	    next if $from eq 'REV';
	    my @missing;
	    my %k;
	    my $counter=0;
	    map{
		$ids_seen{$_}++;
		##################################################
		# Skip if this id has already been mapped	 #
		##################################################
		unless (defined($mapped_ids{$_})) {
		    $counter++;
		    $k{$to_map{$from}{$_}}++; ## avoid duplicates
		    push @missing, $to_map{$from}{$_} unless $k{$to_map{$from}{$_}}>1;
		}
	    }keys(%{$to_map{$from}});
	    next from unless $counter>0;
	    ###################################################################
	    # Check that at least one of the missing IDs is in the flat file  #
	    ###################################################################
	    if ($#missing<=20000) { ##Avoid "Argument list too long" error
		my $list=join("\\|",@missing);
		my $ok=`grep -wm 1 '$list' $flat_file | wc -l`;
		chomp($ok);
		next from unless $ok>0;
	    }
	    ########################################
	    # Print missing ids to tmp file	   #
	    ########################################
	    my $file="$tmp_dir/$species.$$.$from.unimissing";
	    my $M_file="$tmp_dir/$species.$$.$from.unimissed";
	    open(my $fh,">", $file) or die "cannot open unimissing $file for writing : $!";
	    foreach my $missing_id (@missing) {
		unless (defined($mapped_ids{$to_map{REV}{$missing_id}})) {
		    print $fh "$missing_id\n";
		}
	    }
	    close($fh);
	    #######################################
	    # Map from the local flat file	  #
	    #######################################
	    my $vv="";
	    $vv="-v" if $verbose;
	    $vv.="-A" if $from eq 'GeneName'; ##Search aliases as well
	    $vv="-vA" if length($vv)>2;
	    my $cmd="uniprot_parse.pl -RS \"\tFlat $from \" $vv -f $from -t $to -i $file -m $M_file $flat_file";
	    debug("\t$cmd");
	    my $map=`$cmd`;
	    v_say(" ");
	    ########################################################
	    # Read the mapped names into a hash and convert        #
	    ########################################################
	    my %l;
	    map{
		my @a=split(/\t/);
		####################################################
		# UniProt sometimes returns the SAME line twice    #
		####################################################
		push @{$l{$a[0]}},$a[1];
	    }split(/\n/,$map);

	    #######################################################
	    # Now build a hash connecting the original (bad) id   #
	    # to the mapped one				          #
	    #######################################################
	    foreach my $orig_id (keys(%{$to_map{$from}})) {
		defined($l{$to_map{$from}{$orig_id}}) && do {
		    map{
			$mapped_ids{$orig_id}{$_}++;
		    }@{$l{$to_map{$from}{$orig_id}}};    
		};
	    }
	}
    }
    return(\%mapped_ids, \%ids_seen);
}

 

sub usage{
    my $us="[options] <psicquic file>";
    my $desc="This script will take a psicquic file and create a map file with the uniprot accession of each psicquic id.";
    my %opts=(
	      'usage' => $us,
	      'desc' => $desc,
	      "b" => "Map only BINARY interactions",
	      "B" => "Map only BINARY AND IP interactions.",
	      "f" => "Flat file to be used for mapping.",
	      "m" => "File to print MISSING ids into",
	      "h" => "Print this help and exit",
	      "p" => "Create map file (<psicquic file>.pairmap) connecting each psi-pair to its uniprot, including source PMID and DB. ",
	      "s" => "Spacer for STDERR indentation. Useful when $0 is called by a wrapper script.",
	      "t" => "ID type to map to (def: UniProt ID).",
	      "T" => "Directory to store tmp files in (def: /tmp)"
	     );
    print_help_exit(\%opts);

}
