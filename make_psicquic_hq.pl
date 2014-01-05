#!/usr/bin/perl -w
use strict;
use Getopt::Std;
use feature qw(switch say);
use Data::Dumper;


my %opts;
getopts('BhvbdlDT:t:s:p:mP:M:S:L:f:',\%opts);
our $verbose=$opts{v}||undef;
our $debug=$opts{d}||undef;
require "MY_SUBS.pl";
my $binary=$opts{b}||0;
my $map_ids=$opts{m}|| undef;
my $map_file=$opts{M}||undef;
my $min_pubs=$opts{p}||0;
my $min_hq_pubs=$opts{P}||0;
my $tmp_dir=$opts{T}||"/tmp";
my $data_dir=$opts{D}||"./uniprot";
my $species=$opts{s}||undef;
my $spacer=$opts{S}||" ";
my $keep_log=$opts{l}||undef;
my $to_id=$opts{t}||'ID';
##################################################
# If no map file is given (-M), then IDs must be #
# mapped and $map_ids is set to 1		 #
##################################################
unless ($opts{M}) {
    $map_ids=1;
}

usage() if $opts{h};
##############################################################
# Get the list of high quality interaction detection methods #
##############################################################
my %hq;
if ($opts{B}) {
     %hq = %{get_hq_MIs('e')};
     $binary=1;
}
else {
    %hq = %{get_hq_MIs('b')};
}
if ($opts{D}) {
    if ($opts{B} || $opts{b}) {
	print STDERR "## The following methods are considered high quality:\n";
        map{print "$_ : $hq{$_}\n";}keys(%hq);

    }
    else {
	print STDERR "ALL interaction detection methods will be used unless a binary option (-b or -B) is set.\n";
    }
	exit;
}
$keep_log=1 if $opts{L};
my $log_file=$opts{L}||"$ARGV[0]" . ".hq.log";
usage() unless $ARGV[0];


###########################
# Get the species' name   #
###########################
my $taxid;
if($species){
    $species=get_species($species,"-s");
    $taxid=get_species($species,"-s",1);
}
else{
    $species=guess_species($ARGV[0]);
    $taxid=get_species($species,"-s",1);
}
if ($species eq '-1'){die("Please specify a species (-s)\n");}

my $uni_suffix=get_suffix($species)||undef;

my $flat_file=$opts{f}||"$data_dir/$species.flat";

my (%skipped,%pairs, %seen, %pubs, %printed, %map);
my (@tbe_mapped, @pairs);
##################################
# Read the map file, if provided #
##################################
if($map_file){
    open(MAP, $map_file) or die "Could not open map file $map_file: $!\n";
    while(<MAP>){
	chomp;
	my @a=split(/\t/);
	next if $a[1] =~/^\s*deleted\s*$/i;
	next if $a[1] =~/^\s*Fragment\s*$/i;
	$map{Psi2Uni}{$a[0]}=$a[$#a];
	$map{Uni2Psi}{$a[$#a]}=$a[0];
	###############################################
        # If a map file has been given and 	      #
        # we are not mapping to uniprot, unset	      #
	# uni_suffix (used later on)		      #
        ###############################################
	if (defined($uni_suffix) && ($a[$#a]!~/$uni_suffix/)) {
	    $uni_suffix=undef;
	}
    }
} 
############################################
# Create a map file if none has been given #
############################################
elsif($map_ids){%map=%{id_mapper(@tbe_mapped)}}


############################
# Parse the psicquic file  #
############################
my $fh=check_file($ARGV[0],'r' );
while(<$fh>){
    chomp;
    next if /^Total:/;
    next unless /\w/;
    ## Apid has an error where Q9PH20 is reported as Q9PH0, fix
    s/\bQ9PH0\b/Q9PH20/g;

    v_say("Psi file, line : [$.]\r", $spacer, 1) if $. % 100000==0;
    my ($idA, $idB, $detMethod, $orgA, $orgB, $intType, $pub, $sourceDb, $conf) = split(/\t/);
    
    #########################################
    # Skip species we are not interested in #
    #########################################
    $orgA=~s/taxid:(\d+).*/$1/;

    next unless $orgA eq $taxid;
    ############################################################
    # Some intact interactions (maybe others) have information #
    # on the isoform. Remove because very few have this.       #
    ############################################################
    $idA=~s/-\d$//;
    $idB=~s/-\d$//;
    $idA=~s/\.\d$//;
    $idB=~s/\.\d$//;

    ## This should already have been done by 
    ## get_psicquic_interactome.pl. But, just in case...
    $detMethod=~s/.+(MI:\d+).+?$/$1/;
    
    ############################################################
    # Skip interactions for which no detection method is given #
    ############################################################
    next unless $detMethod;
    ###############################################
    # If we have already built the mapfile, skip  #
    # any ids that have not been correctly mapped #
    ###############################################
    my $ll=defined($map{$idA});
    my $kk=defined($map{$idB});
    if($map_file){
    	unless (defined($map{Psi2Uni}{$idA}) && defined($map{Psi2Uni}{$idB})){next}
    }

    #############################################
    # $pair is the pair of interacting proteins #
    #############################################
    my $pair=join("\t", sort {$b lt $a} $idA,$idB);
   
    # ##########################################################
    # # If we want to map IDs, and these are not uniprot ACCs  #
    # ##########################################################
    # if($map_ids && (length($idA)!=6 || length($idB)!=6 )){
    # 	my $changed=0;
    # 	given($idA){
    # 	    when (/Swiss-Prot:([\w\d]+?)/){
    # 		my $s=$1;
    # 		$pair=~s/$idA/$s/;
    # 		$changed++;
    # 	    }
    # 	    when (/uniprot:([\w\d]+?)/){
    # 		my $s=$1;
    # 		$pair=~s/$idA/$s/;
    # 		$changed++;
    # 	    }
    # 	}
    # 	given($idB){
    # 	    when (/uniprot:([\w\d]+?)/){
    # 		my $s=$1;
    # 		$pair=~s/$idB/$s/;
    # 		$changed++;
    # 	    }
    # 	    when (/Swiss-Prot:([\w\d]+?)/){
    # 		my $s=$1;
    # 		$pair=~s/$idB/$s/;
    # 		$changed++;
    # 	    }
    # 	}
    ##############################################
    # Always map ids, to get the most up to date #
    # accession.				 #
    ##############################################
    push @tbe_mapped, $pair;
    # defined($pub) && do {push @{$pairs{$pair}{PUBs}{$pub}}, $detMethod};
    # $pairs{$pair}{MIs}{$detMethod}++;

    ###############################################################
    # Collect the publication this interaction was reported in	  #
    ###############################################################
    defined($pub) && do {
	my @ppubs;
	while($pub=~/pubmed:(\d+)/g){push @ppubs,$1}; 
	my $c=0;
	map{
	    next if $_ eq '-';
#	    push @{$pairs{$pair}{PUBs}{$pub}}, $detMethod;
	    $pairs{$pair}{PUBs}{$_}{$detMethod}++;
	    $c++;
	}@ppubs
    };
    ########################################################
    # Collect the detection method for this interaction	   #
    ########################################################
    if ($binary){
	$pairs{$pair}{MIs}{$detMethod}++ if defined($hq{$detMethod});
    }
    else {
	$pairs{$pair}{MIs}{$detMethod}++;	
    }
    ###############################################
    # Collect the database this pair was found in #
    ###############################################
    $sourceDb=~s/psimi/psi-mi/g;
    $sourceDb=lc($sourceDb);
    if ($sourceDb=~/\|/) {
	my @a=split(/\|/, $sourceDb);
	map{$pairs{$pair}{SOURCE}{$_}++}@a;
    }
    else {
	$pairs{$pair}{SOURCE}{$sourceDb}++;	
    }
}



#####################################################
# Do we have a minimum publication threshold?       #
#####################################################
my $pp=0;
if ($min_pubs) {
    $pp=$min_pubs;
}
elsif ($min_hq_pubs) {
    $pp=$min_hq_pubs;
}
#########################################################
# To be considered HQ, an interaction must have a score #
# greater than this threshold.			        #
#########################################################
my $threshold=0+$binary+$pp;

########################################################
# Main loop, go through each of the interacting pairs  #
# and keep those whose score passes the threshold      #
########################################################
my %seen1;

####################
# Print log header #
####################
my $log;
if ($keep_log) {
    open($log, ">", $log_file)|| die "Could not open logfile $log_file for writing : $!\n";
    print $log "# prot1     \tprot2     \tSupporting_Publications\tDetection_Methods\tSource_DBs\tprot1_IDs_from_psicquic\t   prot2_IDs_from_psicquic\n";
}

pair:foreach my $pair (keys(%pairs)){  
    my $has_bin=0;
    my $hq_pubs=0;
    my $pubs=0;
    #print "PAIRS ($pair): " . Dumper(\%{$pairs{$pair}}) . "\n";
   #######################################
   # If we only want BINARY interactions #
   #######################################
   if($binary){
       mi:foreach my $mi (keys(%{$pairs{$pair}{MIs}})){
	   if (defined($hq{$mi})){
	       $has_bin=1;
	       last mi; 
	   }
       }
   }
   #########################################
   # If we want a minimum of publications  #
   # that identify BINARY interactions	   #
   #########################################
   if($min_hq_pubs){
       foreach my $pub (keys(%{$pairs{$pair}{PUBs}})){
	   my $cc=0;
	   foreach my $mi (keys(%{$pairs{$pair}{PUBs}{$pub}})){
	       $cc++ if defined($hq{$mi});
	   }
	   $hq_pubs++ if $cc>0;
       }
   }
   ############################################################
   # If we just want a minimum of publications, binary or not #
   ############################################################
   if($min_pubs){
       my $kkk=scalar(keys(%{$pairs{$pair}{PUBs}}));
       my @a=(keys(%{$pairs{$pair}{PUBs}}));
       $pubs+=scalar(keys(%{$pairs{$pair}{PUBs}}));
   }
   ####################################
   # Only keep pairs with uniprot IDs #
   ####################################
   my ($idA,$idB)=split(/\t/,$pair);
   my $uni_pair;
   if(defined($map{Psi2Uni}{$idA}) && defined($map{Psi2Uni}{$idB})){
       if (defined($uni_suffix)){
	   ######################################################
           # APID, and maybe others has some sequences 	        #
	   # reported as being of the wrong species. e.g        #
	   # P29288 (PPA5_RAT) is reported as human.	        #
           ######################################################
	   next pair unless $map{Psi2Uni}{$idA}=~/$uni_suffix/;
	   next pair unless $map{Psi2Uni}{$idB}=~/$uni_suffix/;
       }
       $uni_pair=join("\t", sort {$b lt $a} $map{Psi2Uni}{$idA},$map{Psi2Uni}{$idB});
   }
   else{ next pair; }
   ##############################################
   # If this pair's score passes the threshold  #
   # and we have not seen it before, keep.      #
   ##############################################
    my $score=0;
    $binary && do {$score+=$has_bin};
    $min_hq_pubs && do {
	if ($hq_pubs>$pp) {
	    $hq_pubs=$pp;
	}
	$score+=$hq_pubs;
    };
    $min_pubs && do {
	if ($pubs>$pp) {
	    $pubs=$pp;
	}
	$score+=$pubs;
    };
    if($score>=$threshold && not $seen{$uni_pair}){
       print "$uni_pair\n";
       $seen{$uni_pair}++;
       ################################################################
       # For each pair of interacting proteins, collect their	      #
       # original psicquic ID, supporting publications, and source DB #
       ################################################################
       if ($keep_log) {
	   print $log "$uni_pair\t" . join(",",keys(%{$pairs{$pair}{PUBs}})) . "\t";
	   
	   # Just bullshitting to get the format right
	   my @a=keys(%{$pairs{$pair}{MIs}});
	   my $kk=shift(@a);
	   defined($hq{$kk}) ? 
	       print $log "$kk($hq{$kk})" :
		   print $log "$kk" ;
	   map{
	       defined($hq{$_}) ? 
		   print $log ",$_($hq{$_})" :
		       print $log ",$_";
	   }@a;
	   print $log  join(",",keys(%{$pairs{$pair}{MIs}})) . "\t" . 
	       join(",",keys(%{$pairs{$pair}{SOURCE}}))  . "\t";
	   my ($id1,$id2)=split(/\t/,$uni_pair);
	   print $log join(",",@{$map{Uni2Psi}{$id1}}) . "\t" . join(",",@{$map{Uni2Psi}{$id2}}) . "\n";
       }
   }


   $seen1{$uni_pair}++;
   $seen1{$pair}++;

   # ## if these are not uniprot ACCs
   # if(length($idA)!=6 || length($idB)!=6 ){
   #     print "$pair : $idA $idB\n";
   #     defined($map{$idA}) && do {
   # 	   $idA = $map{$idA};
   #     };
   #     defined($map{$idB}) && do {
   # 	   $idB = $map{$idB};
   #     };
   #     if(length($idA)==6 && length($idB)==6 ){
   # 	   $pair=join("\t", sort {$b lt $a} $idA,$idB);
   # 	   if($score>=$threshold && not $seen{$pair}){
   # 	       print "$pair\n";
   # 	       $seen{$pair}++;
   # 	   }
   #     }
   # }
   # else{
   #     if($score>=$threshold && not $seen{$pair}){
   # 	   print "$pair\n";
   # 	   $seen{$pair}++;
   #     }
   # }
}
v_say(" ");
## map ids to uniprot
sub id_mapper{
    #########################################
    # Die if the flat file cannot be found. #
    #########################################
    my $nn=check_file($flat_file, "r");
    do {
	print STDERR "Flat file $flat_file could not be opened, please specify a different file (-f): $!\n" ;
	exit(0);
    } unless $nn;
    
    my (%tmap, %missed_prots, %prot_map);
    my @ids=map{split(/\t/)}@_;
    foreach my $id (@ids){
	my $c=0;
	given ($id){
	    when (/locuslink:(\d+)/){
		$tmap{GeneID}{$id}=$1;
		$c++;
	    }
	    when (/genbank[^:]+:(\d+)/){
		$tmap{GeneID}{$id}=$1;
		$c++;
	    }
	    when (/refseq:([\w\d_]+)/){
		$tmap{RefSeq}{$id}=$1;
		$c++;
	    }
	    when (/emb:([\w\d]+)/i){
		$tmap{EMBL}{$id}=$1;
		$c++;
	    }
	    when (/gb:([\w\d]+)/){
		$tmap{EMBL}{$id}=$1;
		$c++;
	    }
	    when (/\W(ENS.\d+)/){
		$tmap{Ensembl}{$id}=$1;
		$c++;
	    }
	    when (/pdb:([\w\d]+)/i){
		$tmap{PDB}{$id}=$1;
		$c++;
	    }
	    ## Get UniProt ACCs
	    when (length==6){
		$tmap{AC}{$id}=$id;
	    }
	}
    }
    ##########################################################
    # For each type of missing ID found (EMBL, RefSeq etc),  #
    # print the missing id into a tmp file and then map      #
    ##########################################################
    foreach my $from (keys(%tmap)){        
	my @missing;
	my %k;
	map{
	    $k{$tmap{$from}{$_}}++; ## avoid duplicates
	    push @missing, $tmap{$from}{$_} unless $k{$tmap{$from}{$_}}>1;
	}keys(%{$tmap{$from}});
	open(my $fh,">", "$tmp_dir/names.$$.$from") or die "cannot open > $tmp_dir/names.$$.$from: $!";
	local $"="\n";
	print $fh "@missing\n";
	close $fh;
	v_say("Mapping $from..."," ");
	debug("Make_hq: uniprot_parse.pl -Rf $from -t $to_id  -i $tmp_dir/names.$$.$from $flat_file 2>$tmp_dir/missed.$$.$from\n");
#	my $map=`uniprot_parse.pl -f $from -t ID  -i $tmp_dir/names.$$.$from $flat_file 2>$tmp_dir/missed.$$.$from`;
	my $map=`uniprot_parse.pl -vf $from -t $to_id  -i $tmp_dir/names.$$.$from $flat_file `;
	## Read the mapped names into a hash and convert
	my %l;
	map{
	    next if /^\s*$/;
	    my @a=split(/\t/);
	    $l{$a[0]}=$a[1];
	}split("\n",$map);

	## Now build a hash connecting the original (bad) id 
	## to the mapped one
	foreach my $orig_id (keys(%{$tmap{$from}})){
	    $prot_map{Psi2Uni}{$orig_id}=$l{$tmap{$from}{$orig_id}};
	    push @{$prot_map{Uni2Psi}{$prot_map{Psi2Uni}{$orig_id}}},$orig_id;
	}
	## Collect the ids that were NOT mapped
	open($fh, "<", "$tmp_dir/missed.$$.$from");
	my @missed;
	while(<$fh>){
	    /:\s*(.+)/;
	    my $a=$1;
	    @missed=split(/\s+/,$a);
	}
	$missed_prots{$from}=\@missed;
    }
    return(\%prot_map)
}

############################################################
############################################################
sub usage{
    my $us="[options] <psicquic file>";
    my $desc="This script will parse a pscquic interactions file and create a high quality PPI network file. If the psicquic file name that is given is empty, the script will automatically download ALL psicquic interactions from the following online databases: APID, BIND, BioGrid, DIP, I2D, I2D-IMEx, InnateDB, InnateDB-IMEx, IntAct, iRefIndex, MatrixDB, MINT, MolCon, Reactome, Spike, STRING, TopFind\n";

    my %opts=(
	      "usage" => $us,
	      "desc" => $desc,
	      "m" => "Map IDs found in the psicquic file to UniProt IDs (if no file is passed with -M, this is set to true.)",
	      "M" => "Map file linking IDs found in the psicquic file to UniProt IDs",
	      "b" => "Accept only BINARY interactions when building the network.",
	      "B" => "Accept only BINARY AND IP interactions when building the network.",
	      "s" => "The species we are parsing.",
	      "S" => "Spacer. This will be printed at the beginning of error and progress reports. It is useful when this script is part of a larger pipeline",
	      "f" => "The flat file to be parsed for ID mapping. This is only useful in combination with the (-m) option. The script will parse this flat file for any IDs that could not be mapped from UniProt's online service. ",
	      "p" => "Minimum publications mentioning this interaction. Only interactions found in least this many publications (pubmed IDs) will be returned.",
	      "P" => "Minimum high quality publications mentioning this interaction. Only interactions that have been identified by BINARY interaction detection methods in at least this many publications (pubmed IDs) will be returned.",
	      "T" => "Directory to store temporary files in (def:/tmp).",
	      "D" => "List the detection methods (terms) that will be accepted with the current options.",
	      "d" => "Print debugging information to STDERR",
	      "t" => "ID type to map TO. The resulting network will be expressed in this type of ID.",
	      "l" => "Create a log file; For each interaction reported, this log file will contain the supporting publications, source dtabases and IDs that each of the interacting partners were found under in the psiqcuic file. The log file name will be <psicquic filename>.hq.log or whatever is given with option -L.",
	      "L" => "Filename for the log file (def: <psicquic filename>.hq.log).",
	      "h" => "Print this help and exit"
	     );
    print_help_exit(\%opts,0);
}
