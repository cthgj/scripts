#!/usr/bin/perl -w
use strict;
use Getopt::Long qw(:config bundling);
use feature qw(switch say);
# use Smart::Comments;
our $verbose;
our $debug;
our $force;
our $logfile;
require "MY_SUBS.pl";


########################################
# Initialize options to default values #
########################################
my ($help, $binary, $flat_file, $psi_file, $species, $psi2id_passed, $net_file_passed, $dspacer, $missing_file, $missing_GO_file, $GOMAP);
my ($NetUnMapped,$NetMapped,$NetMapped_flat,$NetMapped_uni,$NetMapped_extra)=(0,0,0,0,0);
my $spacer="\t";
my $tmp_dir="/tmp";
my $data_dir=".";
#############################
# Load command line options #
#############################
my @GetOptConf=(
   'flat-file|f=s'=>\$flat_file ,
    'species|s=s' => \$species, 
    'missing-psi|M=s' => \$missing_file , 
    'data-dir|D=s' => \$data_dir , 
    'tmp_dir|T=s' => \$tmp_dir , 
    'missing-go|G=s' => \$missing_GO_file , 
    'id2go-map|g=s' => \$GOMAP , 
    'binary|b' => \$binary , 
    'force|F' => \$force , 
    'help|h' => \$help,
    'verbose|v' => \$verbose,
    'spacer|S' => \$spacer,
    'psi2id|i:s' => \$psi2id_passed , 
    'psi-file|p:s' => \$psi_file , 
   'net-file|n:s' => \$net_file_passed,
   'log-file|l:s' => \$logfile,
    'debug|d' => \$debug 
    );


GetOptions(@GetOptConf) or do {
    exit;
};

$help && do {&usage(@GetOptConf)};
$debug && do {&PrintOptions()};

unless ($flat_file || ($net_file_passed && $psi2id_passed)){
    die "Need a flat file (-f)!\n" unless $flat_file;
}

####################################################
# If no species name has been given, try and guess #
# it from the flat file's name			   #
####################################################
$species||=guess_species($flat_file);

#####################################
# Translate species to default name #
#####################################
$species=get_species($species, "-s");
if($spacer eq "\t"){$dspacer='\t'}

############################################
# Check for mutually incompatible options, #
# files that may be overwritten etc	   #
############################################
my ($previous_map_file, $previous_net_file, $net_file, $psi2id,$LOG)=preFlightCheck();

##################################################
# If a psicquic interactions file has been given #
##################################################
 if(-e $psi_file){
     die("FILE: $psi_file exists but is empty!") if -z $psi_file;

#     p_say(" Parsing psi file: $psi_file", 60, $spacer); 
#     v_say("Counting Psicquic Interactions... ($psi_file) ", $spacer);
#     $TotalPsiInters=`gawk '{print \$1"\\n"\$2}' $psi_file | sort | uniq | wc -l`;
#     $TotalPsiInters=~s/^(\d+).*/$1/s;
 }
###################################################
# If no psicquic interactions file has been given #
# download from psicquic databases		  #
###################################################
else{
    get_PsicquicInteractions($psi_file) unless -e $psi_file;
}

####################################################
# Map all ids in the psicquic file to a UniProt ID #
# unless a map file is given. Generates psi2id     #
####################################################
my ($TotalPsiInters, $TotalPsi,$TotalPsiMapped, $TotalPsiUnMapped)=map_psi2uni($psi_file);


#######################################################
# Parse the psicquic file and build the HQ network    #
#######################################################
my ($NetInters, $NetProts)=build_network($net_file);


#########################################################
#    Build the map file connecting each network file id #
#    to whatever idtype is in this species' go file     #
#########################################################
build_GO_map();
my ($a, $b, $c, $to_id)=build_GO_map();
my %mapped=%{$a};
my %not_mapped=%{$b};
my %names=%{$c};
$NetMapped_uni=scalar keys(%mapped);
# map{v_say("MAPPED: $_ : $mapped{$_}", $spacer);}keys(%mapped);
# map{v_say("NOT_MAPPED: $_ : $not_mapped{$_}");}keys(%not_mapped);

######################################################
# Print any Ids that were not mapped into a tmp_file #
# and map using uniprot_parse.pl.		     #
######################################################
my $file="$tmp_dir/$species.$$.GoNames";
my $mis_file="$tmp_dir/$species.$$.MissingGoNames";
my $unmapped_file="$tmp_dir/$species.$$.UnmappedGoNames";

my @nm=keys(%not_mapped);
debug("MISSED1 : " . scalar @nm , " : @nm\n");

###############################################
# If any IDs were missed, try and find them,  #
# first by parsing the flat file, then by     #
# parsing the extra annotations file	      #
###############################################
if($#nm>=0){
    open(my $fh, ">", "$file") or die " Missing_IDs: Could not open $file\n";
    map{print $fh "$_\n"}keys(%not_mapped);
    close($fh);



#############################################################
# Search for the missing IDs in the flat file using 	    #
# parse_uniprot.pl. IDs that could not be found in 	    #
# the flat file are printed to $mis_file. IDs that 	    #
# were found but could not be mapped are printed to 	    #
# $unmapped_file.					    #
#############################################################
    my $vv="-";
    $vv.="v" if $verbose;
    my $cmd=$vv . "sf ID -t $to_id -i $file -m $mis_file -M $unmapped_file $flat_file";
    my $cmd1= "uniprot_parse.pl -S \"" .  $spacer . "\" " . $cmd; ## this one is run
    my $dcmd="uniprot_parse.pl -S \"" .  $dspacer . "\" " . $cmd; ## this one is printed
    v_say("\tRUNNING(1) $dcmd", $spacer);
    my $map=`$cmd1`;
    #######################################
    # Add the newly mapped ids to %mapped #
    #######################################
    map{
	my @a=split(/\t/);
	$mapped{$a[0]}{$a[1]}++;
	$NetMapped_uni++ unless $mapped{$a[0]}{$a[1]}>1;
    }split(/\n/,$map);

    ###############################################
    # Parse species-specific annotation files to  #
    # get extra ids mapped.			  #
    ###############################################

    my $koko=get_extra_ids($species, $mis_file, $unmapped_file);
    my %extras;
    $koko && do{
	%extras=%{$koko};
	$NetMapped_extra=scalar(keys(%extras));
	v_say("$NetMapped_extra more IDs were mapped from the extra annotations file", $spacer);
    };

    #######################################
    # Add the newly mapped ids to %mapped #
    #######################################
    map{
	foreach my $map_id (keys(%{$extras{$_}})){
	    $mapped{$_}{$map_id}++;
	}
    }keys(%extras);



    ##############################################
    # Print IDs that could not be mapped to a go #
    # annotation name			         #
    ##############################################
    if($missing_GO_file){
	my  $fh=check_file($missing_GO_file,"w") or die "Could not open MissingGoFile : $missing_GO_file : $!\n";
	v_say("Printing IDs with no GO name", $spacer);
	foreach my $map_id (keys(%{$names{ID}})){
	    next if defined($mapped{$map_id});
	    $NetUnMapped++;
	    print $fh "$map_id\n" ;
	}
    }
}## end if($#nm>=0)
############################################
# Print the map file connecting Graph name #
# to go annotation name	into $GOMAP	   #
############################################
## Assign value unless defined($GOMAP)
$GOMAP//="$species.id2go";

my $fh=check_file($GOMAP,"w") or die "Could not open output file : $GOMAP\n";
foreach my $net_id (keys(%mapped)) {    
    $NetMapped++;
    foreach my $map_id (keys(%{$mapped{$net_id}})){
	print $fh "$net_id\t$map_id\n";
    }
}

PrintReport();
exit();

############################# SUBROUTINES ###########################
sub preFlightCheck{
    #################################
    # Check if $tmp_dir is writable #
    #################################
    unless((-d $tmp_dir) && (-w $tmp_dir)) {
	die("ERROR: You do not have write permissions for temp directory $tmp_dir.\n") 
    };

    ###########
    # Logging #
    ###########
    my $LOG=0;
    if(defined($logfile)){
	#print "JOJO ", length($logfile), "\n";
	if ($logfile eq ''){
	    $logfile="$species.$$.log";
	}
	$LOG=check_file($logfile,"w");
    }
    ###############################
    # Check if output files exist #
    ###############################

    $missing_file//="$species.$$.Psi2Uni.missing"; ## assign def val. to $missing_file unless it is defined
    check_file($missing_file,"r") ;
    my ($previous_map_file, $previous_net_file);
    my $psi2id="$species.psi2id";
   ####################################
   # If a map file has been passed as #
   # an option, use it.		      #
   ####################################
    if($psi2id_passed){
	$previous_map_file=$psi2id_passed;
	$psi2id=$psi2id_passed;
	check_file($psi2id,"r");
    }
    ##############################################
    # If no map file was passed as an option and #
    # the default file name exists overwrite it	 #
    ##############################################
    elsif(-e $psi2id){
	print STDERR fold_at_word("\n No psi2id file name was given (-i) and there is already a file with the default name ($psi2id), continuing will overwrite it.", 65,  " "), "\n";
	check_file($psi2id,"w");
    }
    

    my $net_file="$species";
    $binary ? ($net_file.=".binary.gr") : ($net_file.=".gr");
    $net_file_passed && do {$net_file=$net_file_passed };
    if(-e $net_file){
	## If the network file exists, and has been passed as
	## an option, use it.
	if($net_file_passed){
	    $previous_net_file=$net_file;
	    check_file($previous_net_file,"r");

	}
	## If the default network file name exists, and 
	## has NOT been passed as an option, overwrite
	else{check_file($net_file,"w") }
    }
    unless($previous_map_file && $previous_net_file){
	die "Need a psi file (-p) unless both -i and -n are set.\n" unless 
	    $psi_file;
    }
    return($previous_map_file, $previous_net_file, $net_file, $psi2id, $LOG);
}
############################################################
############################################################
sub build_GO_map{
    my $to_id;
    given($species){
	when (/human/){
	    $to_id="ID";
	}
	when (/fly/){
	    $to_id="FlyBase";
	}
	when (/worm/){
	    $to_id="WormBase";
	}
	when (/yeast/){
	    $to_id="SGD";
	}
	when (/mouse/){
	    $to_id="MGI";
	}

    }
#################################################
# Collect the ids to be mapped from the network #
#################################################
    my %names;
    debug("NET_FILE : $net_file");
    my $net;
    if($previous_net_file){
	open( $net, "<", $previous_net_file) or die "Could not open previous net file : $previous_net_file: $!\n";
    }
    else{
	open($net, "<", $net_file) or die "Could not open $net_file\n";
    }
    while(<$net>){
	chomp;
	my @a=split(/\t/);
	########################################################
	# Check that these look like valid UniProt IDs or ACCs #
	########################################################
	if ($a[0]=~/([A-Z\d]{0,8}_[A-Z]{5})/ && $a[1]=~/([A-Z\d]{0,8}_[A-Z]{5})/){
	    $names{ID}{$a[0]}=$names{REV}{$a[0]}=$a[0];
	    $names{ID}{$a[1]}=$names{REV}{$a[1]}=$a[1];
	}
	# elsif((length($a[0])==6 && $a[0]=~/^[A-Z\d]+$/) && 
	# 	  (length($a[1])==6 && $a[1]=~/^[A-Z\d]+$/)) {
	# 	$names{AC}{$a[0]}=$names{REV}{$a[0]}=$a[0];
	# 	$names{AC}{$a[1]}=$names{REV}{$a[1]}=$a[1];
	# }
	else{die("One or both of $a[0] and $a[1] is not a UniProt ID \n");}
    }

    close($net);
###########################################
# Print the network's Ids into a tmp file #
###########################################
# open(my $fh, ">", "$tmp_dir/$species.$$.NetNames") or die "Could not open $tmp_dir/$species.$$.NetNames\n";
# map{print $fh "$_\n"}keys(%names);
# close($fh);
############################
# Map using uniprot_map.pl #
############################
    my ($a, $b)=map_from_UniProt(\%names, $to_id, $species, $tmp_dir, $spacer);
    my %mapped=%{$a};
    my %not_mapped=%{$b};
    
    return($a, $b, \%names, $to_id);
}
############################################################
############################################################
sub map_psi2uni{
    my $psi_file=shift;
    my @mis;
    ###############################################
    # If a map file has been provided, do nothing #
    ###############################################
    if ($previous_map_file && -e $previous_map_file){
	p_say(" Parsing MAP file : $previous_map_file ", 60, $spacer);
	
    }
    else{
	die("FILE: $psi_file is empty!") if -z $psi_file;
	my $cmd="make_psicquic_map.pl -";
	$verbose && do {$cmd.="v"};
	$missing_file//="$species.$$.Psi2Uni.missing";
	$missing_file && do {$cmd.="m $missing_file -"};
	$cmd.="f $flat_file -s $species ";
	my $cmd1=$cmd . "-S \"" .  $spacer . "\" $psi_file"; ## this one is run
	my $dcmd=$cmd . "-S \"" .  $dspacer . "\" $psi_file";## this one is printed
	# v_say(" ");
	p_say(" Mapping psicquic IDs to UniProt IDs ", 60, $spacer);
	debug("CMD1 : $dcmd > $psi2id");
	system("$cmd1 > $psi2id");
	debug("CMD1 : Done");

    }
    ###############
    # Count stuff #
    ###############
    my ($inter_count, $id_count, $map_count, $mis_count)=(0,0,0,0);
    if(-s $psi_file){
	v_say("Counting Psicquic Interactions...($psi_file) ", $spacer,1);
	my $count=`sort $psi_file | uniq`;
	my @aa=split(/\n/, $count);
	my %kk;
	map{
	    chomp();
	    if(/[\w\d]/ && not /^Total:/){
		my @bb=split(/\t/);
		my $pair=join("", sort(@bb));
		$kk{PAIRS}{$pair}++; ## Count unique interactIONS
		$kk{IDs}{$bb[1]}=$kk{IDs}{$bb[0]}=1; ## Count unique interactORS
	    }
	}@aa;
	$inter_count=scalar(keys(%{$kk{PAIRS}}));
	$id_count=scalar(keys(%{$kk{IDs}}));
	v_say(" : $inter_count, $id_count IDs", $spacer);
    }
    ######################################
    # If any IDs were not found	         #
    ######################################
    if(-s $missing_file){
	v_say("Counting missed IDs...($missing_file)", $spacer,1);
#	$mis_count=`sort $missing_file | uniq | wc -l`;
	@mis=split(/\n/, `sort $missing_file | uniq | wc -l`);
	$mis_count=~s/^(\d+).*/$1/s;
	v_say(" : $mis_count", $spacer);
    }
    if(-s $psi2id){
	## Count found IDs
	v_say("Counting found IDs...($psi2id)", $spacer,1);
	$map_count=`gawk '{print \$1}' $psi2id | sort | uniq | wc -l`;
	$map_count=~s/^(\d+).*/$1/s;
	v_say(" : $map_count", $spacer);

    }
    return($inter_count,$id_count,$map_count, $mis_count);
}
############################################################
############################################################
sub build_network{
    my $use_this_net_file=shift();
    if($previous_net_file && -e $previous_net_file){
	p_say(" Parsing NET file : $previous_net_file ", 60, $spacer);
	$use_this_net_file=$previous_net_file;
    }
    else{
	p_say(" Building HQ network ", 60, $spacer);
	my $cmd="make_psicquic_hq.pl -";
	$verbose && do {$cmd.="v"};
	$binary && do {$cmd.="b"};
	$cmd.= "M $psi2id";
	my $cmd1=$cmd . " -S \"" .  $spacer . "\" $psi_file"; ## this one is run
	my $dcmd=$cmd . " -S \"" .  $dspacer . "\" $psi_file"; ## this one is printed
	debug("CMD2 : $dcmd > $use_this_net_file");
	system("$cmd1 > $use_this_net_file");
	debug("CMD2 : Done");
	p_say("", 60, $spacer);
    }
    ## how many interactions are there in the network?
    my $inters=`cat $use_this_net_file | sort | uniq | wc -l`;
    $inters=~s/(\d+).+/$1/s;
    ## how many unique proteins are there in the network?
    my $prots=`gawk '{print \$1"\\n"\$2}' $use_this_net_file | sort | uniq | wc -l`;
    $prots=~s/(\d+).+/$1/s;
    return($inters, $prots)
}
############################################################
############################################################
sub get_extra_ids{
    my $sp=shift;
    my @missing_files=@_;
    my (%not_mapped, %hash);
    ########################
    # Read the missing IDs #
    ########################
    foreach my $f (@missing_files){
	#foreach my $f ("a.aa", "b.aa"){
	debug("EXTRA opening $f\n");
	if(-e $f){
	    open(my $fh, "<", "$f");
	    while(<$fh>){
		chomp;
		$not_mapped{$_}++;
	    }
	}
    }
    my @missing=keys(%not_mapped);
    return(0) unless $#missing>=0;
    my ($file, $TO_ID, $map, $unmapped_file, $mis_file, $wanted_id_pos,$sep);
    $unmapped_file="$tmp_dir/$species.$$.UnmappedExtra";
    $mis_file="$tmp_dir/$species.$$.MissingExtra";
    debug("EXTRA will look for : @missing in $file\n");
    given($sp){
	when ($_ eq 'fly'){
	    ## http://flybase.org/static_pages/downloads/FB2012_02/genes/fbgn_NAseq_Uniprot_fb_2012_02.tsv.gz
	    $file=$data_dir . "/fly.extra";
	    $TO_ID="GeneName";
	    $wanted_id_pos=1;
	    $sep="\t";
	}
	when (/human/){
	    $file=$data_dir . "/";
	    $TO_ID="ID";
	}
	when (/worm/){
	    ## ftp://ftp.wormbase.org/pub/wormbase/species/c_elegans/annotation/geneIDs/c_elegans.current.geneIDs.txt.gz
	    $file=$data_dir . "/worm.extra";
	    $TO_ID="WormBase";
	    $wanted_id_pos=0;
	    $sep=",";

	}
	when (/yeast/){
	    $file=$data_dir . "/yeast.extra";
	    $TO_ID="SGD";
	}
	when (/mouse/){
	    $file=$data_dir . "/mouse.extra";
	    $TO_ID="MGI";
	}
    }
    debug("extra ID type: $TO_ID in file : $file");
    p_say(" Searching for extra IDs ", 60, $spacer);
    ################################################
    # Map the ids to whatever format is present in #
    # the extra annotations file		   #
    ################################################
    if ($TO_ID eq 'ID'){die("SHIIIIIIIIIIIT, IMPLEMENT!\n");}
    else{
	my $tmp_file="$tmp_dir/$species.$$.NamesExtra";
	open(my $fh, ">", "$tmp_file") or die "extra_ids: Could not open $tmp_file\n";
	map{print $fh "$_\n"}@missing;
	close($fh);
        my $vv="-";
	$vv.="v" if $verbose;
	################################################################
        # The -a flag will cause uniprot_parse.pl to return ALL	       #
	# of the matching ID types, synonyms included.		       #
        ################################################################
	my $cmd=$vv . "asf ID -t $TO_ID -i $tmp_file";
	$cmd.=" -m $mis_file -M $unmapped_file $flat_file" ;
	v_say("\tRUNNING(2) uniprot_parse.pl $cmd", $spacer);
    	$cmd="uniprot_parse.pl -S \"\t\" " . $cmd;
	$map=`$cmd`;
	v_say(" ");
    }
    ########################################################
    # Read the mapped names into a hash and convert        #
    ########################################################
    my (%MAP, %WANT);
    my @found_lines=split(/\n/, $map);
    foreach my $line (@found_lines){
	chomp($line);
	my  @a=split(/\t/, $line);
	my $id=shift(@a);
	#################################################
	# There may be one-to-many relationships	#
	# $MAP{$ID}{$to_id}++			        #
	#################################################
	map{
	    $MAP{$id}{$_}++; ## keys(%$MAP{$id} are all the IDs synonyms
	    $WANT{$_}=$id;   ## This synonym maps back to $id
	}@a;
    }
    ##################################################
    # Look for the IDs in the extra annotations file #
    ##################################################
    open(my $fh, "<", $file) or die "extra_ids: Could not open $file\n";
    while(<$fh>){
	next if /^\s*\#/;
	chomp;
	my @a=split(/$sep/);
	foreach my $id (@a){
	    defined($WANT{$id}) && do {
		###############################################################
		# $a[$wanted_id_pos] is the id type we want to extract	      #
		# from the extra annotations file.			      #
		###############################################################
		$hash{$WANT{$id}}{$a[$wanted_id_pos]}++;
	    }
	}
    }
    scalar keys %hash && do {return(\%hash)};
    return(0);
}
############################################################
############################################################

sub PrintOptions(){
    p_say(" OPTIONS ", 60, $spacer);
    my  @a=qw($help $binary $flat_file $species $force $psi2id_passed $net_file_passed $dspacer $missing_file $missing_GO_file $spacer $tmp_dir $data_dir);
    my @b=($help, $binary, $flat_file, $species, $force, $psi2id_passed, $net_file_passed, $dspacer, $missing_file, $missing_GO_file, "-". $spacer. "-", $tmp_dir, $data_dir);
    for my $i (0..$#a){
	my $ll=16-length($a[$i]);
	my $kk=$a[$i] . " " x $ll;
	debug("$kk\t:\t$b[$i]\n") if $b[$i]
    }
    psi2id_passed_say("", 60, $spacer);
    
}

sub get_PsicquicInteractions{
    p_say(" Downloading psicquic interactions ", 60, $spacer);
    my $psi_out=shift();
    ## Check if the output file exists
    check_file($psi_out,"w") ;
    my (@files, @e_files);
    my $prog_file="$tmp_dir/psi.$$.progress";
    system("rm $prog_file") if -e $prog_file ;
    my $db_count=`get_psicquic_interactome.pl -c` or die "Could not count databases: $!\n";
    chomp ($db_count);
    for my $i (0..$db_count-1){
	my $ii=$i+1; ## get_psicquic_interactome.pl takes numes from 1 to $db_count-1
	my $tmp_psi_file=$tmp_dir . "/$$.$ii.psi";
	my $tmp_error_file=$tmp_dir . "/$$.$ii.psi.er";
	$files[$i]=$tmp_psi_file;
	$e_files[$i]=$tmp_error_file;
	my $cmd="get_psicquic_interactome.pl -d $ii > $tmp_psi_file 2>$tmp_error_file && echo \"$ii Done\" >> $prog_file &";
	debug($cmd);
	system($cmd);

   }
    my $finished=0;
    my %hash;
    
    #######################################################
    # $finished will == $db_count when all		  #
    # get_psicquic_interactome.pl instances are finished. #
    #######################################################
    while($finished<$db_count){
	my $done_files;
	-e $prog_file && do {
	    $finished=`wc -l $prog_file`;
	    $finished=~s/(\d+).+/$1/s;
	};
	sleep(2);
	v_say("DOWNLOADED: $finished  of  $db_count\r", $spacer,1);
    }
    v_say("DOWNLOADED: $finished  of  $db_count", $spacer);
    # v_say("Counting Psicquic IDs... $psi_out ", $spacer);
     system("cat @files | grep -v \"^Total:\"> $psi_out");
    # my $cc=`wc -l $psi_out`;
    # $cc=~s/^(\d+).*/$1/s;
    # chomp($cc);
    # return($cc);
}

sub l_say{
    return(0) unless $LOG;
    print $LOG $_, "\n";
}


sub PrintReport{
    if($TotalPsiInters==0){
	$TotalPsiInters=`gawk '{print \$1"\\n"\$2}' $psi_file | sort | uniq | wc -l`;
	$TotalPsiInters=~s/^(\d+).*/$1/s;
    }
    my $rep= <<EndOfMsg;

$spacer---------------------- REPORT ------------------------------
$spacer Psicquic interactions                : $TotalPsiInters, $TotalPsi IDs
$spacer Psicquic IDs mapped/not mapped       : $TotalPsiMapped / $TotalPsiUnMapped
$spacer Interactions in network file         : $NetInters
$spacer Proteins in network file             : $NetProts
$spacer Network IDs mapped to GO_annot IDs   : $NetMapped
$spacer                             UniProt  : $NetMapped_uni
$spacer                           Flat File  : $NetMapped_flat
$spacer                          Extra File  : $NetMapped_extra
$spacer Network proteins with no GO_annot ID : $NetUnMapped
$spacer------------------------------------------------------------
EndOfMsg
v_say($rep);
}
sub usage{
    my $us="[options] <psicquic file>";
    my $desc="This script will take a psicquic file and create a map file\n";
    $desc .="with the uniprot accession of each psicquic id.blah blah";
    my %opts=(
	"usage" => $us,
	"desc" => $desc,
	"f|flat-file" => "Flat file to be used for mapping.",
	"M|missing-psi" => "Missing psicquic file. Psicquic IDs from the input file that could not be linked to a UniProt ID will be printed to this file.",
	"missing-go|G" => "Missing GOfileIDs file. Uniprot IDs from the generated network file that could not be linked to an ID of whatever format used in the GO annotations file of this species, will be printed to this file.",
	"binary|b" => "Accept only BINARY interactions when building the network.",
	"species|s" => "The species we are parsing. If none is given the value is guessed from the flat file",
	"net-file|n" => "Network File. If this file exists, do not parse the psicquic file and create a network this network instead. If it doesn't, the psicquic file will be parsed and the resulting network printed into this file. If no file is given, the psicquic file will be parsed and the network will be printed into <SPECIES_NAME>.<NETWORK_TYPE>.gr.",
	"psi2id|i" => "psi2id map file. If this file exists, do not parse/convert ids from the psicquic input file use this map instead. If it doesn't, the psicquic file will be parsed and the map printed into this file. If no file is given, the psicquic file will be parsed and the map will be printed into <SPECIES_NAME>.psi2id.",
	"id2go|g" => "id2go map file connecting each network file id to whatever idtype is in this species' go file. If no file is given, map will be printed into <SPECIES_NAME>.id2go.",
	"psi-file|p" => "File containing psicquic interaction data. If this file exists, psicquic interactions will be extracted from it. If it doesn't, remote databases will be queried and the interactions downloaded into this file.",
	"log-file|l" => "Create a log file.",
	"data-dir|D" => "Data directory (def: ./)",
	"tmp-dir|T" => "Temp. files will be created here (def: ./tmp)",
	"debug|d" => "Print debugging information to STDERR",
	"force|F" => "Force, overwrite any existing files without prompting",
	# "" => "",
	# "" => "",
	# "" => "",
	# "" => "",
	"help|h" => "Print this help and exit"
	);
    print_help_exit(\%opts,0);
}
