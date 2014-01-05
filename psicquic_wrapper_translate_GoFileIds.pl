#!/usr/bin/perl -w

#####################################################################
# This is the old version of the script, before I found the 	    #
# UniProt GO annotation files. Since these are given with	    #
# UniProt IDs and not whatever was in the original annotation 	    #
# from GO, I no longer need to translate UniProt IDs to GOFile IDs. #
# I am keeping this version just in case.			    #
#####################################################################





use strict;
use Getopt::Long qw(:config bundling);
use feature qw(switch say);
our $verbose;
our $debug;
our $force;
my %files; ## this hash will hold the names of all files created
sub v_say;
require "MY_SUBS.pl";
# unless ($flat_file || ($net_file_passed && $psi2id_passed)){
#     die "Need a flat file (-f)!\n" unless $flat_file;
# }

######################################################
# Parse options, and check for mutually incompatible #
# options, files that may be overwritten etc	     #
######################################################
my ($TEST, $flat_file,$species,$missing_psi_file,$missing_GO_file,$GOMAP,$binary,$help,$psi2id,$psi_file,$get_psi_file,$net_file,$generate_map,$generate_network,$LOG, $tmp_dir, $dspacer, $data_dir, $spacer  )=preFlightCheck();


###################################################
# If no psicquic interactions file has been given #
# download from psicquic databases		  #
###################################################
get_PsicquicInteractions($psi_file);


####################################################
# Map all ids in the psicquic file to a UniProt ID #
# unless a map file is given. Generates psi2id     #
####################################################
my ($TotalPsiInters, $TotalPsi,$TotalPsiMapped, $TotalPsiUnMapped, $kk)=map_psi2uni($psi_file, $psi2id, $missing_psi_file);
my %mapped=%{$kk};
#######################################################
# Parse the psicquic file and build the HQ network    #
#######################################################
my ($NetInters, $NetProts, $ll)=build_network($net_file);
my %network=%{$ll};


PrintReport();
exit();


########################################### OLD START ########################################

###################################################################
# Will use GO annotations from QuickGo, already in uniprot format #
#     therefore, I do not need these steps anymore                #
###################################################################

# ######################################################
# # Build the map file connecting each network file id #
# # to whatever idtype is in this species' go file     #
# ######################################################
# my ($a, $b, $c, $to_id)=map_uni2goid(\%network, \%mapped);
# %mapped=%{$a};
# my %not_mapped=%{$b};
# my %names=%{$c};

# ######################################################
# # Print any Ids that were not mapped into a tmp_file #
# # and map using uniprot_parse.pl.		     #
# ######################################################
# my $file="$tmp_dir/$species.$$.GoNames";
# my $mis_file="$tmp_dir/$species.$$.MissingGoNames";
# my $unmapped_file="$tmp_dir/$species.$$.UnmappedGoNames";

# my @nm=keys(%not_mapped);
# debug("MISSED1 : " . scalar @nm , " : @nm\n");

# ###############################################
# # If any IDs were missed, try and find them,  #
# # first by parsing the flat file, then by     #
# # parsing the extra annotations file	      #
# ###############################################
# my  ($aa, $bb, $NetMapped_extra, $NetUnMapped, $NetMapped, $NetMapped_flat)=0;
# if ( $#nm >= 0) {
#     ($aa, $bb, $NetMapped_extra, $NetUnMapped, $NetMapped, $NetMapped_flat)=get_missing_ids($file, \%mapped, \%not_mapped) if $#nm>=0;
#     %mapped=%{$aa};
#     %not_mapped=%{$bb};
# }


# ############################################
# # Print the map file connecting Graph name #
# # to go annotation name	into $GOMAP	   #
# ############################################
# ## Assign value unless defined($GOMAP)
# $GOMAP//="$species.id2go";
# $files{NAMES}{'gomap'}=$GOMAP;
# $files{CREATED}{'gomap'}="Existing";

# my $fh=check_file($GOMAP,"w") or die "Could not open output file : $GOMAP\n";
# foreach my $net_id (keys(%mapped)) {    
#     $NetMapped++;
#     foreach my $map_id (keys(%{$mapped{$net_id}})) {
# 	print $fh "$net_id\t$map_id\n";        
#     }
# }
# my $NetMapped_uni=scalar keys(%mapped);

# #############################################
# # Create a network containing ONLY IDs that #
# # were mapped to a gofile-type ID           #
# #############################################
# my $NET=check_file($net_file,"r");
# my $gonet_file=$net_file;
# $gonet_file=~s/\.gr/.gonet.gr/;
# my $GO_NET=check_file($gonet_file,"w");
# $files{'NAMES'}{'GOnet'}=$gonet_file;
# $files{'CREATED'}{'GOnet'}='Created';

# v_say("Creating GO Network file: $gonet_file\n", $spacer);
# while (<$NET>) {
#     chomp;
#     my @a=split(/\t/);
#     if (defined($mapped{$a[0]}) && defined($mapped{$a[1]}) ) {
#         print $GO_NET "$_\n"
#     } 
# }

########################################### OLD END ########################################



#####################################################################
############################# SUBROUTINES ###########################
#####################################################################
sub preFlightCheck
{
    ########################################
    # Initialize options to default values #
    ########################################
    my ($TEST,$flat_file,$species,$missing_psi_file,$missing_GO_file,$GOMAP,$binary,$help,$psi2id,$psi_file,$net_file);
    my $spacer="\t";
    my $dspacer=$spacer;
    ## Useful for debugging 
    if ($spacer eq "\t") {
	$dspacer='\t';
    }
    $tmp_dir="/tmp";
    my $data_dir=".";
    #############################
    # Load command line options #
    #############################
    my @GetOptConf=(
		    'flat-file|f=s'=>\$flat_file ,
		    'species|s=s' => \$species, 
		    'missing-psi|M=s' => \$missing_psi_file , 
		    'data-dir|D=s' => \$data_dir , 
		    'tmp_dir|T=s' => \$tmp_dir , 
		    'missing-go|G=s' => \$missing_GO_file , 
		    'id2go-map|g=s' => \$GOMAP , 
		    'binary|b' => \$binary , 
		    'force|F' => \$force , 
		    'help|h' => \$help,
		    'verbose|v' => \$verbose,
		    'spacer|S' => \$spacer,
		    'psi2id|i:s' => \$psi2id , 
		    'psi-file|p:s' => \$psi_file , 
		    'net-file|n:s' => \$net_file,
		    'debug|d' => \$debug,
                    'test|t'=>\$TEST
		   );
    GetOptions(@GetOptConf) or do {
	exit;
    };
    $help && do {usage(@GetOptConf)};
    
    ####################################################
    # If no species name has been given, try and guess #
    # it from the flat file's name		       #
    ####################################################
    $species//=guess_species($flat_file);
    
    #####################################
    # Translate species to default name #
    #####################################
    $species=get_species($species, "-s");
    #################################
    # Check if $tmp_dir is writable #
    #################################
    unless((-d $tmp_dir) && (-w $tmp_dir)) {
	die("ERROR: You do not have write permissions for temp directory $tmp_dir.\n") 
    };
    ##########################################
    # Check if a psi file has been given, if #
    # not, we will need to download it.	     #
    ##########################################
    my $get_psi_file=0;
    if ($psi_file) {
	-z $psi_file && do {die("FILE: $psi_file exists but is empty!")};
	##################################################
        # Download from psicquic databases unless	 #
	# $psi_file exists				 #
        ##################################################
	$get_psi_file=1 unless -e $psi_file;
    }
    ######################################################
    # If no psi file has been given, download it	 #
    ######################################################
    else {
	$get_psi_file=1;
    }
    ###############################
    # Check if output files exist #
    ###############################
    $missing_psi_file//="$species.$$.Psi2Uni.missing"; ## assign def val. to $missing_file unless it is defined
    check_file($missing_psi_file,"r") ;

    ####################################
    # If a map file has been passed as #
    # an option, use it.		      #
    ####################################
    my $generate_map=1;
    if ($psi2id) {
	## if the file exists and is NOT empty
	if (-s $psi2id) {
	    check_file($psi2id,"r");
	    $generate_map=0;	## Do not rebuild psi2id map
	}
	## if the IS empty or nonexistant
	else {
	    check_file($psi2id,"w");	    
	} 
    }
    ##############################################
    # If no map file was passed as an option ,   #
    # use the default file name          	 #
    ##############################################
    else {
	$psi2id="$species.psi2id";
	if (-e $psi2id) {
	    print STDERR fold_at_word("\n No psi2id file name was given (-i) and there is already a file with the default name ($psi2id), continuing will overwrite it.", 65,  " "), "\n";
	}
	check_file($psi2id,"w");
    }
    ########################
    # Get network filename #
    ########################
    my $generate_network=1;
    if ($net_file) {
	## if the file exists and is NOT empty
	if (-s $net_file) {
	    check_file($net_file,"r");
	    $generate_network=0; ## Do not rebuild network
	}
	## if the IS empty or nonexistant
	else {
	    check_file($net_file,"w");
	}     
    }
    ##############################################
    # If no net file was passed as an option ,   #
    # use the default file name          	 #
    ##############################################
    else {
	my $net_file="$species";
	$binary ? ($net_file.=".binary.gr") : ($net_file.=".gr");
	if (-s $net_file) {
	    ## If the network file exists and is NOT empty, warn
	    print STDERR fold_at_word("\n No network file name was given (-n) and there is already a file with the default name ($net_file), continuing will overwrite it.", 65,  " "), "\n";
	}
	check_file($net_file,"w");
    }
    my @options=($TEST,$flat_file,$species,$missing_psi_file,$missing_GO_file,$GOMAP,$binary,$help,$psi2id,$psi_file,$get_psi_file,$net_file,$generate_map,$generate_network,$LOG, $tmp_dir, $dspacer, $data_dir, $spacer);
    $debug && do {PrintOptions(@options)};
    return(@options);
}



############################################################
############################################################
sub get_PsicquicInteractions{
    if ($get_psi_file==0) {
        p_say(" Parsing existing psicquic file : $psi_file ", 60, $spacer);
        $files{NAMES}{'psi'}=$psi_file;
        $files{CREATED}{'psi'}='Existing';
    }
    else {
        p_say(" Downloading psicquic interactions ", 60, $spacer);
        my $psi_out=shift();
        ## Check_File if the output file exists
        check_file($psi_out,"w") ;
        my (@files, @e_files);
        my $prog_file="$tmp_dir/psi.$$.progress";
        system("rm $prog_file") if -e $prog_file ;
        my $db_count=`get_psicquic_interactome.pl -c` or die "Could not count databases: $!\n";
        my $db_names=`get_psicquic_interactome.pl -L` or die "Could not list databases: $!\n";
        chomp($db_names);
        my %dbs=map{split(/\s+:\s+/)} split(/\n/, $db_names); 
        chomp ($db_count);
        for my $i (0..$db_count-1) {
            my $ii=$i+1; ## get_psicquic_interactome.pl takes numes from 1 to $db_count-1
            my $tmp_psi_file=$tmp_dir . "/$$.$dbs{$ii}.psi";
            my $tmp_error_file=$tmp_dir . "/$$.$dbs{$ii}.psi.er";
            $files[$i]=$tmp_psi_file;
            $e_files[$i]=$tmp_error_file;
            my $cmd="get_psicquic_interactome.pl -d $ii > $tmp_psi_file 2>$tmp_error_file && echo \"$ii Done\" >> $prog_file &";
            debug($cmd);
             system($cmd) unless $TEST ;
        }
        my $finished=0;
        $finished=$db_count if $TEST;        
        my %hash;
        #######################################################
        # $finished will == $db_count when all		  #
        # get_psicquic_interactome.pl instances are finished. #
        #######################################################
        while ($finished<$db_count) {
            my $done_files=$finished;
            -e $prog_file && do {
                $finished=`wc -l $prog_file`;
                $finished=~s/(\d+).+/$1/s;
            };
            sleep(2);
            $done_files!=$finished && do {
                v_say("DOWNLOADED: $finished  of  $db_count\r", $spacer,1);
            }
        }
        v_say("DOWNLOADED: $finished  of  $db_count", $spacer);
        $files{NAMES}{'psi'}=$psi_out;
        $files{CREATED}{'psi'}="Created";
        $TEST ?
            v_say("TESTING: cat @files | grep -v \"^Total:\"> $psi_out") :
            system("cat @files | grep -v \"^Total:\"> $psi_out");
    }
}


############################################################
############################################################
sub map_psi2uni{
    my $psi_file=shift;
    my $psi2id=shift;
    my $missing_file=shift;
    my @mis;
    ###############################################
    # If a map file has been provided, do nothing #
    ###############################################
    if ($generate_map==0) {
	p_say(" Parsing MAP file : $psi2id ", 60, $spacer);
        $files{NAMES}{'psi2id'}=$psi2id;
        $files{CREATED}{'psi2id'}='Created';
    } 
    else {
	die("FILE: $psi_file is empty!") if -z $psi_file;
	my $cmd="make_psicquic_map.pl -";
	$verbose && do {$cmd.="v"};
	$missing_file//="$species.$$.Psi2Uni.missing";
	$missing_file && do {$cmd.="m $missing_file -"};
	$cmd.="f $flat_file -s $species ";
	my $cmd1=$cmd . "-S \"" .  $spacer . "\" $psi_file"; ## this one is run
	my $dcmd=$cmd . "-S \"" .  $dspacer . "\" $psi_file"; ## this one is printed
	# v_say(" ");
	p_say(" Mapping psicquic IDs to UniProt IDs ", 60, $spacer);
	debug("CMD1 : $dcmd > $psi2id");
	system("$cmd1 > $psi2id") unless $TEST;
        debug("CMD1 : Done");
        $files{NAMES}{'psi2id'}=$psi2id;
        $files{CREATED}{'psi2id'}="Existing";

    }
    ###############
    # Count stuff #
    ###############
    my ($inter_count, $id_count, $map_count, $mis_count)=(0,0,0,0);
    if (-s $psi_file) {
	v_say("Counting Psicquic Interactions...($psi_file) ", $spacer,1);
	my $fh=check_file($psi_file, "r");
	my %kk;
	my $taxid=get_species($species, "-s", 1);
	while (<$fh>) {
	    chomp;
	    next unless /taxid:$taxid/;
	    my @bb=split(/\t/);
	    my $pair=join("", sort($bb[0],$bb[1]));
	    $kk{PAIRS}{$pair}++; ## Count unique interactIONS
	    $kk{IDs}{$bb[1]}=$kk{IDs}{$bb[0]}=1; ## Count unique interactORS
	}
	$inter_count=scalar(keys(%{$kk{PAIRS}}));
	$id_count=scalar(keys(%{$kk{IDs}}));
	v_say(" : $inter_count, $id_count IDs", $spacer);
    }
    ######################################
    # If any IDs were not found	         #
    ######################################
    if (-s $missing_file) {
	v_say("Counting missed IDs...($missing_file)", $spacer,1);
	#	$mis_count=`sort $missing_file | uniq | wc -l`;
	@mis=split(/\n/, `sort $missing_file | uniq | wc -l`);
	$mis_count=~s/^(\d+).*/$1/s;
	v_say(" : $mis_count", $spacer);
    }
    my %kk;
    if (-s $psi2id) {
	## Count found IDs
	v_say("Counting found IDs...($psi2id)", $spacer,1);
	my $fh=check_file($psi2id, "r");
        while (<$fh>) {
	    chomp;
	    my @bb=split(/\t/);
	    push @{$kk{$bb[0]}{'UniID'}}, $bb[1];
            push @{$kk{$bb[1]}{'PsiID'}}, $bb[0];
	}
	$map_count=scalar keys %kk;
	v_say(" : $map_count", $spacer);
    }
    return($inter_count,$id_count,$map_count, $id_count-$map_count, \%kk);
}
############################################################
############################################################
sub build_network{
    my $net_file=shift();
    if ($generate_network==0) {
	p_say("Parsing NET file : $net_file ", 60, $spacer);
        $files{NAMES}{'net'}=$net_file;
        $files{CREATED}{'net'}='Existing';
    } 
    else {
	p_say(" Building HQ network ", 60, $spacer);
	my $cmd="make_psicquic_hq.pl -";
	$verbose && do {$cmd.="v"};
	$binary && do {$cmd.="b"};
	$cmd.= "M $psi2id -s $species";
	my $cmd1=$cmd . " -S \"" .  $spacer . "\" $psi_file"; ## this one is run
	my $dcmd=$cmd . " -S \"" .  $dspacer . "\" $psi_file"; ## this one is printed
	debug("CMD2 : $dcmd > $net_file");
	system("$cmd1 > $net_file");
	debug("CMD2 : Done");
	p_say("", 60, $spacer);
        $files{NAMES}{'net'}=$net_file;
        $files{CREATED}{'net'}='Created';
    }
    my $nn=check_file($net_file, "r");
    my (%h, %net);
    my $inters=0;
    while (<$nn>) {
	chomp;
	next unless  /\t/;
	my @a=split(/\t/);
	my $pair=join("\t", sort(@a));
	$h{$a[0]}=$h{$a[1]}=1;
        $net{$pair}++;
    }
    ## how many unique proteins are there in the network?
    my $prots=scalar keys %h;
    $inters=scalar(%net);
    return($inters, $prots, \%net)
}
############################################################
############################################################
sub map_uni2goid{
    my %net=%{shift()};
    my %map=%{shift()};
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
    # open($net, "<", $net_file) or die "Could not open $net_file\n";
    # while (<$net>) {
    foreach my $pair (keys (%net)) {
	my @a=split(/\t/, $pair);
	################## ####################################
	# Check that these look like valid UniProt IDs or ACCs #
	########################################################
	if ($a[0]=~/([A-Z\d]{0,8}_[A-Z]{5})/ && $a[1]=~/([A-Z\d]{0,8}_[A-Z]{5})/) {
           
            #$map{ID}{${$map{$a[0]}{'UniID'}}[0]}=$a[0];
            $names{ID}{$a[0]}=$names{REV}{$a[0]}=$a[0];
	    $names{ID}{$a[1]}=$names{REV}{$a[1]}=$a[1];
	}
	# elsif((length($a[0])==6 && $a[0]=~/^[A-Z\d]+$/) && 
	# 	  (length($a[1])==6 && $a[1]=~/^[A-Z\d]+$/)) {
	# 	$names{AC}{$a[0]}=$names{REV}{$a[0]}=$a[0];
	# 	$names{AC}{$a[1]}=$names{REV}{$a[1]}=$a[1];
	# }
	else {
	    die("One or both of $a[0] and $a[1] is not a UniProt ID \n");
	}
    }

    ############################
    # Map using uniprot_map.pl #
    ############################
    my ($a, $b)=map_from_UniProt(\%names, $to_id, $species, $tmp_dir, $spacer);
    return($a, $b, \%names, $to_id);
};
############################################################
############################################################
sub get_extra_ids{
    my $sp=shift;
    my @missing_files=@_;
    my (%not_mapped, %hash);
    ########################
    # Read the missing IDs #
    ########################
    foreach my $f (@missing_files) {
	#foreach my $f ("a.aa", "b.aa"){
	debug("EXTRA opening $f\n");
	if (-e $f) {
	    open(my $fh, "<", "$f");
	    while (<$fh>) {
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
    given($sp){
	when ($_ eq 'fly'){
	    ## http://flybase.org/static_pages/downloads/FB2012_02/genes/fbgn_NAseq_Uniprot_fb_2012_02.tsv.gzip64
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
    debug("extra ID type: $TO_ID, in file : $file");
    p_say(" Searching for extra IDs ", 60, $spacer);
    ################################################
    # Map the ids to whatever format is present in #
    # the extra annotations file		   #
    ################################################
    if ($TO_ID eq 'ID') {
	die("SHIIIIIIIIIIIT, IMPLEMENT!\n");
    }
    else {
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
    foreach my $line (@found_lines) {
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
    while (<$fh>) {
	next if /^\s*\#/;
	chomp;
	my @a=split(/$sep/);
	foreach my $id (@a) {
	    defined($WANT{$id}) && do {
		################################################
		# $a[$wanted_id_pos] is the id type we want    #
		# to extract from the extra annotations file.  #
		################################################
		$hash{$WANT{$id}}{$a[$wanted_id_pos]}++;
	    }
	}
    }
    scalar keys %hash && do {return(\%hash)};
    return(0);
}
############################################################
############################################################
sub get_missing_ids{
    my ($file, $aa, $bb)=@_;
    my %mapped=%{$aa};
    my %not_mapped=%{$aa};
    my ($NetMapped_extra, $NetUnMapped, $NetMapped, $NetMapped_flat)=0;
    
    open(my $fh, ">", "$file") or die " Missing_IDs: Could not open $file\n";
    map{print $fh "$_\n"}keys(%not_mapped);
    close($fh);
    #########################################################
    # Search for the missing IDs in the flat file using     #
    # parse_uniprot.pl. IDs that could not be found in 	    #
    # the flat file are printed to $mis_file. IDs that 	    #
    # were found but could not be mapped are printed to     #
    # $unmapped_file.					    #
    #########################################################
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
    my $uni_mapped=scalar keys %mapped;
    map{
	my @a=split(/\t/);
	$mapped{$a[0]}{$a[1]}++;
	$NetMapped_uni++ unless $mapped{$a[0]}{$a[1]}>1;
    }split(/\n/,$map);
    ## Count the ids that were newly mapped from the flat file
    $NetMapped_flat = (scalar (keys(%mapped))) - $uni_mapped;
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
	foreach my $map_id (keys(%{$extras{$_}})) {
	    $mapped{$_}{$map_id}++;
	}
    }keys(%extras);
    ##############################################
    # Print IDs that could not be mapped to a go #
    # annotation name			         #
    ##############################################
    if ($missing_GO_file) {
	my  $fh=check_file($missing_GO_file,"w") or die "Could not open MissingGoFile : $missing_GO_file : $!\n";
	v_say("Printing IDs with no GO name", $spacer);
	foreach my $map_id (keys(%{$names{ID}})) {
	    next if defined($mapped{$map_id});
	    $NetUnMapped++;
	    print $fh "$map_id\n" ;
	}
    }
return(\%mapped, \%not_mapped, $NetMapped_extra, $NetUnMapped, $NetMapped, $NetMapped_flat);
}


############################################################
############################################################
sub PrintReport{
    if ($TotalPsiInters==0) {
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
$spacer------------------- FILES CREATED --------------------------
$spacer Psicquic Interactions file           : $files{NAMES}{psi} ($files{CREATED}{psi})
$spacer Psicquic ID to Uniprot ID map file   : $files{NAMES}{psi2id} ($files{CREATED}{psi2id})
$spacer Uniprot ID to GO-file ID map file    : $files{NAMES}{gomap} ($files{CREATED}{gomap})
$spacer Network file                         : $files{NAMES}{net} ($files{CREATED}{net})
$spacer GO Network file                      : $files{NAMES}{GOnet} ($files{CREATED}{GOnet})
EndOfMsg
    v_say($rep);
}

############################################################
############################################################
sub PrintOptions{
    p_say(" OPTIONS ", 60, $_[$#_]);
    my  @a=qw($flat_file $species $missing_psi_file $missing_GO_file $GOMAP $binary $help $psi2id $psi_file $get_psi_file $net_file $generate_map $generate_network $LOG  $tmp_dir  $dspacer $data_dir $spacer );

    for my $i (0..$#a) {
        my $ll=16-length($a[$i]);
    	my $kk=$a[$i] . " " x $ll;
	my $lo=$_[$i];
	$lo="-". $lo . "-" if $a[$i] eq '$spacer';
    	debug("$kk\t:\t$lo\n") if defined($_[$i])
    }
    p_say("", 60, $_[$#_]);
    
}

############################################################
############################################################
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
              "test|t" => "Test run",
	      # "" => "",
	      # "" => "",
	      # "" => "",
	      # "" => "",
	      "help|h" => "Print this help and exit"
	     );
    print_help_exit(\%opts,0);
}
