#!/usr/bin/perl -w
use strict;
use Getopt::Long qw(:config bundling);
use feature qw(switch say);
use File::Path;
our $verbose;
our $debug;
our $force;
my %files; ## this hash will hold the names of all files created
sub v_say;
require "MY_SUBS.pl";

usage() unless $ARGV[0];
######################################################
# Parse options, and check for mutually incompatible #
# options, files that may be overwritten etc	     #
######################################################
my ($TEST, $flat_file,$species,$missing_psi_file,$binary,$help,$psi2id,$psi_file,$get_psi_file,$net_file,$generate_map,$generate_network,$LOG, $tmp_dir, $dspacer, $data_dir, $dodablast, $blast_file, $get_all_five, $get_all_species, $spacer  )=preFlightCheck();


###################################################
# If no psicquic interactions file has been given #
# download from psicquic databases		  #
###################################################
get_PsicquicInteractions($psi_file);


#########################################################
# Map all ids in the psicquic file to a UniProt ID      #
# unless a map file is given. Generates $psi2id file    #
#########################################################
my ($TotalPsiInters, $TotalPsi,$TotalPsiMapped, $TotalPsiUnMapped, $kk)=map_psi2uni($psi_file, $psi2id, $missing_psi_file);
my %mapped=%{$kk};



#######################################################
# Parse the psicquic file and build the HQ network    #
#######################################################
my ($NetInters, $NetProts, $ll)=build_network($net_file);
my %network=%{$ll};


#################################################################
# Attempt to remove redundancy. Retrieve FASTA seqs for         #
# each of the sequences, blast against each other and 	        #
# remove those those sequences with 100%id, where 	        #
# query!=sbjct and where the entire sequence of the subject     #
# is aligned perfectly against the query for its ENTIRE length. #
#################################################################
my ($nr_NetInters, $nr_NetProts)=(0,0);
$dodablast && do {($nr_NetInters, $nr_NetProts)=run_network_blast($net_file,$flat_file, $blast_file)};

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
    my ($TEST,$species,$missing_psi_file,$binary,$help,$psi2id,$psi_file,$net_file, $dodablast, $blast_file, $get_all_five, $get_all_species);
    my $spacer="\t";
    my $dspacer=$spacer;
    ## Useful for debugging 
    if ($spacer eq "\t") {
	$dspacer='\t';
    }
    $tmp_dir=".psi_tmp";
    my $data_dir="./data";
    #############################
    # Load command line options #
    #############################
    my @GetOptConf=(
		    'all-five|a' => \$get_all_five,
		    'all-species|A' => \$get_all_species,
		    'flat-file|f=s'=>\$flat_file ,
		    'species|s=s' => \$species,
		    'missing-psi|M=s' => \$missing_psi_file , 
		    'data-dir|D=s' => \$data_dir , 
		    'tmp_dir|T=s' => \$tmp_dir , 
		    'binary|b' => \$binary , 
		    'blast|B' => \$dodablast,
		    'blast-file:s' => \$blast_file,
		    'force|F' => \$force , 
		    'help|h' => \$help,
		    'verbose|v' => \$verbose,
		    'spacer|S' => \$spacer,
		    'psi2id|i:s' => \$psi2id , 
		    'psi-file|p:s' => \$psi_file , 
		    'log-file|l' => \$LOG, 
		    'net-file|n:s' => \$net_file,
		    'debug|d' => \$debug,
                    'test|t'=>\$TEST
		   );
    GetOptions(@GetOptConf) or do {
	exit;
    };
    $help && do {usage(@GetOptConf)};

    ###########################################################
    # If no tmp dir was given by the user, use default.	      #
    # First check if it exists and if it does not, create it. #
    ###########################################################
    eval { mkpath($tmp_dir) };
    if ($@) {
	die "Couldn't create $tmp_dir: $@";
    }

    ####################################################
    # If no species name has been given, try and guess #
    # it from the flat file's name		       #
    ####################################################
    $flat_file && do {$species//=guess_species($flat_file)};
    
    #####################################
    # Translate species to default name #
    #####################################
    $species=get_species($species, "-s");

    #####################################################
    # Carp if both -a and -A options have been passed   #
    #####################################################
    if ($get_all_five && $get_all_species) {
	print STDERR "ERROR: The options -a and -A are mutually exclusive\n";
	exit(1);
    }
    #####################################################
    # If no flat file has been given, look for default. #
    #####################################################
    $flat_file//="$data_dir/$species.flat";

    #####################################################
    # If no blast file has been given, use for default. #
    #####################################################
    $blast_file//="$tmp_dir/$$.out";

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
	-z $psi_file && do {die("FILE: $psi_file exists but is empty! Please delete it and run again.")};
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
    check_file($missing_psi_file,"w") ;

    ####################################
    # If a map file has been passed as #
    # an option, use it.	       #
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
	    print STDERR fold_at_word("\n No psi2id file name was given (-i) and there is already a file with the default name ($psi2id), continuing will overwrite it.", 65,  " "), "\n" unless $force;
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
    my @options=($TEST,$flat_file,$species,$missing_psi_file,$binary,$help,$psi2id,$psi_file,$get_psi_file,$net_file,$generate_map,$generate_network,$LOG, $tmp_dir, $dspacer, $data_dir, $dodablast, $blast_file, $get_all_five, $get_all_species, $spacer);
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
	## Trick to allow me to change the species name for the
	## get_psi_interactome call only
	my $spc = $species;
	if ($get_all_five) {
	    $spc='five';
	}
	elsif ($get_all_species) {
	    $spc='all';
	}
        for my $i (0..$db_count-1) {
            my $ii=$i+1; ## get_psicquic_interactome.pl takes nums from 1 to $db_count-1
            my $tmp_psi_file=$tmp_dir . "/$$.$species.$dbs{$ii}.psi";
            my $tmp_error_file=$tmp_dir . "/$$.$species.$dbs{$ii}.psi.er";
            $files[$i]=$tmp_psi_file;
            $e_files[$i]=$tmp_error_file;
            my $cmd="get_psicquic_interactome.pl -d $ii -s $spc > $tmp_psi_file 2>$tmp_error_file && echo \"$dbs{$ii} Done\" >> $prog_file &";
            debug($cmd);
	    system($cmd) unless $TEST ;
        }
        my $finished=0;
        $finished=$db_count if $TEST;        
        my %hash;
        #######################################################
        # $finished will == $db_count when all		      #
        # get_psicquic_interactome.pl instances are finished. #
        #######################################################
	my $last_finished=0;
	while ($finished<$db_count) {
            my $done_files=$finished;
            -e $prog_file && do {
                $finished=`wc -l $prog_file`;
                $finished=~s/(\d+).+/$1/s;
            };
            sleep(2);
            if (($done_files!=$finished) && ($finished!=$last_finished))  {
                v_say("DOWNLOADED: $finished  of $db_count\r", $spacer,1);
		$last_finished=$finished;
            }
	}
    
        v_say("DOWNLOADED: $finished  of $db_count", $spacer);
        $files{NAMES}{'psi'}=$psi_out;
        $files{CREATED}{'psi'}="Created";
        $TEST ?
            v_say("TESTING: cat @files | grep -v \"^Total:\"> $psi_out") :
		
            system("cat @files | grep -v \"^Total:\"> $psi_out");
	###########################
        # Remove tmp files	  #
        ###########################
	map{unlink($_)}(@files,@e_files, $prog_file);

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
        $files{CREATED}{'psi2id'}='Existing';
    } 
    else {
	die("FILE: $psi_file is empty!") if -z $psi_file;
	my $cmd="make_psicquic_map.pl -pT $tmp_dir -";
	$verbose && do {$cmd.="v"};
	$debug && do {$cmd.="d"};
	$binary && do {$cmd.="b"};
	$missing_file//="$tmp_dir/$species.$$.Psi2Uni.missing";
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
        $files{CREATED}{'psi2id'}="Created";

    }
    ###############
    # Count stuff #
    ###############
    my ($inter_count, $id_count, $map_count, $mis_count)=(0,0,0,0);
    if (-s $psi_file) {
	v_say("Counting Psicquic Interactions...($psi_file) ", $spacer,1);
	my $fh=check_file($psi_file, "r");
	my %ll;
	my $taxid=get_species($species, "-s", 1);
	while (<$fh>) {
	    chomp;
	    next unless /taxid:$taxid/;
	    my @bb=split(/\t/);
	    my $pair=join("", sort($bb[0],$bb[1]));
	    $ll{PAIRS}{$pair}++; ## Count unique interactIONS
	    $ll{IDs}{$bb[1]}=$ll{IDs}{$bb[0]}=1; ## Count unique interactORS
	}
	$inter_count=scalar(keys(%{$ll{PAIRS}}));
	$id_count=scalar(keys(%{$ll{IDs}}));
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
	    next if /^\s*$/;
	    my @bb=split(/\t/);
	    push @{$kk{$bb[0]}{'UniID'}}, $bb[1];
            push @{$kk{$bb[1]}{'PsiID'}}, $bb[0];
	    $map_count++;
	}
	#$map_count=(scalar keys %kk)/2;
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
	p_say(" Building HQ network $psi_file", 60, $spacer);
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
    $inters=scalar keys %net;
    return($inters, $prots, \%net)
}
############################################################
############################################################
sub run_network_blast{
    my $net_file=shift;
    my $flat_file=shift;
    my $blast_file=shift||"$tmp_dir/$$.netnames.pep.cdh";
    
    my $fh=check_file($net_file,"r");
    open(my $fh2, ">","$tmp_dir/$$.netnames")||die("Could not open $tmp_dir/$$.netnames for writing : $!\n");
    debug("Writing netnames : $tmp_dir/$$.netnames\n");
    my %NAMES;
    while (<$fh>) {
	chomp;
	next if /^\s*$/; 
	my @a=split(/\t/);
	$NAMES{$a[0]}++;
	$NAMES{$a[1]}++;
	print $fh2 "$a[0]\n" unless $NAMES{$a[0]}>1; 
	print $fh2 "$a[1]\n" unless $NAMES{$a[1]}>1; 
    }
    close($fh2);
    v_say("Fetching FASTA seqs from UniProt", $spacer);
    ## Various problems when downloading from uniprot, parse flat instead
    #my $cmd="uniprot_fetch.pl $tmp_dir/$$.netnames > $tmp_dir/$$.netnames.pep";
    my $cmd="uniprot_parse.pl -Ff ID -i  $tmp_dir/$$.netnames $flat_file > $tmp_dir/$$.netnames.pep";
    debug("BLAST1: $cmd");
    `$cmd`;
    v_say("Blasting...", $spacer);
    my $nr_net;
    if ($net_file=~/^(.+)\.(.+?)$/) {
	$nr_net=$1 . ".nr." . $2;
    }
    else {
	$nr_net=$net_file .  ".nr";
    }
    $cmd="network_blast.pl -v";
    $cmd.="d" if $debug;
    $psi2id=~/(.+).psi2id/;
    my $blast_map=$1 . ".blastmap";
    $cmd.="f $flat_file -M $blast_map -b $blast_file -F $tmp_dir/$$.netnames.pep $net_file > $nr_net";
    debug("BLAST2: $cmd");
    `$cmd`;
    close($fh);
    $fh=check_file($nr_net,"r");
    my (%ll,%kk);
    while (<$fh>) {
	chomp;
	next if /^\s*$/; 
	my @a=split(/\t/);
	my $pair=join("\t", sort(@a));
	$kk{$a[0]}=$kk{$a[1]}=1;
	$ll{$pair}++;	
    }
    $files{NAMES}{nr}=$nr_net;
    $files{CREATED}{nr}="Created";
    ## how many unique proteins are there in the network?
    my $prots=scalar keys %kk;
    my $inters=scalar keys %ll;
    ## Cleanup
   # system("rm $tmp_dir/$$.netnames*");
    return($inters, $prots);
	
    
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
$spacer Interactions in NR network file      : $nr_NetInters 
$spacer Proteins in NR network file          : $nr_NetProts
$spacer------------------- FILES CREATED --------------------------
$spacer Psicquic Interactions file           : $files{NAMES}{psi} ($files{CREATED}{psi})
$spacer Psicquic ID to Uniprot ID map file   : $files{NAMES}{psi2id} ($files{CREATED}{psi2id})
$spacer Network file                         : $files{NAMES}{net} ($files{CREATED}{net})
$spacer NR Network file                      : $files{NAMES}{nr} ($files{CREATED}{nr})
EndOfMsg
    v_say($rep);
}

############################################################
############################################################
sub PrintOptions{
    p_say(" OPTIONS ", 60, $_[$#_]);
    my  @a=qw($TEST $flat_file $species $missing_psi_file $binary $help $psi2id $psi_file $get_psi_file $net_file $generate_map $generate_network $LOG  $tmp_dir  $dspacer $data_dir $spacer );

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
    my $desc="This script will parse a pscquic interactions file and create a high quality PPI network file. If the psicquic file name that is given is empty, the script will automatically download ALL psicquic interactions from the following online databases: APID, BioGrid, IntAct, InnateDB, DIP, MINT, MatrixDB, iRefIndex, BIND, STRING, Reactome, InnateDB-IMEx, MolCon, I2D-IMEx, I2D, Spike, TopFind\n";

    my %opts=(
	      "usage" => $us,
	      "desc" => $desc,
	      "all-five|a" => "Download interactions for human, mouse, yeast, worm and fly.",
	      "all-species|A" => "Download interactions for ALL species.",
	      "f|flat-file" => "Flat file to be used for mapping.",
	      "M|missing-psi" => "Missing psicquic file. Psicquic IDs from the input file that could not be linked to a UniProt ID will be printed to this file.",
	      "binary|b" => "Accept only BINARY interactions when building the network.",
	      "blast|B" => "Retrieve FASTA seqs for each of the sequences in the created network, blast against each other and identify those sequences with 100%id, where query!=sbjct and where the entire sequence of the subject is aligned perfectly against the query for its ENTIRE length. If such a sequence pair consists of a TrEMBL prot that is a subsequence of a SWISPROT one, remove the trEMBL ID from the network, replacing it with the ID of the  SWISPROT protein.",
	      "blast-file" => "Blast out file to parse. If none is given, the sequences are retrieved and blasted (see -b)",
	      "species|s" => "The species we are parsing. If none is given the value is guessed from the flat file",
	      "net-file|n" => "Network File. If this file exists, do not parse the psicquic file and create a network this network instead. If it doesn't, the psicquic file will be parsed and the resulting network printed into this file. If no file is given, the psicquic file will be parsed and the network will be printed into <SPECIES_NAME>.<NETWORK_TYPE>.gr.",
	      "psi2id|i" => "psi2id map file. If this file exists, do not parse/convert ids from the psicquic input file use this map instead. If it doesn't, the psicquic file will be parsed and the map printed into this file. If no file is given, the psicquic file will be parsed and the map will be printed into <SPECIES_NAME>.psi2id.",
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
