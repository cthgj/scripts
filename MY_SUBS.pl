use Carp;
use Data::Dumper;
use feature qw(switch say);
$Data::Dumper::Indent = 1;


#####################################################
# Translate GO terms GO descriptions and vice versa #
#####################################################
sub terms2gos{
    my (%terms_to_GOs,%GOs_to_terms, %ret);
    my $terms_file=shift||"$ENV{HOME}/research/testing/new/data/GO.terms_alt_ids";
    my @terms=@_;
    my %res;
    my $T=check_file($terms_file, 'r');
    while (<$T>) {
	next if /^\!/ || /^\s+$/; 
	chomp;
	s/^\s*//;
	my @a=split(/\t/);
	push my @gos, $a[0];
	## If we have alternative IDs
	if ($a[1]=~/./) {
	    map{push @gos,$_}split(/\s+/,$a[1]);
	}
	$a[2] =~s/\s/_/g;
	map{
	    $terms_to_GOs{$a[2]}=$_;
	    $GOs_to_terms{$_}=$a[2];
	}@gos;
    }
    foreach (@terms) {
	if (/^\s*GO:\d+\s*$/){ 
	    $res{$_}= $GOs_to_terms{$_};    
	}
	else {
	    $res{$_}= $terms_to_GOs{$_};
	}
    }
return(%res);
}


##############################################################
# Run uniprot_map.pl, splitting large queries into 	     #
# smaller ones. The first argument must be a hash of	     #
# hashes whose keys are types of ids and whose values	     #
# are ids: $to_map{RefSeq}{NP_12345}.			     #
# Returns two hashes:		 			     #
#     %{$mapped{$id}}, %{$not_mapped{$id}}	 	     #
##############################################################
sub map_from_UniProt{
    my $hash=shift();
    my %to_map=%{$hash};
    $hash=shift;
    my %mapped=%{$hash};
    my $to_id=shift()||"ID";
    my $species=shift()||"unknown_spc";
    my $tmp_dir=shift()||"/tmp";
    my $spacer=shift();
    $spacer or do {$spacer="\t"};
    my (%mapped_ids, %not_mapped); 
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
    from:foreach my $from (@froms){## For each type of ID
	next if $from eq 'REV';
	next if $from eq 'UNK';
	next if $from eq 'PDB';
	my @missing;
	my %k;
	my $counter=0;
	map{
	    ##################################################
            # Skip if this id has already been mapped	     #
            ##################################################
	    #next from if defined($map{$_});
	    unless (defined($mapped{$_})){
		$counter++;
		$k{$to_map{$from}{$_}}++; ## avoid duplicates
		push @missing, $to_map{$from}{$_} unless $k{$to_map{$from}{$_}}>1;
	    }	    
	}keys(%{$to_map{$from}});
	next unless $counter>0;
	my @files=("$tmp_dir/$species.names.$$.$from.0");
	open(my $fh,">", "$files[$#files]") or croak "cannot open $files[$#files] > : $!";
	$counter=0;
	my $cc=0;

	foreach my $missing_id (@missing){
	    #####################################################
            # Skip this id if it has already been mapped        #
            #####################################################
	    
	    unless (defined($mapped_ids{$to_map{REV}{$missing_id}})){
		print $fh "$missing_id\n" ;
		$counter++;
		#########################################################
		# If there are MANY ids, break the uniprot query        #
		# into smaller ones.				        #
		#########################################################
		if($counter>1 && $counter % 2000 == 0){
		    $cc++;
		    close $fh;
		    debug("Opening $missing_id : $counter $tmp_dir/$species.names.$$.$from.$cc");
		    push @files, "$tmp_dir/$species.names.$$.$from.$cc";
		    open($fh,">", "$files[$#files]") or croak "cannot open $files[$#files] > : $!";
		    $counter=0;
		}
	    }
	}
	close $fh;
	####################
	# Get from UniProt #
	####################
	my $tt=$#files+1;
	for my $num (1..$tt){
	    my $file=$files[$num-1];
	    my $cmd="uniprot_map.pl -l";
	    #$verbose && do {$cmd.="v"};
	    $cmd.="f $from -t $to_id $file";
#	    my $cc=$num-1;
	    v_say("$from : uniprot_map.pl -lf $from -t $to_id $file\n\t\t",$spacer, 1);
	    my $map=`$cmd`;
	    ########################################################
	    # Read the mapped names into a hash and convert        #
	    ########################################################
	    my %l;
	    map{
		my @a=split(/\t/);
		####################################################
		# UniProt sometimes returns the SAME line twice	   #
		####################################################
		push @{$l{$a[0]}},$a[1] if $a[0];
	    }split("\n",$map);
	    ##########################################################
	    # Now build a hash connecting the original (bad) id      #
	    # to the mapped one					     #
	    ##########################################################
	from:foreach my $from (@froms){## For each type of ID
		foreach my $orig_id (keys(%{$to_map{$from}})){   
		    if(defined($l{$to_map{$from}{$orig_id}})) {
			map{
			    $mapped_ids{$orig_id}{$_}++;
			}@{$l{$to_map{$from}{$orig_id}}};
		    }
		    elsif (defined($mapped{$orig_id})) {
			##Assume the 1st ID is the good one...
			my @a=keys(%{$mapped{$orig_id}});
			$mapped_ids{$orig_id}{$a[0]}++;
			#print STDERR "\$mapped{$orig_id} : $mapped{$orig_id}"; die();
		    }
		    else{$not_mapped{$orig_id}++}
		}
	    }
	} ## end for my $num (1..$tt)
    } ## end foreach my $from
    ##########################################################
    # Now build a hash connecting the original (bad) id      #
    # to the mapped one					     #
    ##########################################################
 from:foreach my $from (@froms){## For each type of ID
	foreach my $orig_id (keys(%{$to_map{$from}})){   
	    if(defined($l{$to_map{$from}{$orig_id}})) {
		map{
		    $mapped_ids{$orig_id}{$_}++;
		}@{$l{$to_map{$from}{$orig_id}}};
	    }
	    elsif (defined($mapped{$orig_id})) {
		##Assume the 1st ID is the good one...
		my @a=keys(%{$mapped{$orig_id}});
		$mapped_ids{$orig_id}{$a[0]}++;
	    }
	    else{$not_mapped{$orig_id}++}
	}
    }
    
    return(\%mapped_ids, \%not_mapped);
}

##############################################################
# Run uniprot_fix_ids.pl, which will request FASTA seqs from #
# UniProt and then parses the FASTA header to extract IDs.   #
# This is useful because UniProt's mapping service does not  #
# find obsolete IDs but UniProt's sequence retrieval service #
# does. Go figure...					     #
##############################################################
sub run_uniprot_fix{
    my %ids=%{shift()};
    my %mapped_ids=%{shift()};
    my %n_mapped;
    my %new_map;
    my $tmp_file=shift() ||"/tmp/$$";
    my $spacer=shift() || "\t\t";
    my %to_fix;
    ####################################################
    # This will collect the IDs that were successfully #
    # mapped to UniProt IDs			       #
    ####################################################
    foreach my $orig_id (keys(%ids) ) {
	if (defined($mapped_ids{$orig_id})){
	    #my @a=keys(%{$mapped_ids{$orig_id}});
	    map{
		$to_fix{$_}=$orig_id;
	    }keys(%{$mapped_ids{$orig_id}});
	}
	else {
	    $n_mapped{$orig_id}++;
	}
    }
    ##########################################################
    # If any of those not mapped is an (obsolete) UniProt ID #
    # try mapping it as well				     #
    ##########################################################
    foreach my $mis (keys(%n_mapped)) {
	if ($mis=~/^[\w\d]{6}$/ || $mis=~/_\w{5}$/) {
	    $to_fix{$mis}=$mis;
	} 
    }
    my @ff=keys(%to_fix);
    if ($ff[0]) {
        ###################################
	# Split request into smaller ones #
	###################################
	my $fh;
	open( $fh,">", "$tmp_file.0") or croak "cannot open tmp_file $tmp_file.0 : $!";
	$counter=0;
	my $cc=0;
	
	print $fh "$ff[0]\n";
	for my $ii (1..$#ff) {
	    print $fh "$ff[$ii]\n";
	    if ($ii % 500 == 0) {
		$counter++;
		close($fh);
		my $kk="$tmp_file.$counter";
		open($fh,">", "$kk") or croak "cannot open $kk : $!";
	    }
	}
	my $finished=0;
	my $prog_file="$tmp_file.progress";
	v_say("Running uniprot_fix_ids.pl $counter..." );
	for (my $i=0; $i<=$counter; $i++) {
	    my $cmd="uniprot_fix_ids.pl ";
	    #$verbose && do {$cmd.="-v "};
#	    $cmd.=" $tmp_file.$i | gawk 'BEGIN{OFS=\"\\t\"}{print \$1,\$3}' > $tmp_file.$i.out && echo \"$i Done\" >> $prog_file &";
	    $cmd.=" $tmp_file.$i  > $tmp_file.$i.out && echo \"$i Done\" >> $prog_file &";
	    debug("$cmd",$spacer);
	    system($cmd);
	}
	########################################################
        # $finished will == $counter when all		       #
	# instances of uniprot_fix_ids.pl are finished.	       #
        ########################################################
	my $last_finished=0;
	while ($finished<=$counter) {  
	    my $done_files=$finished;
            -e $prog_file && do {
                $finished=`wc -l $prog_file`;
                $finished=~s/(\d+).+/$1/s;
            };
            sleep(2);
            if (($done_files!=$finished) && ($finished!=$last_finished))  {
		my $aa=$counter+1;
                v_say("FIXED: $finished  of $aa\r", $spacer,1);
		$last_finished=$finished;
            }	    
	}
	my $map=`cat $tmp_file.*.out | gawk 'BEGIN{OFS=\"\\t\"}{print \$1,\$3}'`;
	map{
	    ## The output format of uniprot_fix_ids.pl is:
	    ## Requested_ID\tRetrieved_AC\tRetrieved_ID
	    my @a=split(/\t/);
	    ## $to_fix{$a[0]} is the psicquic ID of $a[0]
	    $new_map{$to_fix{$a[0]}}{$a[1]}++;
	}split("\n",$map);
	
	###########################
	# Check, just in case     #
	###########################
	foreach my $id (keys(%ids)) {
	    if ( not defined($new_map{$id})) {
		if (defined($mapped_ids{$id})) {
		    my @a=keys(%{$mapped_ids{$id}});
		    map{$new_map{$id}{$_}++}@a;
		}
		
	    # iff (not defined($new_map{$id})){
		# 	print STDERR "Aaaaaaaaaaaaaargh : $id has no new_map \n" ;
		# }
	    }
	}
    }

    #####################################################
    # Now remove any that were mapped in the last step  #
    # from %not_mapped					#
    #####################################################
    my %koko;
    map{$koko{$_}++ unless defined($new_map{$_});}keys(%ids);
    return(\%new_map, \%koko);
}
###################################################################
###################################################################

###############################################################
# Check if output files exist, prompt user for action.        #
# USAGE: check_file(FILE_NAME, MODE). mode can be r,w or a    #
# for read, write or append. It returns the opened filehandle #
###############################################################
sub check_file{
    my $f=shift();
    croak("No file name given for check_file subroutine : $!\n") unless $f;
    my $mode=shift()|| die("No mode given for sub check_file\n");
    if($mode eq 'r'){
        $mode="<";
        ##############################
        # Check for compressed files #
        ##############################
        if ($f=~/\.gz$/) {
            $f="zcat $f";
            $mode="-|";
        }
    }
    elsif($mode eq 'a'){$mode=">>";}
    elsif($mode eq 'w'){$mode=">";}
    

    else{die("Invalid mode for check_file ($mode) must be one of w,r,a\n");}
    my $fh;
    ######################
    # If the file exists #
    ######################
    if(-e $f){ 
	## $force is a global our
	if($force){
	    open($fh, "$mode", $f)||return(undef);
	    return(*$fh);
	}
	else{
	    if($mode eq ">"){
		print STDERR fold_at_word("\n There is already a file called $f, continuing will overwrite it.", 65,  " "), "\n";
		while(0==0){
		    print STDERR fold_at_word("\n Type 'y' to overwrite the file, 'r' to rename it and 'n' to exit: [Y/n/r]: ", 65,  "  ");
		    my $resp=<STDIN>;
		    ############################
		    # If the user says yes     #
		    ############################
		    if($resp eq "\n" || $resp =~ /^\s*[yY]\s*$/){
			open($fh, "$mode", $f)||return(undef);
			return(*$fh);
		    }
		    ##########################
		    # If the user says no    #
		    ##########################
		    elsif($resp =~ /^\s*[nN]\s*$/){exit;}
		    ###################################
		    # If the user renames the file    #
		    ###################################
		    elsif($resp =~ /^\s*[rR]\s*$/){
			print STDERR "New filename:";
			$resp=<STDIN>;
			if($resp=~/[\w\d]/){
			    chomp($resp);
			    $fh=check_file($resp, "$mode");
			    return(*$fh);
			}
			########################################
			# If no response was given, try again  #
			########################################
			else{
			    check_file($f, "$mode");
			}
		    }
		}
	    }
	    else{      
		open($fh, "$mode", $f)||return(undef);
		return(*$fh);	
	    }
	    return(*$fh);
	}
    }
    ##############################
    # If the file does not exist #
    ##############################
    else{
	if ($mode eq '<') {
	    croak("Could not open file $f for reading: $!");
	    
	}
	open($fh, "$mode", $f)||return(undef);
	return(*$fh);
    }
}

###############################################
# Will print the 1st ARG passed to STDERR iff #
# global $debug is set			      #
###############################################
sub debug{ 
    my $msg=shift;
    chomp($msg);
    #$msg.="\n" unless $msg=~/\n$/;
    print STDERR "\tDEBUG: $msg \n" if $debug; 
}
#################################################
# p_say (p for pretty) will print the message   #
# to STDERR surrounded by a row of dashes added #
# until the length of $msg reaches $lim.        #
#################################################
sub p_say{
    my $msg=shift;
    my $lim=shift;
    my $spacer=shift;
    ($lim-length($msg)) % 2 != 0 && do {
	$msg.=" ";
     };
    my $ll=(($lim-length($msg))/2);
    my $kk="-" x $ll . $msg . "-" x $ll;
    v_say($kk, $spacer,);
}
###############################################################
# v_say can take up to 3 arguments. The first arg will be     #
# printed to STDERR with a \n appended unless one is already  #
# there. Any STRING passed will be prepended to the output.   #
# Any NUMBER passed will cause v_say to NOT append a \n to    #
# its output.						      #
###############################################################
sub v_say{
    return() unless $verbose;
    my $msg=shift();
    my ($spacer, $nl);
    ##################################
    # If there are 2 more arguments. #
    ##################################
    if($#_>0){
	$nl=1;
	# If $_[0] is a num and $_[1] is a string
	if (!($_[0] & ~$_[0]) && ($_[1] & ~$_[1])){
	    $spacer=$_[1];
	}
	# If $_[1] is a num and $_[0] is a string
	elsif(!($_[1] & ~$_[1]) && ($_[0] & ~$_[0])){
	    $spacer=$_[0];
	}
	else{croak "If 3 args are given to v_say, two must be strings and one must be a number.\n"  }
    }
    ########################################################
    # If there is only 1 more argument, if the arg is 	   #
    # a string it will be the spacer. If it is >=1 it will #
    # cause v_say to NOT add a \n at the end of $msg.	   #
    ########################################################
    else{
	####################################################
        # If there are no more arguments, 	 	   #
	# print and exit.				   #
        ####################################################
	unless ($_[0]){
	    $msg.="\n" unless $msg=~/\n$/;
	    print STDERR "$msg" ;
	    return(); 
	}
	## If $_[0] is a number
	if (!($_[0] & ~$_[0])){
	    $nl=1;
	}	
	else{$spacer=$_[0];}
    }
    !$nl && do{$msg.="\n" unless $msg=~/\n$/;};
    $msg=~s/^\t//;
    $spacer||="" ;
    print STDERR $spacer, "$msg"; 
    return(); 
}
##########################################################
# Will print the 1st ARG passed to STDERR piped through  #
# more and then exit. It will also trim $0 to get the 	 #
# program name.	If no 2nd arg is passed, it will pipe    #
# through more, otherwise through less.                  #
##########################################################
sub print_help_exit{
    my %msg=%{shift()};
    my $ll=shift;
    ## If a 2nd arg is given, we do want to pipe through less.
    my $H;
    if ($ll) {
	open ($H, "| less") or croak "MY_SUBS.pl, print_help_exit: Could not pipe to less: $!\n";
    }
    else {
	open ($H, ">&", \*STDERR) or croak "MY_SUBS.pl, print_help_exit: Could not open STDERR: $!\n";
    }

    $0=~/.+?([^\/]+)$/;
    my $p_name=$1;
    my $INDENT="  ";
    $msg{'desc'}=~s/$0/$p_name/;
    print $H fold_at_word("\nDESCRIPTION:\n$INDENT$msg{desc}\n", 80, "$INDENT");
    my @opts=sort { lc $a cmp lc $b } grep(!/desc|usage/, keys(%msg));
    print $H <<EndOfHelp;

USAGE: 
\t$p_name $msg{'usage'}

OPTIONS:
EndOfHelp
    my %to_print;
## Parse option/description pairs
    foreach my $opt (@opts){
	## Are there synonyms?
	if($opt =~/\|/){
	    my @oo=sort {length($a) >= length($b)} split(/\s*\|\s*/, $opt);
	    my $first_opt=shift(@oo);
	    ## 1st option
	    length($first_opt) > 1 ? 
		($first_opt="--"  . $first_opt ):
		($first_opt="-"  . $first_opt);
	    ## Synonyms
	    map{
		length() > 1 ? 
		    ($_="--"  . $_ ):
		    ($_="-"  . $_ );
	    }@oo;
	    my $new_opt=$first_opt;
	    map{$new_opt.= ", $_"}@oo;
	    $to_print{$opt}=$new_opt;
	}
	else{
	    length($opt) > 1 ? 
		($new_opt="--"  . $opt ):
		($new_opt="-"  . $opt );
	    $to_print{$opt}=$new_opt;
	}

	## Parse the descriptions
	$desc=$msg{$opt};
	$desc=~s/$0/$p_name/g;
	my $SPACER=" ";
	while(length($to_print{$opt})<10){
	    $to_print{$opt}.=" ";
	}
	if(length($desc)>80){
            ## Sometimes I format it by hand, so do not
            ## fold if the line contains a newline
            unless ($desc =~/\n/) {
                $desc=fold_at_word($desc, 50, " " x 26 );
            }
            
	}
	print $H "  $to_print{$opt} $SPACER $desc\n";
    }
    exit;

}

##############################################
# This will return a species name if known.  #      
##############################################
sub get_species{
    my $species=shift();
    my $spc_opt=shift()|| "-s";
    my $return_num=shift();
    ## Check options
    if ($spc_opt =~/^\d*$/){
	$return_num=$spc_opt;
	$spc_opt="";
    }
    if   (($species eq "human") || 
	  ($species eq "9606") || 
	  ($species eq "hum") || 
	  ($species eq "hs")){
	$species='human';
	$num=9606;
    }
    elsif(($species eq "fly")   || 
	  ($species eq "7227")  || 
	  ($species eq "dmel")  || 
	  ($species eq "fb")  || 
	  ($species eq "dm") || 
	  ($species eq "dro")){
	$species='fly';
 	$num=7227;
   }
    elsif(($species eq "worm")  || 
	  ($species eq "6239")  || 
	  ($species eq "wb")  || 
	  ($species eq "ce") || 
	  ($species eq "ele")){
	$species='worm';
 	$num=6239;
   }
    elsif(($species eq "mouse") || 
	  ($species eq "mg")  || 
	  ($species eq "10090")  || 
	  ($species eq "mm") || 
	  ($species eq "mus")){
	$num=10090;
	$species='mouse';
    }
    elsif(($species eq "yeast") || 
	  ($species eq "4932") || 
	  ($species eq "sgd") || 
	  ($species eq "scc")){
	$num=4932;
	$species='yeast';
    }
    elsif ($species eq "6334") {
	$num=6334;
	$species="trichinella";
    }
    elsif($species eq "unknown_spc"){
	$num=-1;
    }
    else{

	
print STDERR   <<EndOfMsg;
Species ($spc_opt) must be "unknown_spc" or one of the following:
         "human" || "hum" || "hs" || "9606"
	 "fly"   || "fb"  || "dm" || "dro" || "dmel" || "7227"
	 "worm"  || "wb"  || "ce" || "ele" || "6239"
	 "mouse" || "mg"  || "mm" || "mus" || "10090"
	 "yeast" || "sgd" || "scc"|| "4932"
EndOfMsg
exit(0);
    }
    $return_num ? return($num) : return($species);
}
##################################################################
# This sub will return the UniProt Suffix for the given species. #
##################################################################
sub get_suffix {
    my $species=shift;
    $species=get_species($species);
    my %SUF=(
	     'human' => 'HUMAN' ,
	     'fly'   => 'DROME' ,
	     'worm'  => 'CAEEL' ,
	     'mouse' => 'MOUSE' ,
	     'yeast' => 'YEAST' ,
	    );
    $SUF{$species} && do{return($SUF{$species})}

}
####################################################
# This will try and guess a species from a string. #
####################################################
sub guess_species{
    ###########################################################
    # The string (e.g. a filename) to search for a species in #
    ###########################################################
    my $fname=shift(); 
    if (($fname=~/\bhuman\b/) || 
	 ($fname=~/\bhum\b/) || 
	  ($fname=~/\bhs\b/)){
	      $species='human';
	  }
    elsif(($fname=~/\bfly\b/)   || 
	  ($fname=~/\bfb\b/)  || 
	  ($fname=~/\bdm\b/) || 
	  ($fname=~/\bdro\b/)){
	$species='fly';
    }
    elsif(($fname=~/\bworm\b/)  || 
	  ($fname=~/\bwb\b/)  || 
	  ($fname=~/\bce\b/) || 
	  ($fname=~/\bele\b/)){
	$species='worm';
    }
    elsif(($fname=~/\bmouse\b/) || 
	  ($fname=~/\bmg\b/)  || 
	  ($fname=~/\bmm\b/) || 
	  ($fname=~/\bmus\b/)){
	$species='mouse';
    }
    elsif(($fname=~/\byeast\b/) || 
	  ($fname=~/\bsgd\b/) || 
	  ($fname=~/\bscc\b/)){
	$species='yeast';
    }
    elsif($fname=~/([^\/]+).flata/){
	$species=$1;
    }
    elsif($fname=~/([^\/]+).psi/){
	$species=$1;
    }
    else{$species='unknown_spc';}
    ## If we want the taxid
    if ($_[0]) {
	$species=get_species($species,1);
    }
    return($species)
}

######################################################
# Fold a long string at specified length, respecting #
# word boundaries				     #
######################################################
sub fold_at_word{
    my $s=shift()|| croak "&fold_at_word needs 2 arguments";
    my $l=shift()|| croak "&fold_at_word needs 2 arguments";
    my $spacer=shift()|| "";
    return($s) if length($s)<=$l;
    my $ret="";
    my @a=split(/ /, $s);
    my $L=0;
    for my $c (0..$#a-1){
	$ret.="$a[$c] ";
	$L+=length($a[$c]);
	$L > $l && do {
	    $ret.="\n$spacer";
	    $L=0;
	}
    }
    $L+=length($a[$#a]);
    $L > $l && do {
	$ret.="\n$spacer";
	$L=0;
    };
    $ret.="$a[$#a] ";
    return($ret);
}


##############################################################
# Get a list of accepted MIs (interaction detection methods) #
##############################################################
sub get_hq_MIs{
    my $mode=shift||'b';
    my (%ehq,%hq);
    ## Extra MIS wanted by Christine
    if ($mode eq 'e') {
	$mode='b'; ## get the rest
	%ehq=(
	      "MI:0019" => "coimmunoprecipitation",
	      "MI:0006" => "anti bait coimmunoprecipitation",
	      "MI:0007" => "anti tag coimmunoprecipitation",
	      "MI:0858" => "immunodepleted coimmunoprecipitation",
	      "MI:0096" => "pull down",
	      "MI:0963" => "interactome parallel affinity capture",
	      "MI:0676" => "tandem affinity purification"
	     );
    }
    if($mode eq 'b'){
	%hq=(
	    "MI:0008" => "array technology",
	    "MI:0009" => "bacterial display",
	    "MI:0010" => "beta galactosidase complementation",
	    "MI:0011" => "beta lactamase complementation",
	    "MI:0012" => "bioluminescence resonance energy transfer",
	    "MI:0013" => "biophysical",
	    "MI:0014" => "adenylate cyclase complementation",
	    "MI:0016" => "circular dichroism",
	    "MI:0017" => "classical fluorescence spectroscopy",
	    "MI:0018" => "two hybrid",
	    "MI:0020" => "transmission electron microscopy",
	    "MI:0030" => "cross-linking study",
	    "MI:0031" => "protein cross-linking with a bifunctional reagent",
	    "MI:0034" => "display technology",
	    "MI:0040" => "electron microscopy",
	    "MI:0041" => "electron nuclear double resonance",
	    "MI:0042" => "electron paramagnetic resonance",
	    "MI:0043" => "electron resonance",
	    "MI:0047" => "far western blotting",
	    "MI:0048" => "filamentous phage display",
	    "MI:0049" => "filter binding",
	    "MI:0051" => "fluorescence technology",
	    "MI:0052" => "fluorescence correlation spectroscopy",
	    "MI:0053" => "fluorescence polarization spectroscopy",
# "MI:0054" => "fluorescence-activated cell sorting",
	    "MI:0055" => "fluorescent resonance energy transfer",
	    "MI:0065" => "isothermal titration calorimetry",
	    "MI:0066" => "lambda phage display",
	    "MI:0073" => "mrna display",
	    "MI:0081" => "peptide array",
	    "MI:0084" => "phage display",
	    "MI:0089" => "protein array",
	    "MI:0090" => "protein complementation assay",
	    "MI:0091" => "chromatography technology",
	    "MI:0092" => "protein in situ array",
	    "MI:0095" => "proteinchip(r) on a surface-enhanced laser desorption/ionization",
	    "MI:0097" => "reverse ras recruitment system",
	    "MI:0098" => "ribosome display",
	    "MI:0099" => "scintillation proximity assay",
	    "MI:0107" => "surface plasmon resonance",
	    "MI:0108" => "t7 phage display",
	    "MI:0111" => "dihydrofolate reductase reconstruction",
	    "MI:0112" => "ubiquitin reconstruction",
	    "MI:0114" => "x-ray crystallography",
	    "MI:0115" => "yeast display",
	    "MI:0226" => "ion exchange chromatography",
	    "MI:0227" => "reverse phase chromatography",
	    "MI:0231" => "mammalian protein protein interaction trap",
	    "MI:0232" => "transcriptional complementation assay",
# "MI:0254" => "genetic interference",
	    "MI:0255" => "post transcriptional interference",
# "MI:0256" => "rna interference",
# "MI:0257" => "antisense rna",
	    "MI:0369" => "lex-a dimerization assay",
	    "MI:0370" => "tox-r dimerization assay",
	    "MI:0397" => "two hybrid array",
	    "MI:0398" => "two hybrid pooling approach",
	    "MI:0399" => "two hybrid fragment pooling approach",
	    "MI:0400" => "affinity technology",
	    "MI:0401" => "biochemical",
	    "MI:0405" => "competition binding",
	    "MI:0406" => "deacetylase assay",
	    "MI:0410" => "electron tomography",
	    "MI:0411" => "enzyme linked immunosorbent assay",
	    "MI:0415" => "enzymatic study",
	    "MI:0416" => "fluorescence microscopy",
	    "MI:0419" => "gtpase assay",
	    "MI:0420" => "kinase homogeneous time resolved fluorescence",
	    "MI:0423" => "in-gel kinase assay",
	    "MI:0424" => "protein kinase assay",
	    "MI:0425" => "kinase scintillation proximity assay",
	    "MI:0426" => "light microscopy",
	    "MI:0428" => "imaging technique",
# "MI:0430" => "nucleic acid uv cross-linking assay",
	    "MI:0432" => "one hybrid",
	    "MI:0434" => "phosphatase assay",
	    "MI:0435" => "protease assay",
	    "MI:0437" => "protein three hybrid",
# "MI:0438" => "rna three hybrid",
#"MI:0439" => "random spore analysis",
	    "MI:0440" => "saturation binding",
#"MI:0441" => "synthetic genetic analysis",
	    "MI:0508" => "deacetylase radiometric assay",
	    "MI:0509" => "phosphatase homogeneous time resolved fluorescence",
	    "MI:0510" => "homogeneous time resolved fluorescence",
	    "MI:0511" => "protease homogeneous time resolved fluorescence",
	    "MI:0512" => "zymography",
	    "MI:0513" => "collagen film assay",
	    "MI:0514" => "in gel phosphatase assay",
	    "MI:0515" => "methyltransferase assay",
	    "MI:0516" => "methyltransferase radiometric assay",
#"MI:0588" => "three hybrid",
	    "MI:0655" => "lambda repressor two hybrid",
	    "MI:0657" => "systematic evolution of ligands by exponential enrichment",
#"MI:0663" => "confocal microscopy",
	    "MI:0678" => "antibody array",
	    "MI:0695" => "sandwich immunoassay",
	    "MI:0696" => "polymerase assay",
# "MI:0697" => "dna directed dna polymerase assay",
# "MI:0698" => "dna directed rna polymerase assay",
# "MI:0699" => "rna directed dna polymerase assay",
# "MI:0700" => "rna directed rna polymerase assay",
	    "MI:0726" => "reverse two hybrid",
	    "MI:0727" => "lexa b52 complementation",
	    "MI:0728" => "gal4 vp16 complementation",
	    "MI:0809" => "bimolecular fluorescence complementation",
	    "MI:0813" => "proximity enzyme linked immunosorbent assay",
	    "MI:0824" => "x-ray powder diffraction",
	    "MI:0825" => "x-ray fiber diffraction",
	    "MI:0827" => "x-ray tomography",
	    "MI:0841" => "phosphotransferase assay",
	    "MI:0870" => "demethylase assay",
	    "MI:0872" => "atomic force microscopy",
	    "MI:0879" => "nucleoside triphosphatase assay",
	    "MI:0880" => "atpase assay",
	    "MI:0887" => "histone acetylase assay",
	    "MI:0889" => "acetylase assay",
	    "MI:0892" => "solid phase assay",
	    "MI:0894" => "electron diffraction",
	    "MI:0895" => "protein kinase A complementation",
	    "MI:0899" => "p3 filamentous phage display",
	    "MI:0900" => "p8 filamentous phage display",
	    "MI:0905" => "amplified luminescent proximity homogeneous assay",
	    "MI:0916" => "lexa vp16 complementation",
	    "MI:0920" => "ribonuclease assay",
	    "MI:0921" => "surface plasmon resonance array",
	    "MI:0946" => "ping",
	    "MI:0947" => "bead aggregation assay",
	    "MI:0949" => "gdp/gtp exchange assay",
	    "MI:0953" => "polymerization",
	    "MI:0968" => "biosensor",
	    "MI:0969" => "bio-layer interferometry",
	    "MI:0972" => "phosphopantetheinylase assay",
	    "MI:0976" => "total internal reflection fluorescence spectroscopy",
	    "MI:0979" => "oxidoreductase assay",
	    "MI:0984" => "deaminase assay",
	    "MI:0989" => "amidase assay",
	    "MI:0990" => "cleavage assay",
	    "MI:0991" => "lipid cleavage assay",
	    "MI:0992" => "defarnesylase assay",
	    "MI:0993" => "degeranylase assay",
	    "MI:0994" => "demyristoylase assay",
	    "MI:0995" => "depalmitoylase assay",
	    "MI:0996" => "deformylase assay",
	    "MI:0997" => "ubiquitinase assay",
	    "MI:0998" => "deubiquitinase assay",
	    "MI:0999" => "formylase assay",
	    "MI:1000" => "hydroxylase assay",
	    "MI:1001" => "lipidase assay",
	    "MI:1002" => "myristoylase assay",
	    "MI:1003" => "geranylgeranylase assay",
	    "MI:1004" => "palmitoylase assay",
	    "MI:1005" => "adp ribosylase assay",
	    "MI:1006" => "deglycosylase assay",
	    "MI:1007" => "glycosylase assay",
	    "MI:1008" => "sumoylase assay",
	    "MI:1009" => "desumoylase assay",
	    "MI:1010" => "neddylase assay",
	    "MI:1011" => "deneddylase assay",
	    "MI:1016" => "fluorescence recovery after photobleaching",
	    "MI:1019" => "protein phosphatase assay",
	    "MI:1024" => "scanning electron microscopy",
	    "MI:1026" => "diphtamidase assay",
	    "MI:1030" => "excimer fluorescence",
	    "MI:1031" => "protein folding/unfolding",
#"MI:1034" => "nuclease assay",
#"MI:1035" => "deoxyribonuclease assay",
	    "MI:1036" => "nucleotide exchange assay",
	    "MI:1037" => "Split renilla luciferase complementation",
	    "MI:1038" => "silicon nanowire field-effect transistor",
	    "MI:1087" => "monoclonal antibody blockade",
	    "MI:1088" => "phenotype-based detection assay",
	    "MI:1089" => "nuclear translocation assay",
	    "MI:1111" => "two hybrid bait or prey pooling approach",
	    "MI:1112" => "two hybrid prey pooling approach",
	    "MI:1113" => "two hybrid bait and prey pooling approach",
	    "MI:1137" => "carboxylation assay",
	    "MI:1138" => "decarboxylation assay",
	    "MI:1142" => "aminoacylation assay",
	    "MI:1145" => "phospholipase assay",
	    "MI:1147" => "ampylation assay"
	    );
	## Add the rest from mode e. Keys will
	## be empty unless the mode was originally 'e'
	map{
	    $hq{$_}=$ehq{$_};
	}keys(%ehq);
    }

    return(\%hq);
}

sub print_hash{
    my %h=%{$_[0]};
    if($#_>0){
	map{print STDERR "$_\t$h{$_}\n"}keys(%h);
    }
    else{
	map{print STDOUT "$_\t$h{$_}\n"}keys(%h);
    }
}

sub this_script_name{
    $0=~/.+?([^\/]+)$/;
    return($1);
}



1;
