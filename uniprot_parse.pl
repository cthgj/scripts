#!/usr/bin/perl
use strict;
use warnings;
use Data::Dumper;
use Getopt::Std;
use SWISS::Entry;
# use SWISS::GeneGroup;
# use SWISS::KW;
use Switch;

my %opts;
getopts('gARwsavhIrFef:n:t:i:m:M:S:',\%opts) || do { print "Invalid option\n"; exit(1); };
&usage() if $opts{h};
&print_ID_types if $opts{I};
&usage() unless $ARGV[0];
my $toID=$opts{t}||"AC";
my $fromID=$opts{f}|| 'unknown';
my $synonyms=$opts{A}||undef;
my $flat=$ARGV[0]||die("Need a flat file\n");
my $missing_file=$opts{m}||undef;
my $not_mapped_file=$opts{M}||undef;
my $names = $opts{n}||undef;
my $IDlist;
my $verbose=$opts{v}||undef;
my $all_names=$opts{a}||undef;
my %entries;
my $return_fasta=$opts{F}||undef;
my $looking_for;
my $allow_fragments=$opts{r}||undef;
my $mark_fragments=$opts{R}||undef;
my $return_entry=$opts{e}|| undef;
my $spacer=$opts{S}||"";
$opts{s}++ if $opts{M};
my $strict_name=$opts{g}||undef;
my %ids=%{&read_ids()};

#print "\n"; ##for the webserver. If no results are found then the outfile is empty and php can't open it...

my %pdbs;

###################################################
# Marking fragments only makes sense when mapping #
###################################################
if ($mark_fragments && ($return_entry || $return_fasta)) {
    print STDERR "\tERROR : The mark fragments option (-R) cannot be used with the -F or -e flags\n";
    exit;
}
my @map;

##############
## Change the line termination string so we read an 
## entire entry at a time
##############
local $/ = "\n//\n";

if($flat=~/\.gz$/){
    open(A,"zcat $flat |")|| die("cannot open flat file $flat : $! \n");
}
else{
    open(A,"$flat")|| die("cannot open flat file $flat : $!: $@\n");
}
#########################################################
# Gene names, especially CG12345-style FlyBase names,   #
# can sometimes be found as other terms. If found, they #
# will be stored in $gene_name.			        #
#########################################################
my $gene_name;

my $species;
while (<A>) {
    $gene_name=undef;
    &print_report() if $looking_for==0;
    ##############
    ## Read in all the entries and fill %entries 
    ##############
    my $entry = SWISS::Entry->fromText($_);
    #  print "AA : " . $entry->DEs->elements()->category . "\n";
    
   
    # my $kw = new SWISS::KW;
    # $kw->text('Gelelectrophoresis');
    # $entry->KWs->add($kw);
    # foreach my $kw ($entry->KWs->elements) {
    # print $kw->text, ", ";
    # }die();
    ###################
    # Get the species #
    ###################
    $.==1 && do {
	my @A1=$entry->OXs->NCBI_TaxID()->elements();
	$species=$A1[0]->text;
	$species == 7227 &&  do {
	    ## Warn about CG preference
	    if ($toID eq 'GeneName' && $verbose){
	     print STDERR $spacer, "\n\tWARNING: For fly, precedence is given to CG1234-style names unless a Gene Name is specifically given in the flat file\n" ;
	 };
	};
    };
    print STDERR $spacer, "$looking_for/$.\r" if $verbose;
    
    #####################################
    # If this is a fly flat file, check #
    # for CG style names	        #
    #####################################
    if ($species == 7227 && /(CG\d{5,8})/){
	$gene_name||=$1 ;
    }; 


    my @want=();
    ##############
    ## Read entry's ID
    ##############
    my $id=${$entry->IDs->list}[0];

    ##############
    ## If we know what we are looking for
    ## or if we are looking for an ID
    ##############
    if($fromID eq 'ID'){
	next unless defined($ids{WANT}{$id});
    }
    
    my $frag=0;
    if(defined($ids{WANT}{$id})){
	###############################
	# Deal with protein fragments #
	###############################
	if (defined($entry->DEs->hasFragment)) {
	    $frag=1;
	    $ids{FOUND}{$id}="Fragment";

	}

	$return_entry && do {
	    $ids{FOUND}{$id}++;
	    $ids{MAPPED}{$id}++;
	    print STDOUT $entry->toText;
	    $looking_for--;
	    next;
	};
	$return_fasta && do {
	    &return_fasta($id,$entry,'bob');
	    $looking_for--;
	    next;
	};
	##############
	## Return whatever it is we 
	## are looking for.
	##############
	print_output($entry,$toID,$id,$id,$frag);

	# my @kk=($entry,$toID,$id,'bob');
	# push @map, \@kk; 

    }
    elsif($fromID ne 'unknown'){
	
	my %kk=%{&check_entry($entry,$fromID, $gene_name)};
	#@want=&check_entry($entry,$fromID);
	@want=keys(%kk);
	if($want[0] eq 'bob'){
	    next;
	}
	else{
	    next if $fromID eq 'PDB';
	    ###############################
	    # Deal with protein fragments #
	    ###############################
	    if (defined($entry->DEs->hasFragment)) {
		if ($mark_fragments) {
		    print "$id\tFragment\n";
		    $looking_for-- unless defined($ids{FOUND}{$id});
		}
		else {next unless $allow_fragments;}	
	    }
	    map{print_output($entry,$toID,$_,$id,$frag);}@want;
	    # map{
	    # 	my @kk=($entry,$toID,$_,$id);
	    # 	push @map, \@kk; 
	    # }@want;
	}
    }
    else{
	##############
	## Check whether we want
	## this entry and, if so,  print
	##############
	
	my %kk=%{&check_entry($entry,$fromID, $gene_name)};
	@want=keys(%kk);
	###############################
	# Deal with protein fragments #
	###############################
	if (defined($entry->DEs->hasFragment)) {
	    if ($mark_fragments) {
		print "$id\tFragment\n";
		$looking_for-- unless defined($ids{FOUND}{$id});
	    }
	    else {next unless $allow_fragments;}	
	}
	map{print_output($entry,$toID,$_,$id,$frag)}@want unless $want[0] eq 'bob';
	# map{
	#     my @kk=($entry,$toID,$_,$id);
	#     push @map, \@kk; 
	# }@want unless $want[0] eq 'bob';
	
    }
}

###########################################################
# PDBs often map to multiple, unrelated proteins. 	  #
# Therefore, only return a PDB match if it is the	  #
# ONLY match.					          #
###########################################################
if ($fromID eq 'PDB') {
    my @multi;
    foreach my $pdb (keys %pdbs) {
	if (scalar(@{$pdbs{$pdb}})==1) {
	    my $entry=${$pdbs{$pdb}}[0];
	    if (defined($entry->DEs->hasFragment)) {
		if ($mark_fragments) {
		    print "$pdb\tFragment\n";
		}
		else {next unless $allow_fragments;}	
	    }
	    print_output($entry,$toID,$pdb,'bob');
	}
	else {
	    push @multi, $pdb;
	}	
    }
    print STDERR scalar @multi  . " PDBs mapped to multiple proteins: @multi\n";
}
else {
    &print_report();    
}

############################################################
############################################################

sub print_report{
## Clean up
my (@not_mapped,@missing,@frags);
foreach my $id (keys(%{$ids{WANT}})){
    if ($id=~/FOUND/ || $id=~/MAPPED/){next}
    if(defined($ids{FOUND}{$id})){
	if ($ids{FOUND}{$id} eq "Fragment") {
	    push @frags, $id;
	}
	else {
	    push @not_mapped, $id unless defined($ids{MAPPED}{$id});	    
	}
    }
    else{
	push @missing, $id;	
    }

}
unless(scalar(@frags)==0){
    scalar(@frags) > 1 ?
	print STDERR "$spacer       " . scalar(@frags) . " IDs map to fragments " :
	print STDERR "$spacer       " . scalar(@frags) . " ID maps to a fragment " ;
    print STDERR "\t: @frags\n" unless($opts{s});
}

my $mm=scalar(@not_mapped);
unless($mm==0){
    $toID='UniProt AC' if $toID eq 'AC';
    $toID='UniProt' if $toID eq 'ID';
    $mm > 1 ?
	print STDERR "$spacer       $mm IDs had no $toID ID : " :
	print STDERR "$spacer       $mm ID had no $toID ID : ";

    ## Print to file if -M
    if($not_mapped_file){
	print STDERR "\n";
	open(my $fh , ">", $not_mapped_file) or die "Could not open not mapped file: $not_mapped_file: $!\n";
	map{
	    print $fh "$_\n";
	}@not_mapped;
	close($fh);
    }
    ## Print to STDERR unless -M
    else{
	unless($opts{s}){
	    print STDERR " : \t";
	    map{
		print STDERR  "$_ ";
	    }@not_mapped;
	}
	print STDERR  "\n";
    }
}
my $c=0;
$mm=scalar(@missing);
unless($mm==0){
    if ($verbose) {
	$mm > 1 ?
	    print STDERR $spacer, "       $mm of " . scalar keys(%{$ids{WANT}}) . " IDs were not found":
		print STDERR $spacer, "       $mm of " . scalar keys(%{$ids{WANT}}) . " ID was not found";
    }

    if($missing_file){
	open(MM, ">$missing_file") or die "\nCould not open $missing_file for writing: $!\n";
	map{
	    print MM "$_\n";
	}@missing;
    }
    else{
	unless($opts{s}){
	    if ($verbose) {   
		print STDERR" :\t";
		map{
		    print STDERR "$_ ";
		}@missing;
	    }
	}
    }
    print STDERR "\n";
}

exit();
}

############################################################
############################################################
sub print_output{
    my $entry=shift;  ## SWISS::Entry object
    my $toID=shift;   ## This is the ID type I'm requesting
    my $id=shift;     ## This is the ID I will return
    my $real_id=shift;## This is the ID I asked for
    my $is_frag=shift;## is 1 if this is a fragment
    ## Some IDs (eg flybase) can correspond to many
    ## entries. see FBgn0261565
    $looking_for-- unless defined($ids{FOUND}{$id});
    #######################
    # Is this a fragment? #
    #######################
    my $frag="";
    $is_frag!=0 && do {$frag="\tFragment"; };
    ########################################################
    # If this ID has already been found, do not print this #
    # 2nd (or Nth) id if it is a trEMBL id		   #
    ########################################################
    if (defined($ids{FOUND}{$id}) && ($ids{FOUND}{$id} ne "Fragment") ) {
	return if $entry->IDs->dataClass eq 'Unreviewed';
    }

    $ids{FOUND}{$id}++;
    $return_entry && do {
	$ids{MAPPED}{$id}++;
	print STDOUT $entry->toText;
	return;
    };
    $return_fasta && do {
	&return_fasta($id,$entry,$real_id);
	return;
    };
 
    if($toID eq 'AC'){	 
	$ids{MAPPED}{$id}++;
	my @names=@{$entry->ACs->list};
	$all_names ?
	    print "$id\t@names$frag\n" :
	    print "$id\t$names[0]$frag\n";
	return;  
	
    }
    ##############
    ## Do we want a gene name?
    ##############
    elsif($toID eq 'GeneName'){
	if(defined($entry->GNs->text)){
	    $ids{MAPPED}{$id}++;
	    #########################################################
            # Sometimes two gene names have no space between	    #
            # them. Split'em.					    #
            #########################################################
	    my $aa=$entry->GNs->text;
	    ###########################################################
	    # Sometimes $entry->GNs->text returns names		      #
	    # enclosed in parentheses. e.g.: (SSX2 OR SSX2A) and SSX2B#
	    # Remove the parentheses.				      #
	    ###########################################################
	    $aa=~s/([\w\d])\;([\w\d])/$1 OR $2/g;
	    $aa=~s/\sand\s/ /i; 
	    $aa=~s/\(//g;
	    $aa=~s/\)//g;
	    my @names=(split(/\s+OR\s+/,$aa));

	    if($all_names){
		local $"="\t";
	    	print "$id\t@names$frag\n";
	    }
	    else{print "$id\t$names[0]$frag\n";}
	 
	}
	####################################################
        # If the entry has no GN line BUT it has a 	   #
	# CG1234-style FlyBase name, return that	   #
        ####################################################
	elsif($species == 7227 && $gene_name){
	   $ids{MAPPED}{$id}++;
	   print "$id\t$gene_name$frag\n";

	}

	#########################
	# Check AltNames        #
	#########################
	elsif(defined($entry->DEs->text) && not $strict_name){
	    #print "IDDDD : $id\n", Dumper($entry->DEs);
	    foreach my $de ($entry->DEs->elements){
		if($de->category eq 'RecName'){
		    print "$id\t" . $de->text   . "$frag\n" ;	
		    $ids{MAPPED}{$id}++;
		}
		elsif($de->category eq 'AltName'){
		    print "$id\t" . $de->text   . "$frag\n" ;	
		    $ids{MAPPED}{$id}++;
		}
		elsif($de->category eq 'SubName'){
		    print "$id\t" . $de->text   . "$frag\n" ;	
		    $ids{MAPPED}{$id}++;
		}
	    }
	    # my @names=split(/\s*OR\s*/,$entry->DEs->text);
	    # map{
	    # 	foreach my $de ($entry->DEs->elements){
	    # 	    defined($ids{WANT}{$de->text}) && do {
	    # 		$ret{$de->text}++; 
	    # 		$found++;
	    # 	    };
	    # 	}
	    # }@names;
	}
	## Some entries have no gene name
	return;  
    }
    elsif($toID eq 'ID'){
	print "$id\t", ${$entry->IDs->list}[0], "$frag\n";
	$ids{MAPPED}{$id}++;
    }
    else{
	##############
	## Check DR lines
	##############
	if($entry->DRs->get($toID)>0){
	    $ids{MAPPED}{$id}++;

	    foreach my $list ($entry->DRs->get($toID)){
		  print "$id";
		shift(@{$list});
		if($all_names){
		    print "\t@{$list}$frag\n" ;
		}
		else{
		    print "\t${$list}[0]$frag\n" ;
		}
	    }
	}
	else{
	    ## Sometimes not FBG Flybase ID is given as a FlyBase DR line BUT
	    ## such an ID is present anyway under EnsemblMetazoa
	    if($toID =~/flybase/i){
		if($entry->DRs->get("EnsemblMetazoa")>0){
		    foreach my $list ($entry->DRs->get("EnsemblMetazoa")){
			print "$id";
			shift(@{$list});
			if($all_names){
			    print "\t@{$list}$frag\n";
			}
			else{
			    map{
				if(/^FBgn/){
				    print "\t$_$frag\n";
				    $ids{MAPPED}{$id}++;
				    return();
				}
			    }@{$list};
			}
		    }
		}
	    }
	    
	}
    }
}
############################################################
############################################################
sub check_entry{
    my $entry=shift;
    my $id=${$entry->IDs->list}[0];
    my $from=shift;

    ## Sometimes, different query names correspond to
    ## the same entry. Therefore, use a hash to store
    ## the results instead of return()ing directly.
    my %ret; 
    my $found=0;
    my %null_hash=('bob' => 1);
    ## Do we know what we have been given?
    if($from ne 'unknown'){
	if($from eq 'GeneName') {
	    if(defined($entry->GNs->text)){
		## go through the requested names and
		## check if any of them are present in this entry
		# map{
		#     if ($entry->GNs->isPresent($_)){
		# 	$ret{$_}++;
		# 	$found++;
		#     }
		# }keys(%{$ids{WANT}});
		my @names=split(/\b\s+(OR|AND)\s+\b/i,$entry->GNs->text);
		##############################################
                # If we are also interested in synonyms (-s) #
                ##############################################
		if ($synonyms) {
		    map{
			###############################################################
                        # Sometimes $entry->GNs->text returns names		      #
			# enclosed in parentheses. e.g.: (SSX2 OR SSX2A) and SSX2B    #
			# Remove the parentheses.				      #
                        ###############################################################
			s/\(//g;
			s/\)//g;
			s/\s//g;
			#print "NN : \$ids{WANT}{$_} \n";
			if(defined($ids{WANT}{$_})){
			    $ret{$_}++;
			    $found++;
			}
		    }@names;
		}
		else {
		    if(defined($ids{WANT}{$names[0]})){
			    $ret{$names[0]}++;
			    $found++;
			}
		}
	    }
	    #########################
            # Check AltName	    #
            #########################
	    unless ($strict_name) {
		if(defined($entry->DEs->text)){
		    my @names=split(/\s*OR\s*/,$entry->DEs->text);
		    map{
			foreach my $de ($entry->DEs->elements){
			    defined($ids{WANT}{$de->text}) && do {
				print "AA ", $de->text ,"\n";
				$ret{$de->text}++; 
				$found++;
			    };
			}
		    }@names;
		}
	    }

	    $found>0 ? return(\%ret) : return(\%null_hash);
	}
	elsif($from eq 'AC'){
	    my @names=@{$entry->ACs->list};
	    map{
		$ret{$_}++  if(defined($ids{WANT}{$_}));
		## deal with refseq versions
		s/(.._\d+)\.\d+/$1/;
		##deal with genbank versions
		s/\.\d+//;
		defined($ids{WANT}{$_}) && do {
		    $ret{$_}++; 
		    $found++;
		};
	    }@names;
	     $found>0 ? return(\%ret) : return(\%null_hash);
	}
	elsif ($from eq 'PDB') {
	    ###########################################################
            # PDBs often map to multiple, unrelated proteins. 	      #
	    # Therefore, only return a PDB match if it is the	      #
            # ONLY match.					      #
            ###########################################################
	    foreach my $array ($entry->DRs->get($from)){
		shift(@{$array});
		map{
		    defined($ids{WANT}{$_}) && do {
			push @{$pdbs{$_}}, $entry;
		    };
		}@{$array};
	    
	    }
	}
	else{ 
	    foreach my $array ($entry->DRs->get($from)){
		## First element is DR class
		## eg, KEGG
		shift(@{$array});
		map{
		    defined($ids{WANT}{$_}) && do {
			$ret{$_}++;
			$found++;
		    };
		    ## deal with refseq versions
		    s/(.._\d+)\.\d+/$1/;
		    ## deal with genbank versions
		    s/\.\d+//;
		    defined($ids{WANT}{$_}) && do {
			$ret{$_}++;
			$found++;
		    };
		    ################################################################
                    # Deal with CG1912-style dmel names. They are sometimes	   #
		    # present in DR lines.     					   #
                    ################################################################
		    if($species == 7227 && /(CG\d{5,8})/){
		    	defined($ids{WANT}{$1}) && do {
		    	    $ret{$1}++;
		    	    $found++;
		    	};	
		    }
		    ##deal with mgi IDs
		    if($from eq "MGI"){
			my $ww=0;
			s/mgi://i;
			defined($ids{WANT}{$_}) && do {
			    $ret{$_}++;
			    $found++;
			};
		    }
		}@{$array}
	    }
	    $found>0 ? return(\%ret) : return(\%null_hash);
	}
    }	
    #############################################
    # If we do not know what we have been given #
    #############################################
    else{
	## Check for gene name
	if(defined($entry->GNs->text)){
	    my @names=split(/\s*OR\s*/,$entry->GNs->text);
	    map{
		defined($ids{WANT}{$_}) && do {
		    $ret{$_}++;
		    $found++;
		};
	    }@names;
	}
	#####################
	# Check AltName	    #
	#####################
	if(defined($entry->DEs->text)){
	    my @names=split(/\s*OR\s*/,$entry->DEs->text);
	    map{
		foreach my $de ($entry->DEs->elements){
		    defined($ids{WANT}{$de->text}) && do {
			$ret{$de->text}++; 
			$found++;
		    };
		}
	    }@names;
	}
	## foreach DR line
	foreach my $array ($entry->DRs->elements){
	    ## First element is DR class
	    ## eg, KEGG
	    shift(@{$array});
	    map{
		defined($ids{WANT}{$_}) && do {
		    $ret{$_}++;
		    $found++;
		};
		## deal with refseq versions
		s/(.._\d+)\.\d+/$1/;
		##deal with genbank versions
		s/\.\d+//;
		defined($ids{WANT}{$_}) && do {
		    $ret{$_}++;
		    $found++;
		};
		############################################################
		# Deal with CG1912-style dmel names. They are sometimes	   #
		# present as "DR   KEGG" lines:				   #
		# DR   KEGG; dme:Dmel_CG42703; -.	   	   	   #
		############################################################
		if($species == 7227 && /(CG\d{5,8})/){
		    defined($ids{WANT}{$1}) && do {
			$ret{$1}++;
			$found++;
		    };	
		}
	    }@{$array}
	}
    }
    $found>0 ? return(\%ret) : return(\%null_hash);
}

############################################################
############################################################

sub read_ids{
    my %ids;
    ## get desired ids
    if($opts{i}){
	$IDlist=$opts{i};
        ## Read desired IDs. Expecting a text file, one ID per line
	open(ID,$IDlist)||die("Cannot open ID list $IDlist :$!\n");
	while(<ID>){
	    next if /^\#/o;
	    next if /^\s*$/o;
	    chomp;
	    s/\s+$//o;
	    s/^\s+//o;
	    $ids{WANT}{$_}=$_;
	    $looking_for++;
	}
	die("Need a list of desired IDs. The file $opts{i} is empty!\n") unless scalar(keys(%{$ids{WANT}}))>0;
    }
    if($names){
	map{
	    $ids{WANT}{$_}=$_;
	    $looking_for++;
	}split(/\s+/,$names);
    }
    die("Need a list of desired IDs. Either as a text file (-i) or as a space separated list (-n)\n") unless scalar(keys(%{$ids{WANT}}))>0;

return(\%ids)
}
############################################################
############################################################

sub return_fasta{
    my $id=shift;
    my $entry=shift;
    my $real_id=shift; ## This is the ID I asked for
    my $l = length($entry->SQs->seq);
    my $i = 0;
    $ids{FOUND}{$id}++;
    $ids{MAPPED}{$id}++;
    $id .="|$real_id" unless $real_id eq 'bob';
    print ">$id \n";
    while($i<$l)   ## convert to FASTA
    {
	print substr($entry->SQs->seq,$i,60) . "\n"; 
	$i=$i+60;
    }
}
############################################################
############################################################

sub get_all_ids{
    my $entry=shift;
    my $ID=shift;
    ## foreach DR line
    foreach my $array ($entry->DRs->elements){
	## foreach DR element in DR line
	my $toID=${$array}[0];
	if(defined($ids{WANT}{$$array[1]})){
	    $ids{WANT}{$ID}=$$array[1];
	}
    }
    foreach my $array ($entry->ACs->list){
	if(defined($ids{WANT}{$$array[0]})){
	    $ids{WANT}{$ID}=${$array}[0];
	}

    }
}

############################################################
############################################################

sub usage{ 
    $0=~s/.+?([^\/]+)$/$1/;
    print STDERR <<EndOfHelp;

$0 parses UniProt Flat Files and returns either entire entries, or an entry\'s 
FASTA sequence or translates between ID types. For a list of available ID 
types, type : '$0 -I'. If no input ID type is specified, 'unknown' is assumed and $0 will search for EVERYTHING. If no return ID type is specified, UniProt accessions are returned. 

USAGE:  
        $0 [options] <FLAT FILE>

EXAMPLES:
  
    $0 -n CATA_HUMAN  -t AC human.flat
    $0 -n P04040  -t ID human.flat
    $0 -i names.txt -t AC human.flat
    $0 -Fn CATA_HUMAN human.flat


OPTIONS:
        -a : Return ALL names (inc. synonyms)
        -A : Search for matches in ALL names, including synonyms.
 	-e : Return entire entry in Flat File format.
	-f : What ID type the input IDs are ((f)rom).
	-F : Return FASTA sequence.
	-g : Strict gene names, only useful when giiving gene names as input. When this is set,
             the script will not check AltNames (DE) and will only return exact name matches 
             (GN lines, eg "GN Name=FDXR")
	-h : Print this help and exit.
	-i : List of IDs to look for.
	-I : Print the (long) list of IDs $0 can cope with.
	-m : Filename into which MISSING ids will be printed.
	-M : Filename into which the IDs that were found in 
	     the flat file but could not be mapped to the 
	     desired ID type will be printed.
	-n : What follows is a QUOTED, SPACE separated list of the names desired.
	-r : Allow protein fragments. If this is not set the script will skip any
             entries that describe protein fragments.
	-R : Mark protein fragments. If a fragment is found, the script will return
             <query_prot_name>\t'Fragment'
	-t : What ID type the input IDs should be converted (t)o.
	-s : Silent, do not print missed IDs to STDERR.
	-S : Spacer for STDERR indentation. Useful when $0 is called
	     by a wrapper script.

EndOfHelp
    exit();
}
############################################################
############################################################

sub print_ID_types{
    $0=~s/.+?([^\/]+)$/$1/;
    open(H, "| less");
    print H  <<EndOfHelp;
These are the various database identifiers that $0 can deal with. 
By default, if no other options are given, it will translate 
between UniProt accessions and ids. In addition to the options below,
'AC' and 'ID' are also valid FROM/TO values.

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
SGD                     :       S000005232
SMR			:	P04040
STRING			:	P04040
SUPFAM			:	SSF56634
SWISS-2DPAGE		:	P04040
UCSC			:	uc001mvm.1
UniGene			:	Hs.502302
WormBase                :       WBGene00003920
FlyBase                 :       FBgn0010339
MGI                     :       MGI:98578 || 98578
GeneName                :       par-5

AC                      :       P04637
ID                      :       P53_HUMAN
EndOfHelp
close(H);
exit();
}

