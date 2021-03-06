#!/usr/bin/perl -w
use strict;
use Getopt::Std;
use Data::Dumper;



my $skipped_hits = 0;
my $skipped_hitsA = 0;
my $bit=0;
my $c = 0;
my $counter = 0;
my $score_line = undef;
my ($qhit_position, $is_aln_file, $curr_hit_position, $query_type, $subject_type, $blast_type, $length, $idval, $evalue, $score, $query_search_space, $bits, $identities, $subject_ID, $id, $info, $q_match, $hsp_start, $hsp_end, $middle_line,$previous_mid_line, $read_next_line,$hit_line,$HSP_length) = "";
my (%letters,%shit_position, %qhit_position, %hit_position, %subjects,%query,%subj,%hits,%scores, %species_hash,%evalue,%previous_mid_line, %hit_lines, %nohits, %wanted,%qqq, %unwanted);
my ( @smatches, @ids, @frames, @subjects, @hsps, @infos, @algnments,@hit,@q_matches,@species, @plants );
my $gotahit = 0;
my $frame = ""; ## initialise frame to empty string, in case working with blastn or blastp
my $query_counter = 0;
my $best_printed = 0;  ## yet another counter
my $specific_hits = 0; ## yes,still another counter
my $first_run_of_sub_gff = 1; ## counter to run smthng jus once at sub gff_output
my $specific = undef; ## Have we asked for a specific list of queries/subjects?
my $previous_line = 0; ## Used to check the previous line while reading the blast file;
my $query_ended = 10;
my $subj_ended  = 10;
my $printed_hits = 0;
my $prev_query = "tistheiassou";
my $query_ID = "tisdikiassou";
my %opts;
getopts('c:e:i:l:m:M:r:R:q:s:I:g:Q:T:f:u:U:pnFAXpVCSdvkhbBLx',\%opts) || do { print "Invalid option, try 'alignthingie.pl -h' for more information\n"; exit(1); };
&usage() if $opts{h};
unless ($ARGV[0])
{
    &usage("Need a blast file!!");
    
    exit(1);
    
}
if ($opts{b} && $opts{B}){&usage("ERROR : Options -b and -B are mutually exclusive");  exit(1);}

my $filename = $ARGV[0] || die("Cannot open file $ARGV[0] : $!\n");

#die("$filename.gff exists from a previous instance, please remove\n") if -r "$filename.gff";
### Set command line options and thresholds	 ###
my $want_all        = $opts{A} || undef;  # Want all HSPs which pass thresholds irrespective of aligned residues (*-* or whatever)
my $best            = $opts{b} || undef;  # Only return a single hit with the lowest evalue
my $single_best     = $opts{B} || undef;  # Return all hits which share the lowest evalue
my $CONS_ALL        = $opts{c} || 6;      # Minimum number of conserved residues (or '+') allowed  in the 2 $SECAROUND
my $check_for_C     = $opts{C} || undef;  # Check for aligned *-C as well
my $debug           = $opts{d} || undef;  # debugging info
my $EVALUE          = $opts{e} || 10;     # maximum evalue  >>
my $no_self         = $opts{f} || undef;  # No sel(f) : Skips subjects whose name matches (case-INsensitive) the value passed. (string)
my $no_same_seq     = $opts{F} || undef;  # Do not return hits against same sequence if present in database
my $gff             = $opts{g} || undef;  # Produce gff output for sub or query
my $IDMAX           = $opts{I} || 100;    # maximum    >>      >>
my $IDENTITY        = $opts{i} || 1;      #   >>    identity allowed
my $chck_cons       = $opts{k} || 1;      # Check conservation?
my $CONS_LEFT       = $opts{l} || 3;      # Minimum number of conserved residues (or '+') allowed on the left side of the matched residue
my $most_likely     = $opts{L} || undef;  # print most likely hits
my $SCOREMIN        = $opts{m} || 0;      # minimum score allowed
my $SCOREMAX        = $opts{M} || 10000;  # maximum   >>    >>
my $nohit_queries   = $opts{n} || undef;
my $no_plants       = $opts{P} || undef;
my $get_position    = $opts{p} || undef;
my $QMATCH          = $opts{q} || '\*';   # aa in the query seq.
my $query_list      = $opts{Q} || undef;  # Return only those HSP whose query is in $query_list (filename passed at cmndline)
my $CONS_RIGHT      = $opts{r} || 3;      # Minimum number of conserved residues (or '+') allowed on the right side of the matched residue
my $SECAROUND       = $opts{R} || 6;      # Length of region around matched residue to check for conservation ($SECAROUND-matched residue-$SECAROUND)
my $SMATCH          = $opts{s} || '\*';   # aa in subjct string 
my $strict          = $opts{S} || undef;  # Strict options, see if($strict) below
my $sbjct_list      = $opts{T} || undef;  # Return only those HSPs whose subject is in $sbjct_list (filename passed at cmndline)
my $unaligned_res   = $opts{u} || undef;  # How many unaligned residues are allowed (query length - hsp_length <= $unaligned_res)
my $UNWANTED        = $opts{U} || undef;  # List of queries I am NOT interested in
my $verbose         = $opts{v} || 0;      # Verbosity
my $no_same_species = $opts{X} || undef;
my $want_redox      = $opts{x} || 0;      # am I looking for redox boxes CXXC?


$verbose = 1 if $opts{V};
if($opts{u}){
    $unaligned_res = 0.1 if $opts{u} == 0;
}

# when you give "0" as cmd line par, perl takes it as empty and so, EVALUE==10

my $left_cons;
my $right_cons;
my $whohoo = 0;
# This is necessery for the functions
$best = 1 if $single_best;


if($gff){
die("****** ERROR ******\n-g option value must be \"s\" or \"q\"\n****** ERROR ******\n") unless $gff eq "s" || $gff eq "q";
}
if($strict)
{
    $EVALUE     =  0.01 unless $opts{e};
    $SCOREMIN   =  50 unless $opts{m};
    $IDENTITY   =  65 unless $opts{i};
    $CONS_LEFT  =  4 unless $opts{l};
    $CONS_RIGHT =  4 unless $opts{r};
    $SECAROUND  =  8 unless $opts{R};
    $CONS_ALL   =  10 unless $opts{c};
}


if ($opts{s})
{
 
    @smatches = split /\s+/, $SMATCH ;
}

else {push @smatches, $SMATCH;}
## Check for *-Cs as well
push @smatches, 'C' if $check_for_C;

### remove '\' from SMATCH for the eq operator. (usefull if pattern is '\*')
foreach my $ss (@smatches)
{
    if ($ss  =~ /\\/)
    {
	$ss =~ /\\(.)/;
	$ss = $1;
    }
}
 
#==================================================================#
# Define a list of species names (from ncbi's est database)	   #
# in order to avoid same-species hits				   #
#==================================================================#
&avoid_same_species() if $no_same_species;

## Get list of words to be avoided
my @no_selfs = split /\s+/, $no_self if $no_self;

#========================================================#
# Set $no_self to 1 to use the right subroutines 	 #
# if opts{F} has been given. Done after the 		 #
# above step since this option is to skip hits		 #
# against seqs with the same name and not 		 #
# avoid a specific key words as in $no_self		 #
#========================================================#
$no_self = 1 if $no_same_seq || $no_same_species || $no_plants;

if ($no_plants){&no_plants()}

## Get query or sbjct lists if present
if ($sbjct_list)
{
    $specific = 1;
    if (-e $sbjct_list )
    {
	open (A, "$sbjct_list");
	while(<A>)
	{
	    /^(.*)\s/;
	    $wanted{$1}++;
	}
	close(A);
    }
    else {$wanted{$sbjct_list}++; }
}
if ($query_list)
{
    $specific = 1;
    if (-e $query_list )
    {
	open (B, "<$query_list") || die("cannot open file $query_list\n");
	
	while(<B>)
	{
	    /^(.+?)\s?\n/;
#	    /^(.+?)\|/;
	    $wanted{$1}++;
	    my $a=$1;
	}
	close(B);
    }
    else {$query_list =~ /^(.+?)\s?$/ || die("cant match desired query\n"); $wanted{$1}++;}
}

# get list of unwanted queries
if($UNWANTED)
{
    if(-e $UNWANTED )
    {
	open (B,"$UNWANTED");
	while(<B>)
	{
	    /^(.+?)\s?\n/;
	    $unwanted{$1}++;	   
	}
    }
    else
    {
	my @tmp = split(/\s+/,$UNWANTED);
	map{$unwanted{$_}++}@tmp;
    }
}


my $qlength;
if ($QMATCH =~ /\\/){ my $qm = $QMATCH;  $qm =~ s/\\//; $qlength = length($qm); } 
else { $qlength = length($QMATCH)}

## Check that the options are all right
my @jj = keys(%opts);
my @numbers = ($SCOREMIN, $SCOREMAX, $EVALUE, $IDENTITY, $SECAROUND, $CONS_LEFT, $CONS_RIGHT, $CONS_ALL);

foreach my $number (@numbers)
{
  die("*** PROBLEM : \"$number\" is not a valid value, must use a number\ntry 'alignthingie.pl -h' for more information\n")  unless $number =~ /^[\d\.]+$/;
}

##  Do not check conservation
$chck_cons = 0 if $opts{k};
$chck_cons = 0 if $want_all;
## counts how many hsps hit
my $actual_hits = 0;


if ($debug)
{
    print STDERR "\$SCOREMIN      = $SCOREMIN\n";
    print STDERR "\$SCOREMAX      = $SCOREMAX \n";
    print STDERR "\$EVALUE        = $EVALUE\n";
    print STDERR "\$IDENTITY      = $IDENTITY\n";
    print STDERR "\$IDMAX         = $IDMAX\n";
    print STDERR "\$SECAROUND     = $SECAROUND\n";
    print STDERR "\$CONS_LEFT     = $CONS_LEFT\n";
    print STDERR "\$CONS_RIGHT    = $CONS_RIGHT \n";
    print STDERR "\$CONS_ALL      = $CONS_ALL\n";
    print STDERR "\$QMATCH        = \'$QMATCH\'\n";
    print STDERR "\$SMATCH        = \'$SMATCH\'\n";
    print STDERR "\@smatches      = @smatches\n";
    print STDERR "\$chck_cons     = $chck_cons \n";
    print STDERR "\$debug         = $debug\n";
    print STDERR "\$verbose       = $verbose\n";
    print STDERR "\$best          = $best\n" if $best;
    print STDERR "\$single_best   = $single_best\n" if $single_best;
    print STDERR "\$nohit_queries = $nohit_queries\n" if $nohit_queries;
    print STDERR "\$unaligned_res = $unaligned_res\n"  if $unaligned_res;
    print STDERR "\$no_same_seq   = $no_same_seq\n";
}



my $isblast = undef;
if ($filename =~ /\.gz$/) {
	open(FILE,"zcat $filename |");
    } else {
	open(FILE,"< $filename");
    };
my $skip_irrelevant_lines = 0;
my $skip_other_irrelevant_lines = 0;
my $first_hsp_for_sbjct = 0;
while(<FILE>)
{
    next if /^\s+$/;
    if($skip_irrelevant_lines==1)
     { 
 	unless(/^Query=/o)
 	{
 	    next;
 	}
     }
    if($skip_other_irrelevant_lines==1)
    { 
 	if(/^>/o || /^Query=/o)
 	{
 	   $skip_other_irrelevant_lines = 0;
 	}
	else{ next;}
     }

    unless($isblast)
    {
	die("***\nFile $filename is not a blast outfile!\n***\n") unless /BLAST/;
	if(/\#-([T]?BLAST[PNX]?)/)
	{
	    $is_aln_file = 1;
	}
	elsif(/([T]?BLAST[PNX]?)/)
	{
	    $is_aln_file = 0;
	}
	else{die("unkown file type\n");}
	$blast_type = $1;
	## trick so as I can run the script
	## on its own output
	print STDOUT "######################-$blast_type-#####################\n"; 
	if($get_position)
	{
	    if ($1 eq 'BLASTN'){
		$query_type   = 'd';
		$subject_type = 'd';
	    }
	    elsif($1 eq 'BLASTP'){
		$query_type   = 'p';
		$subject_type = 'p';
	    }
	    elsif($1 eq 'TBLASTN'){
		$query_type   = 'p';
		$subject_type = 'd';
	    }
	    elsif($1 eq 'BLASTX'){
		$query_type   = 'd';
		$subject_type = 'p';
	    }
	    elsif($1 eq 'TBLASTX'){
		$query_type   = 'd';
		$subject_type = 'd';
	    }
	    else{
		die("unknown blast flavor\n");
	    }
	   
	}
	
	    
	$isblast = 1;
    }  
      	
#=======================================#
# Get query name and empty	        #
# @algnments from previous query        #
#=======================================#
    if (/^Query=/)
    {	

	$query_ended = 0;
	$query_counter++;
	if ($verbose == 1)
	{  
	    print STDERR "." if $query_counter % 10 == 0; 
	    print STDERR "\t[$query_counter]\n" if $query_counter % 1000 ==0; 
	}	
	$whohoo = ".";
	
	if($query_ID)  ## If this is not the first query
	{
	    @algnments = ();
	}	
	/Query=\s(.*?)[\s]*\n/;
	$prev_query = $query_ID unless $query_ID eq "tisdikiassou"; ## unless this is the first, so no prev query
	$query_ID = $1;
	print "QQ : $query_ID\n";
	if($is_aln_file == 1 && $prev_query ne $query_ID && $prev_query ne "tistheiassou")
	{
	    if ($best){ &print_best_hit(@{$hits{$prev_query}}) if ${$hits{$prev_query}}[0]} 
	    else { &print_hits(@{$hits{$prev_query}})  if ${$hits{$prev_query}}[0]; }
	}
	$query_ID =~ /^(.+?)[\s\|]/;
#	$query_ID =~ /^(.+?)\|/;
	my $koko = $1;
	## Skip hits against same sequence
	## usefull when blasting the database seqs
	## against themselves and want to skip
	## self hits (opts{F})
	$query_ID =~ /^(.+?)\s/;
	print "QQ : $query_ID\n";
	$no_selfs[0] = $1 if $no_same_seq;
	&debug("###################################\nQuery : $query_ID\n");	
	## a trick to only look at the queries we are interested in 
	if($UNWANTED)
	{
	    if($unwanted{$query_ID})
	    {
		$skip_irrelevant_lines = 1;
		next;
	    }
	    else
	    {
		$skip_irrelevant_lines = 0;
	    }
	}
	else{$skip_irrelevant_lines = 0;}
	if($query_list)
	{
	    if ($wanted{$koko})
	    {
		$skip_irrelevant_lines = 0;
	    }
	    else
	    {
		$skip_irrelevant_lines = 1;
		next;
	    }
	}
	next;
    }
#===========================================#
# Get query length and put into hash	    #
# $letters = length		            #
#===========================================#
    if (/letters\)/)
    {
	/(\d+)/;
	$letters{$query_ID} = $1 || " ";
	$query_ended=1;
	$nohits{$query_ID} = 0;
	if($query_ID =~ /^(.+?)similar/i)
	{
	    $query_search_space = $1;
	}
	else {$query_search_space = $query_ID;}

	if($no_same_species)
	{
	    my $mm = 0;
	    for (my $i=0; $i<scalar(@species); $i++)
	    {
		if($query_search_space =~ /$species[$i]/i)
		{
		    $no_selfs[0] = $species[$i];
		    $mm=1;
		}
		$i=scalar(@species) if $mm>0;
	    }
	}
	next;
    }
    if($query_ended==0)
    {
	$query_ID = $query_ID . " " . $_;
	$query_ID =~ s/\s+/ /g;
	## initialise $nohits to 0, if it stays at 0 there are hits
	
    }
    if (/^>/) 
    {     
        ## @hsps = list of lines of last hsp
	## @alignments = list of hsps for last subject
	if($subject_ID)  ### if this is NOT the first sbjct for this query
	{	
	    push @algnments, [@hsps]; 
	    ## Make $a = unique name for this hsp
	    my $a = $subject_ID . $query_ID;   
	    $subjects{$a} = [@algnments]; 
	    &check_for_hits($a,"1") unless $nohit_queries;
	    @hsps = ();
	    @algnments=();
	   
	}

	#===================================================#
        # this is used so as if the evalue for this 	    #
	# HSP is already higher than the limits		    #
	# we know that all folloing subjects will	    #
	# also have a tpp-high evalue, and so, we 	    #
	# can directly jump to the next sbjct		    #
	# it is put to 0 at the score line.		    #
        #===================================================#
	$first_hsp_for_sbjct = 1;

	#===================================================#
        # Counter counts number of hsps per sub/quer	    #
	# so, set to 0 because we just got a new sub	    #
        #===================================================#
	$counter=0;
	$subj_ended = 0;
	chomp $subject_ID if $subject_ID;
	$subject_ID = $_;
	## @subjects = list of subject ids
	push @subjects, $subject_ID;
	next;	    
    }	
    ## self explanatory
    if (/No\shits\sfound/)
    {
	$nohits{$query_ID} = 1;
	next;
    }
    ## Get length 
    if (/Length\s?=/)
    {
	/=\s?(\d*)/ || die("Cannot match Length : $_\n");
	$length = $1;
	$unaligned_res = $length unless $opts{u};
	$subj_ended=1;
    }
    if($subj_ended == 0)
    {
	$subject_ID = $subject_ID . $_;
	$subject_ID =~ s/\s+/ /g;
    }

#========================================================================#
#  Get score line. The /Length/ and /Sbjct/ conditions			 #
#  are because for some reason, some EST deflines decide		 #
#  to include something like similar to : blah blah, Score = 100	 #
#  that screwed up the script, so I made it more specific		 #
#========================================================================#
    if((/^\s?Score\s+=/) && ($previous_line =~ /Length\s?=/ || $previous_line =~ /Sbjct/ ) )
    {  	
	&debug("Subject : $subject_ID");
	if ($counter>0)  ## if this is not the first HSP for this sub-quer pair
	{
	    ### usefull to give each hsp of this sub-quer pair unique name
	    my $c = $counter -1;
	    my $a = $subject_ID . $query_ID . $c;
	    push @algnments, [@hsps];   ## @algnments contains the list of lines of the last HSP (query, subject, middle)
	    @hsps = ();                   ## @hsps back to () in readiness for the next one.
	}
	## Counter counts number of hsps per sub/quer
	$counter++;
	$score_line = $_;
	chomp $score_line;
	/=\s+(.*?)\s.*=\s+([e\-\d]*)/ || die("score or eval problem $. :  $_\n"); 
	$score = $1;
	$evalue = $2;
	if ($evalue =~ /^\s?e/)
	{
	    $evalue = "1" . $evalue;
	}
	$skip_other_irrelevant_lines=1 if $evalue > $EVALUE;
	$skip_other_irrelevant_lines=1 if $score < $SCOREMIN;
	# if this is the first subject
	if($first_hsp_for_sbjct == 1)
	{
	    $skip_irrelevant_lines=1 if $evalue > $EVALUE;
	    $first_hsp_for_sbjct = 0;
	}
	## see line ~262 for explanation 
	$bit = 1;

	/^.*?=\s(.*?)\sbits/;
#	$bits = $1;  ### get bit score, am not using it at the moment....
	
    }
### Get Identities line
    if (/Identities/)
    {
	$identities = $_;
	chomp $identities;
	## this will select the id value irrispective of whether
	## this is an .out file or the result of a previous run
	## of alignthingie.pl. The line will only match "Score"
	## in the latter case.
	if (/Identities/)
	{
	    if($is_aln_file == 0)
	    {
		/Identities\s*=\s*\d+\/(\d+)\s*\(((\d{1,3})%)\)/ ||  die("identities problem1 (line $.) : $_ \n");
		$HSP_length = $1;
		$id = " Identities : " . $2;  ## Actual identities %age
		$idval = $3;
	    }
	    else
	    {
		/Identities\s*:\s*((\d+)%).*\s?(\d+)\n/ || die("identities problem2 (line $.) : $_ \n");
		$HSP_length = $3;
		$id = "$1";  ## Actual identities %age
		$idval = $2;
	    }
	    
	}
	$skip_other_irrelevant_lines=1 if $idval < $IDENTITY;
	next;
    }

### Get frame 

    if (/Frame/)
    {
	next if /MatID/;
	/Frame\s=\s(..)/;
	$frame = $_;#" Frame : " . $1 || ""; 
	next;
    }


       #====================================================================================#
       # This is the part where we parse the HSPs. That last $_ !~ /Query=/ is to 	    #
       # avoid problems that came up when attempting to run on blast results copied 	    #
       # from ncbis online blast remove if necessary...					    #
       #====================================================================================#

   
    if((/^Query/ || /^Sbjct/ || ($previous_line =~ /^Query/ )) && ($_ !~ /Score/) && ($nohits{$query_ID} == 0) && ($_ !~ /Query=/))
    { 
	my $c = $counter -1;
	my $a = $subject_ID ."**" . $query_ID . "**" . $c || die("$subject_ID ."**" . $query_ID . "**" . $c\n$_");
	#==============================================================#
        # Since this will loop through all the lines of the hsp	       #
	# and I only want to collect the scores once,		       #
	# I set bit=1 as soon as I get the score and only allow	       #
	# loop if bit is 1. At the end of the 1st iteration	       #
	# is 0 again so it will not continue.			       #
        #==============================================================#
	if ($bit == 1) 
	{
	    ## $info has all information for this HSP
	    if($is_aln_file == 0)
	    {
		$info = "(" . "Length = " . $length . ")\n\n\n" .  $score_line . "," . $id . " ===hsp : " . $HSP_length ."===". "\n" . $frame . "\n"  . "\n";
	    }
	    else
	    {
		$info = "(" . "Length = " . $length . ")\n\n\n" . $score_line . "\n" . $frame . "\n"  . "\n";
	    }
	    push @infos, $info;
	    die("counter : $counter\n") if $counter == 0;  ## $counter should never be 0 at this point.
	    
            #===================================================#
            # last HSP's scores and clear list		        #
	    # can access the scores by using 		        #
	    # unique key : scores{sub**quer**hsp_number}        #
            #===================================================#
	    $scores{$a} = [@infos];
	    @infos = ();
	    $bit = 0;
	}
	## Push current line => @hsps
	push @hsps, $_; 
	push @hsps, "\n" if /^Sbjct/;

	$previous_mid_line = "";
	
    }
    # If all we want is those queries with NO hits
    if ((/Matrix/) && ($nohit_queries))
    {
	next;
    }
    
    ## if we are at the end of the query's entry 
    ## and we have had at least 1 hit (and we want the hits)
    elsif((/Matrix/ || (/\#\#\#\#\#\#\#\#\#\#\#/ && $_ !~ /BLAST/)) && ($nohits{$query_ID} == 0)) 
    {
	## @algnments contains the list of lines of the last HSP (query, subject, middle)
	push @algnments, [@hsps];
	@hsps=();

	## make $a => unique name for quer/sub : subject+query
	my $a = $subject_ID . $query_ID;
	$subjects{$a} = [@algnments];
	&check_for_hits($a,"2");
	## %subjects contains all hsps 
	## for specific query/sub pair

	## make $a => unique name for specific hsp
	## so scores{a} will give all info for this hsp
	$a = $subject_ID ."**" . $query_ID . "**" . $c; 
	$scores{$a} = [@infos] unless exists $scores{$a};
	## empty alignments for next query
	@algnments = ();
	#===================================================#
        # If we have had at least one good 		    #
	# alignment, with qmatch and smatch aligned.	    #
	# See sub check_for_hits			    #
        #===================================================#
	if ($hit[0])
	{
	    #============================================#
            # @hit is the list of names			 #
	    # of those hsps of this sub/quer pair	 #
	    # which had a hit				 #
            #============================================#
	    if ($is_aln_file == 1) { 
		push @{$hits{$query_ID}},$hit[0] if $is_aln_file ==1;
	    }
	    else{	    
		if ($best){ &print_best_hit(@hit) if $hit[0]}
		else { &print_hits(@hit); }
	    }
	    ## Get ready for the next
	    @hit = ();
	    	    
	}

	## Else clear the hashes to save memorey usage
	else
	{
	    %subjects = () unless $is_aln_file == 1;
	    %scores =() unless $is_aln_file == 1;
	}
	$subject_ID = "";
    }   
    $previous_line = $_;

}
if($is_aln_file == 1)
{
    if ($best){ 
	&print_best_hit(@{$hits{$query_ID}})if defined($hits{$query_ID})} 
    else { &print_hits(@{$hits{$query_ID}})if defined($hits{$query_ID}); }
}
elsif ($hit[0])
{
    if ($best){ &print_best_hit(@hit) }
    else { &print_hits(@hit); }
}

#my @ahit = keys(%hits);

## number of sub/quer pairs with at least one hsp hit
#my $hit_number = @ahit;
print STDERR "[$query_counter] queries have been processed" if $verbose;

if ($nohit_queries)
{
    
    
    my $number_of_no_hit_queries = 0;
    print STDERR "\n";
    foreach my $query_name(keys(%nohits))
    {
	next if $nohits{$query_name} == 0;
	$number_of_no_hit_queries++;
	print STDOUT "$query_name\n";
    }

 
    print STDERR "$number_of_no_hit_queries sequences returned no HSPs\n";
    exit(1);
    
}

## No hits == no good... lets die.
if ($actual_hits == 0)
{
    print STDERR "\nSorry 0 hits found!\n";
    exit(1);
	
}


#================================================#
# Print the number of hits found/skipped 	 #
# depending on the options passed		 #
#================================================#

if ($best && $verbose)
{

    if($specific)
    {
	my @a = keys(%wanted);
	print STDERR "\n$actual_hits hits found, printed $specific_hits HSPs for requested query sequence(s) \n" if $query_list;
	print STDERR "\n$actual_hits hits found, printed $specific_hits HSPs for requested subject sequence(s) : @a\n" if $sbjct_list;
	exit(1);
    }

    print STDERR "\n$actual_hits hits found, $printed_hits hits printed\n";
    if ($no_self){
	$skipped_hits < 2 ? print STDERR "$skipped_hits self-hit ignored, " . ($actual_hits - $skipped_hits) . " HSPs printed \n" : print STDERR "$skipped_hits self-hits ignored, $best_printed hits printed (printing best hit for each query)\n" if $single_best;
	$skipped_hits < 2 ? print STDERR "$skipped_hits self-hit ignored, " . ($actual_hits - $skipped_hits) . " HSPs printed \n" : print STDERR "$skipped_hits self-hits ignored, $best_printed hits printed (printing best hits for each query)\n" unless $single_best;

	exit(0);
    }
    print STDERR "$best_printed hits printed (printing best hit for each query)\n" if $single_best;
    print STDERR "$best_printed hits printed (printing best hits for each query)\n" unless $single_best;
    exit(0);
}
elsif($specific)
{
    my @a = keys(%wanted);
    if($no_self)
    {
	print STDERR "\n$actual_hits hits found. $skipped_hits self-hits ignored and " . ($actual_hits - $skipped_hits) . " HSPs printed for requested query sequence(s) \n" if $verbose && $query_list;
    }
    else
    {
	print STDERR "\n$actual_hits hits found for requested query sequence(s)\n" if $verbose && $query_list;
	print STDERR "\n$actual_hits hits found for requested subject sequences\n" if $verbose && $sbjct_list;
    }
    
    exit(1);
}
elsif($verbose && $no_self)
{
    print STDERR "\n$actual_hits hits found, $skipped_hits self-hits ignored, " . ($actual_hits - $skipped_hits) . " HSPs printed \n";
    exit(1);
}
else
{
    print STDERR "\n$actual_hits hits found\n" if $verbose;
    exit(1);

}



############################################################################################
################                        FUNCTIONS                         ##################
############################################################################################

#=====================================================================================================#
# This is where it gets a little complicated... @ahit is the list of keys of %hits. 		      #
# So, this is a list of query names. We iterate through this list				      #
# using each query name to access the list of hits for that query. ( @{$hits{$query_key}}  ). 	      #
# Then foreach of those we print what we want to print.					              #
#=====================================================================================================#
sub print_hits
{
    
    my $kk=0;
    &debug("Printing hits...\n");
  hit:foreach my $hit (@_)
  {$kk++;
   unless($single_best){$best_printed++ if $best;}
   ## Split at ** to retrieve subject, query and number
   my ($SUB, $QUERY, $NUM) = split /\*\*/, $hit;
   &debug("\$SUB, \$QUERY, \$NUM : $SUB, $QUERY, $NUM");
   # print only desired hits
   if($specific)
   {
       my ($Q, $S) = ($QUERY, $SUB);
       if($query_list)	{$Q = $1 if $Q =~ /^(.+?)\|/; next hit unless defined $wanted{$Q}; }
       if($sbjct_list)	{$S =~ s/>//; $S = $1 if $S =~ /^(.*?)\s/;  next hit unless defined $wanted{$S};   }

   }
   $specific_hits++;
   
   # don't print self hits
   my $has_already_matched_another_of_the_selfs=0;	    
   if ($no_self)
   {
       &debug("noselfs :: @no_selfs");
       if($no_same_seq)
       {
	   $SUB =~ />(.+?)\s/;
	   my $kk = $1;
	   $QUERY =~ /^([^\s]*)\s*.*$/ || die();
	   
	   $1 eq $kk && do { $skipped_hits++ ; next hit;}
       }
       else
       {
	 self:foreach my $self (@no_selfs)
	 {	
	     if ($has_already_matched_another_of_the_selfs == 1){next hit;}
	     if ($SUB =~ /($self)/i)
	     {  
		 $skipped_hits++;
		 $has_already_matched_another_of_the_selfs=1;
		 if($kk < scalar(@_)){next hit}
		 else{return}
	     }
	     else
	     {
		 next self;
	     }
	 }
       }
   } 
   #==========================================================#
   # Print everything we need. Look at the comments for	      #
   # each of the data structures above and below	      #
   # to see what's going on				      #
   #==========================================================#
   my $a = $SUB . $QUERY;
#    if($get_position)
#    {
# #       ($SUB, $QUERY) = (&get_hit_position($SUB, $QUERY, $NUM, $hit));
#    }
   chomp $SUB;
   &debug("hit : $hit"); 
   $get_position ? print STDOUT "Query= $QUERY $qhit_position{$hit}\n\t($letters{$QUERY} letters)\n\n$SUB $shit_position{$hit}\n\t" : print STDOUT "Query= $QUERY\n\t($letters{$QUERY} letters)\n\n$SUB\n\t"; 
   ${$scores{$hit}}[0] =~ s/\s?===hsp\s:\s(.+\d+).*===/, $1/ || die("cannot match score : ${$scores{$hit}}[0] : $!\n") unless $is_aln_file == 1;
   &debug("\${\$scores{\$hit}}[0] : $hit");
   print STDOUT "@{$scores{$hit}}\n";
   map{print STDOUT "$_"; } @{$subjects{$a}[$NUM]}; 
   &gff_output(@{$subjects{$a}}[$NUM],@{$scores{$hit}},$QUERY,$SUB,$hit ) if $gff;
   print STDOUT "\n######################################################\n\n";
   $printed_hits++;
}
    ## Empty hashes to save memory usage
    %subjects = ();
    %scores =();
    
}
########################################################################################
########################################################################################

## Get the position of the matched residue in the
## query and/or subject sequences
sub get_hit_position
{
    &debug("getting hit position");
    my ($ff,$oo,$mm);
    my $hit_pos;
    my ($qhit_pos,$shit_pos);
    my ($hit, $qhit_position,$sbjct_line,$query_line,@array) = @_;
    my ($SUB, $QUERY, $NUM) = split(/\*\*/, $hit);
    my $k = 0;
    $query_line =~ /(\d+)\s*([^\d]*)\s*\d/ || die("cannot find query start position : $!\n");
    my $qstart = $1;
    my $qseq = substr($2,0,$qhit_position);
    my $qgaps = ($qseq =~ s/-//g);
    $sbjct_line =~ /(\d+)\s*([^\d]*)\s*\d/ || die("cannot find sbjctstart position : $!\n");
    my $sstart = $1;
    my $sseq = substr($2,0,$qhit_position);
    my $ss;
    my $sgaps = ($sseq =~ s/-//g);
    print STDOUT "\$qhit_position was $qhit_position\n"if $SUB =~ /TC202351/;
    my $shit_position = $qhit_position - $sgaps;
    $qhit_position = $qhit_position - $qgaps;
    print STDOUT "\$shit_position is $shit_position\n" if $SUB =~ /TC202351/;
    print STDOUT "\$qhit_position is $qhit_position\n"if $SUB =~ /TC202351/;
    print STDOUT "";
    if($blast_type eq 'TBLASTX')
    {
	map{
	    next unless /Frame/; 
	    /([+-]).*([+-])/; 
	    $ff = $1 . $2
	    } @{$scores{$hit}}; 
    }
    elsif($blast_type eq 'BLASTX' || $blast_type eq 'TBLASTN')
    {
	map{next unless /Frame/; /=\s?([-+])/; $ff = $1}@{$scores{$hit}}; 
    }
    else{$ff=0;}
    
    if ($blast_type eq 'BLASTN')
    {
	$qhit_pos =  "\(" . ($qstart + $qhit_position) . "\)";
	$shit_pos =  "\(" . ($sstart + $shit_position)  . "\)";
    }
    elsif($blast_type eq 'BLASTP')
    {
	$qhit_pos =  "\(" . ($qstart -1 + $qhit_position) . "\)" ;
	$shit_pos =  "\(" . ($sstart -1 + $shit_position)  . "\)";
    }
    elsif($blast_type eq 'BLASTX')
    {
	$shit_pos = "\(" . ($sstart -1 + $shit_position) . "\)" ;
	if($ff eq '-'){
	    $qhit_pos = "\(" . (($qstart +1) - ($qhit_position * 3)) . "\)";
	}
	else{
	    $qhit_pos = "\(" . (($qstart -3) + ($qhit_position * 3)) . "\)";
	}
    }
    elsif($blast_type eq 'TBLASTN')
    {
	$qhit_pos = "\(" . ($qstart -1 +  $qhit_position) . "\)";
	if($ff eq '-'){	    
	    $shit_pos = "\(" . (($sstart +3) - ($shit_position * 3)) . "\)";
	}
	else{				
	    $shit_pos = "\(" . (($sstart -3) + ($shit_position * 3)) . "\)" ;
	}
	
    }
    elsif($blast_type eq 'TBLASTX')
    {
	if($ff eq '--'){
	    $qhit_pos = "\(" . (($qstart +1) - ($qhit_position * 3)) . "\)" . "(($qstart -3) + ($qhit_position * 3))";
	    $shit_pos = "\(" . (($sstart +3) - ($shit_position * 3)) . "\) (($sstart -3) + ($shit_position * 3))\n";
	}
	elsif($ff eq '++'){
	    $qhit_pos = "\(" . (($qstart -3) + ($qhit_position * 3)) . "\)" . "(($qstart -3) + ($qhit_position * 3))" ;
	    $shit_pos = "\(" . (($sstart -3) + ($shit_position * 3)) . "\) (($sstart -3) + ($shit_position * 3))\n from : $sbjct_line\n" ;
	}
	elsif($ff eq '+-'){
	    $qhit_pos = "\(" . (($qstart -3) + ($qhit_position * 3)) . "\)" . "(($qstart -3) + ($qhit_position * 3))";
	    $shit_pos = "\(" . (($sstart -3) - ($shit_position * 3)) . "\)(($sstart -3) + ($shit_position * 3))" ;	  
	}
	elsif($ff eq '-+'){				  
	    $qhit_pos = "\(" . (($qstart +1) - ($qhit_position * 3)) . "\)" . "(($qstart -3) + ($qhit_position * 3))";
	    $shit_pos = "\(" . (($sstart -3) + ($shit_position * 3)) . "\)(($sstart -3) + ($shit_position * 3))" ;
	}
	else{die("cannot understand frames : *$ff*\n")}
    }
    else { die("unknown blast flavor : $!\n");}


    return($qhit_pos,$shit_pos);
}


########################################################################################
########################################################################################

## Print only the best hit (lowest eval) foreach query
sub print_best_hit
{
    &debug("Printing best hits...\n");
    my $evalue = 100;
    my $lowest_evalue = 100;
    my %best_hit;
    my @best_hits;
    my ($SUB, $QUERY, $NUM);
    my ($best_SUB,$best_NUM);
    my $kk=0;
  hit:foreach my $hit (@_)
  {
	$kk++;
	($SUB, $QUERY, $NUM) = split /\*\*/, $hit;
	my $has_already_matched_another_of_the_selfs=0;	    
	if ($no_self)
	{
	  self:foreach my $self (@no_selfs)
	  {
	    
	      next hit if $has_already_matched_another_of_the_selfs == 1;
	      if ($SUB =~ /$self/i)
	      {  
		  $skipped_hits++;
		  $has_already_matched_another_of_the_selfs=1;
		  if($kk < scalar(@_)){ next hit}
		  else{return}
	      }
	      else
	      { 
		  next self;
	      }
	  }
	}

	${$scores{$hit}}[0] =~ /Expect[^=]+=\s+(.*?),/ || die("Cannot match evalue(1) $hit:: ${$scores{$hit}}[0]  $.\n");
	$evalue = $1;
	if ($1 =~ /^\s?e/)
	{
	    $evalue = "1" . $evalue;
	}
	&debug("current evalue : $evalue lowest evalue : $lowest_evalue");
	## Debugging stuff, never mind
	if ($debug){$evalue < $lowest_evalue ? print STDERR "Smaller than smallest ? YES\n" : print STDERR "Smaller than smallest ? NO\n";}
	if($evalue <= $lowest_evalue)
	{
	    if($single_best)
	    {
		unless(defined($best_hit{$QUERY})){$best_hit{$QUERY} = $hit unless $has_already_matched_another_of_the_selfs;}
		$best_SUB = $SUB;
		$best_NUM = $NUM;
	    }	
	    push @best_hits, $hit;
	    $lowest_evalue = $evalue;
	}
	else
	{
	    next;
	}
    }
    if($single_best)
    {
	$best_printed++;
#	if (defined($best_hit{$QUERY})){die("$QUERY :: $best_hit{$QUERY}\n")}
	&print_hits($best_hit{$QUERY}) if defined($best_hit{$QUERY});
	my $a = $best_SUB . $QUERY;
 	&gff_output(@{$subjects{$a}}[$best_NUM],@{$scores{$best_hit{$QUERY}}},$QUERY,$best_SUB,$best_hit{$QUERY} ) if $gff;
    }
    else
    {   
	&print_hits(@best_hits);
    }
}

########################################################################################
########################################################################################

sub gff_output
{
    &debug("writing gff...");
  
    my @hsp = @{$_[0]};
    my $l = @hsp;

    $hsp[0] =~ /:\s(\d+)\s/ || die("No match for gff start, hsp line :\n$hsp[0]\n");
    my $query_start = $1;

    $hsp[$l-4] =~ /\s(\d+)\s?\n/;
    my $query_end = $1;

    my $QUERY = $_[2];
    $QUERY =~ s/\s+//g;

    my $SUBJECT = $_[3];
    $SUBJECT =~ s/>//;
    $SUBJECT = $1 if $SUBJECT =~ /^(.*?)\s/;

    $hsp[2] =~ /:\s?(\d+)\s*?\w/ || die("No match for gff start, hsp line :\n$hsp[2]\n");
    my $sbjct_start = $1;

    $hsp[$l-2] =~ /\s(\d+)\s?\n/ || die("No match for gff end, hsp line :\n$hsp[$l-2]\n");
    my $sbjct_end = $1;   

    my $info = $_[1];
    $info =~ /Frame\s+=\s+([-+])(\d+)/ || die("hoho : \n@$info\nargv : @_\n");
    my $frame = $2;
    my $strand   = $1;
    my $score    = '1';
    my $realname = '.';
    my $hit = $_[4];


    if($gff eq "s")
    {
	
	## Remove gff file from previous instance if there
	if ((-e "$QUERY.sbjcts.gff") && ($gff) && $first_run_of_sub_gff)
	{
	    print STDERR "$QUERY.sbjcts.gff exists from a previous instance, remove(y/n)?\n";
	    my $apan = <STDIN>;
	    chomp $apan;
	    if($apan eq "y"){unlink("$QUERY.sbjcts.gff");}
	    else{die("Please remove or rename $QUERY.sbjcts.gff to continue\n");}
	}

	open (GFF, ">>$QUERY.sbjcts.gff ") || die ("Cannot open $QUERY\.gff: $!");
	print  GFF  "$SUBJECT\talignthingie\tHSP\t$sbjct_start\t$sbjct_end\t$score\t$strand\t$frame\tQuery \"$QUERY\"\t.\n" if $strand eq '+'; 
	print  GFF "$SUBJECT\talignthingie\tHSP\t$sbjct_end\t$sbjct_start\t$score\t$strand\t$frame\tQuery \"$QUERY\"\t.\n" if $strand eq '-';
	close (GFF);

    }
    elsif($gff eq "q")
    {

	if ((-r "$QUERY.gff") && ($gff) && $first_run_of_sub_gff)
	{
	    print STDERR "$QUERY.gff exists from a previous instance, remove(y/n)?\n";
	    ## What??? Apantisi is answer in greek...
	    my $apan = <STDIN>;
	    chomp $apan;
	    if($apan eq "y"){unlink("$QUERY.gff");}
	    else{die("Please remove or rename $QUERY.gff to continue\n");}
	}
	
	open (GFF, ">>$QUERY.gff ") || die ("Cannot open $filename\.gff: $!");
	print  GFF  "$QUERY\talignthingie\tHSP\t$query_start\t$query_end\t$score\t$strand\t$frame\tTarget \"$SUBJECT\"\t.\n" if $strand eq '+'; 
	print  GFF "$QUERY\talignthingie\tHSP\t$query_end\t$query_start\t$score\t$strand\t$frame\tTarget \"$SUBJECT\"\t.\n" if $strand eq '-';
	close (GFF);
    }
   
     $first_run_of_sub_gff = undef;

}

########################################################################################
########################################################################################


sub check_for_hits
{
    my $key = $_[0];
    my $temp = $_[1];
    &debug("xxxxxxxxxxxxxxxxxxxxxx checking hits $temp xxxxxxxxxxxxxxxxxxxxxx");
    my $c = $counter -1;
    my $number = @{$subjects{$key}};
    if($want_all)
    {
	for (my $i=0; $i<$number;$i++)
	{
	    
	  array:	foreach my $array (@{$subjects{$key}}[$i] )
	  { 
	      my $ba = $subject_ID ."**" . $query_ID . "**" . $i;
	      &debug("\$a : $ba");
	      my $array_counter = 0;
	      my $prev_query="";
	      foreach my $hsp_line (@{$array})
	      {
		  $array_counter++;
		  next if $hsp_line =~ /^\s+$/;
		  &debug("Checking cons 1\n");
		  &check_cons($ba); 
		  if($gotahit)
		  {
		      $whohoo = "!";
		      print STDERR "!" if $verbose && $opts{V};
		      ## Count the hits
		      $actual_hits++;
		      my $a = $subject_ID ."**" . $query_ID . "**" . $i;
		      $hit_lines{$a} = $hsp_line;
		      push @hit, $a ;
		      $gotahit = 0;
		      next array;
		  }
	      }
	  }
	    
	}
     }
     else
     {	
	for (my $i=0; $i<$number;$i++)
	{
	  array:foreach my $array (@{$subjects{$key}}[$i] )
	  { 
	      my $ba = $subject_ID ."**" . $query_ID . "**" . $i;
	      &debug("\$a : $ba");
	      my $array_counter = 0;
	      my $prev_query="";
	       my $query_line;
	      foreach my $hsp_line (@{$array})
	      {
		  $array_counter++;
		  next if $hsp_line =~ /^\s+$/;
		 
		  if ($hsp_line =~ /^Query/) {  
		      &check_query($hsp_line,$array, $prev_query);
		      $prev_query= $hsp_line;
		      $query_line = $hsp_line;
		      }
		  elsif ($hsp_line =~ /^\s/)  {  &check_middle($hsp_line,$array,$ba);  }
		  else {  &check_subj($hsp_line,$array,$array_counter,$ba);  }		
		  if($gotahit)
		  {
		      $whohoo = "!";
		      print STDERR "!" if $verbose && $opts{V};
		      ## Count the hits
		      $actual_hits++;
		      my $a = $subject_ID ."**" . $query_ID . "**" . $i;
		      if ($get_position)
		      {
			  ($qhit_position{$a},$shit_position{$a}) = (&get_hit_position($a,$qhit_position,$hsp_line,$query_line,@{$array}));
		      }			  
#		      $hit_position{$a}{query} = $query_position;
#		      $hit_position{$a}{sbjct} = $sbjct_position;		      
		      $hit_lines{$a} = $hsp_line;
		      push @hit, $a ;
		      $gotahit = 0;
		      next array;
		  }
	      }
	  }
	    &debug("--------------- Checked for hits... --------------- ");
	    
	}
    }
}
########################################################################################
########################################################################################

sub check_query
{
    &debug("Checking query...\n \$hsp_line: $_[0]\n"); 
    @q_matches = (); 
    my @array = @{$_[1]};
    ## Am only interested in the aa/nt part of the line
    $_[0] =~ /y:?\s+\d+\s+(.*?)\s\d/;
    my $line = $1 or die("cannot match pattern : $_\n") ;
    ## used later on to get the right bit of the middle line
    $hsp_start = index($_[0], $line) ;
    $hsp_end = length($line);
    ## will remain -1 if no hits
    $q_match = -1;
    %qqq=();
    #=======================================================#
    # go through query line looking for $QMATCH    	    #
    # and push all matchihng positions to @q_matches	    #
    #=======================================================#
    while($line =~ /$QMATCH/ig)  
    {
	$q_match =  pos($line);
	$qqq{$q_match}++;
	if($most_likely)
	{
	    my $kl = $line;
	    if($want_redox == 1)
	    {
		if($q_match<3){ # if is too close to start of line, concat prev line
		    chomp($_[2]);
		    $kl = $_[2] . $line
		    }
		## find redox : CXXU
		my $redox = substr($kl, $q_match-4, $qlength);
		my $kkk=0;
		if ($redox eq "C")
		{
		    push @q_matches, $q_match;
		    $kkk=1;
		}
		## get UXXC as well
		if($kkk==0)
		{
		    if(length($kl)>$q_match+2)
		    {
			$redox = substr($kl, $q_match+2, $qlength);
			push @q_matches, $q_match if $redox eq "C";
		    }
		}
	    }
	    else
	    {
		push @q_matches, $q_match;
	    }
	    
	}
	else{	    
	    push @q_matches, $q_match;
	}
	
    }
}
&debug("Checked query.\n");
    

########################################################################################
########################################################################################

sub check_subj
{ 
    &debug("Checking subject...");
    my $s_match;
    my @array = @{$_[1]};
    my $array_counter = $_[2];
    my $a = $_[3];
    ## Am only interested in the aa/nt part of the line
    $_[0] =~ /\d+\s*([^\s]+?)\s/ || die("cannot match line :  $_[0]\n");
    my $line = $1;
    my %jj;  
    foreach my $sss (@smatches)
    {
	my $SS = quotemeta($sss);
	while($line =~ /$SS/ig)  
	{
	    $jj{pos($line)}++;
	}
    }
    my @kkk = keys(%jj);
    my @qqq = keys(%qqq);
    
    $q_match = -1;
    foreach my $j (@kkk)
    {
	$q_match = $j if defined($qqq{$j});
    }

#     if($unaligned_res)
#     {
# 	$array[0] =~ /:\s?(\d+)/;
# 	$Start = $1;
# 	for (my $i=scalar(@array)-1; $i--;)
# 	{
# 	    if($array[$i] =~ /Query:/)
# 	    {
# 		$array[$i] =~ /\s?(\d+)$/;
# 		$End = $1;
# 		$i=0;
# 	    }
# 	}
#     }
    
    ## If we found QMATCH, see what it alignes to
    if (($q_match != -1) && ($gotahit == 0))  ## do it only once
    { 
	unless($gotahit)
	{
	  Q: foreach my $q (@q_matches)
	  {
	      
	      &debug("q : $q");
	      next Q if $gotahit;
	      ## Get position aligned to matching position in query line
	      ## Use qlength to allow for longer patterns
	      $s_match =  substr($line, $q-$qlength, $qlength) || die ("match problem\n\$_  : $_[0] line : -$line-\nq : $q\nqmatches : @q_matches\n");

	    sbj:foreach my $subj_match (@smatches)
	    {
		if ($s_match eq $subj_match)
		{
		    $qhit_position = $q;
		    #$curr_hit_position = $q;
		    if($want_redox == 1)
		    {    
			my $redox = substr($line, $q-4, $qlength);
			my $redox1 = substr($line, $q+2, $qlength);
			next sbj unless ( $redox eq "C" || $redox1 eq "C");
		    }
		    if($most_likely)
		    {
			
			my $prev_subj_line = '';
			## concatenate all subject lines up to this one and discard if /*/
			if($array_counter > 2){ # if this is NOT the 1st sbjct line
			    my $ik = 2;
			    while($ik < $array_counter-1)
			    {
				$array[$ik] =~ /\d\s*?([^\d]+)\s+\d/ || die("arr : $array[$ik]");
				$prev_subj_line =  $prev_subj_line  . $1;
				$ik = $ik+4;
			    }
			}
			my $sub_upto_match =  $prev_subj_line . substr($line, 0, $q-1);
			next sbj if $sub_upto_match =~ /\*/;
		    }
		    &debug("matched (query)$QMATCH to $subj_match(sbjct) ");
		    if ($chck_cons == 0)
		    {
			&debug("Checking cons 2\n");
			&check_cons($a); 
			next Q;
		    }
		    else
		    {
			#=================================================================================#
			# These are the regions where we want to check for conservation :	          #
			# Immediately to the left and right of QMATCH and then in the surrounding area    #
			#=================================================================================#
			my ($left, $right, $region);
			#=========================================================================#
			# If the target QMATCH is too close to the beginning of the HSP's	  #
			# line, the left region to be checked for conservation will be too	  #
			# short, so we need to concatenate the region from the preceeding 	  #
			# line and then check the whole thing for conservation		          #
			#                                                                         #
			# The curent line is $array[$array_counter-1]                             #
			#=========================================================================#
			if ($q <=$SECAROUND )
			{  
			    #===============================================================================#
			    # If there isn't a previous middle line in the array, and			    #
			    # by extension no more conservation, goto next match 			    #
			    # (unless we are not checking conservation, in which case this is a hit)	    #
			    #===============================================================================#
			    
			    if ($array_counter-6 <0)
			    {  
				next Q if $chck_cons;
			    }
			    ## Else if there is...
			    else
			    {
				$previous_mid_line = $array[$array_counter-6];
				my $temp = substr($middle_line, 0, $q-1);
				my $len  = $SECAROUND-$q;
				$left = substr($previous_mid_line, -$len) . $temp;
			    }
			}
			else
			{
			    $left= substr($middle_line,$q-$SECAROUND-1,$SECAROUND );
			}
			$right = substr($middle_line, $q,$SECAROUND) || substr($middle_line, $q);
			$left_cons = 0;
			$right_cons = 0;
			my $region_cons = 0;
			#=========================================================================#
			# If the target QMATCH is too close to the end of the HSP's line	  #
			# the right  region to be checked for conservation will be too short      #
			# so we need to concatenate the region from the following line and 	  #
			# then check the whole thing for conservation			          #
			#=========================================================================#
			if (length($right) < $SECAROUND)
			{
			    my $temp_line;
			    if ($array[$array_counter+1])
			    { 
				#=============================================================#
				# Sometimes the blast output has 2 \n and sometimes 3.        #
				# This step is to get the right line			      #
				#=============================================================#
				if (($array[$array_counter+1]) && ($array[$array_counter+1] =~ /Query/))
				{
				    #==================================================================================#
				    # But what if, as has happened, the middle line is itself completely empty?	       #
				    # In that case the line following "query" will be "sbjct". This also	       #
				    # means that there is no conservation so we can go to next Q. Very special 	       #
				    # case but worth dealing with   						       #
				    #==================================================================================#
				    if($array[$array_counter+2] =~ /Sbjct/ )
				    {
					
					next Q;
				    }
				    else
				    {
					$temp_line = substr($array[$array_counter +2], $hsp_start );
#				      $temp_line = $array[$array_counter +2];
					&debug("tline : $temp_line\n");
				    }
				} 
				elsif (($array[$array_counter+2]) && ($array[$array_counter+2] =~ /Query/))
				{
				    $temp_line = $array[$array_counter +3];
				}
				## If this is the end of the HSP, next...
				else
				{
				    next Q;
				}
				#====================================================#
				# We now have the right line, so extract the         #
				# relevant regions and check for conservation	     #
				#====================================================#
				my $temp_hsp_start = index($array[$array_counter +2], $temp_line);
				my $temp_hsp_end = length($temp_line);
				my $new_line = $right . substr($array[$array_counter +2], $temp_hsp_start, $temp_hsp_end);
				my $aline = substr($array[$array_counter +2], $temp_hsp_start, $temp_hsp_end);
				$right = $right . substr($temp_line, 0, $SECAROUND - length($right));
				chomp $right;
				$region = $left . $right;
				&debug("l : -$left-\nr : -$right-\nregion : -$region-"); 
				while($left =~ /[A-Z\+]/g){$left_cons++}
				while($right =~ /[A-Z\+]/g){$right_cons++}
				while($region =~ /[A-Z\+]/g){$region_cons++}
			    }
			    # Special case, such as TRs, where TGA is just
			    # at the end of the query sequence. 
			    
			    elsif($letters{$query_ID} - $HSP_length <=4 )
			    {
				$right_cons = $region_cons = 30;	
				while($left =~ /[A-Z\+]/g){$left_cons++}
			    }
			    else
			    {
				next Q;
			    }

			}
			else
			{
			    $region = $left . $right;
			    &debug("l : $left\nr : $right\nregion : -$region-"); 
			    while($left =~ /[^\s]/g){$left_cons++}
			    while($right =~ /[^\s]/g){$right_cons++}
			    while($region =~ /[^\s]/g){$region_cons++}
			}
			
			
			# $previous_ine = $line;
			if($want_redox == 1)
			{    
			    my $redox = substr($line, $q-4, $qlength);
			    my $redox1 = substr($line, $q+2, $qlength);
			    next sbj unless ( $redox eq "C" || $redox1 eq "C");
			}
			&debug("Checking cons 3\n");

			&check_cons($left_cons,$right_cons, $region_cons,$a); 
		    }
		}
		## If $s_match is not equal to $SMATCH
		else
		{
		    next sbj;
		}
	    }
	  }
	}
    }
    
}
########################################################################################
########################################################################################
# &check_middle($hsp_line,$array,$ba)
sub check_middle
{
    &debug("Checking middle...\n \$_[0], \$hsp_start, \$hsp_end : -$_[0]-, -$hsp_start-, -$hsp_end-\n");
    $middle_line = substr($_[0], $hsp_start, $hsp_end);
}
########################################################################################
########################################################################################
sub check_cons
{
    my $missing = 10;
    if ($chck_cons == 1)
    {
	&debug("Checking conservation...");
	my ($left_cons,$right_cons,$region_cons, $a,$q) = @_;
	&debug("scores : ${$scores{$a}}[0]"); 
	${$scores{$a}}[0] =~ /.*Expect[^=]*=\s*(.*?)[\s\n,]/ || die("Cannot match evalue(2)\n ${$scores{$a}}[0], a : $a\n");  
	my $evalue = $1;
	if ($evalue =~ /^\s?e/)
	{
	    $evalue = "1" . $evalue;
	}
	${$scores{$a}}[0] =~ /Score\s=\s+(.*?)\s/;
	my $score =$1;
	${$scores{$a}}[0] =~ /Identities\s:\s+((.*?)%)/;
	my $idval = $2;
	$id = " Identities : " . $1;  ## Actual identities %age
	${$scores{$a}}[0] =~ /hsp\s:\s(\d+).*===/;
	my $hsp_length = $1;
	## how many query residues are not aligned?
	$opts{u} ?  $missing = abs($letters{$query_ID} - $hsp_length) :  $missing = 10 ;

	my @temp_array = ($left_cons,$right_cons,$region_cons, $score, $idval);
	foreach my $number (@temp_array)
	{	    
	    die("*** PROBLEM : left : $left_cons right : $right_cons region : score :  $score.\n $number is not a valid value, must use a number\ntry 'alignthingie.pl -h' for more information\n")  unless $number =~ /^[\d\.e+-]+$/;
	}
	&debug("eval $EVALUE >= $evalue, id : $IDMAX >= $idval >= $IDENTITY, score : $score\nlft cons : $left_cons >= $CONS_LEFT, rght cons : $right_cons >= $CONS_RIGHT");
	&debug("$missing <= $unaligned_res $a") if defined($unaligned_res);
	if ( ($left_cons >= $CONS_LEFT )  && ($right_cons >= $CONS_RIGHT) && ($region_cons >= $CONS_ALL) && ($score >= $SCOREMIN) && ($score <= $SCOREMAX ) && ($evalue <= $EVALUE ) && ($idval >= $IDENTITY) && ($idval <= $IDMAX) && ($missing <= $unaligned_res)) 
	{
	    if($most_likely)
	    {
		$gotahit++ unless $right_cons - $left_cons < -2;
	    }
	    else
	    {
		$gotahit++;
		&debug("reg cons  : $region_cons >= $CONS_ALL\nline $_\nxaxaxa left : $left_cons <= $right_cons");
	    }
	    
	}
	#	else {die("id : $idval, $evalue,$score, $a ")}
    }
    else
    {	&debug("Not checking conservation...");

	my $a=$_[0];
	${$scores{$a}}[0] =~ /Score\s=\s+(.*?)\s/;
	&debug("scores : ${$scores{$a}}[0]");
	my $score =$1;
	${$scores{$a}}[0] =~ /Expect[^=]+=\s+(.*?),/ || die("Cannot match evalue(3) ${$scores{$a}}[0], $a, $!\n");
	my $evalue = $1;
	${$scores{$a}}[0] =~ /hsp\s:\s(\d+).*===/|| die("${$scores{$a}}[0] :: $1 \n");
	my $hsp_length = $1;
	if ($evalue =~ /^\s?e/)
	{
	    $evalue = "1" . $evalue;
	}
	
	## how many query residues are not aligned?
	if(defined($opts{u}))
	    {
		$missing = abs($letters{$query_ID} - $hsp_length);
	    }
	else{$missing = 0;}
	${$scores{$a}}[0] =~ /Identities\s:\s+((.*?)%)/ || die("Regex problem(identities) sub check_cons\n ");
	my $idval = $2;
	$id = " Identities : " . $1;  ## Actual identities %age
	
	&debug("eval $EVALUE >= $evalue, id : $IDMAX >= $idval >= $IDENTITY, score : $score\n(missing)$missing <= $unaligned_res(unal)= abs($letters{$query_ID} - $hsp_length) $a");
	$gotahit++ if (($score >= $SCOREMIN) && ($score <= $SCOREMAX ) && ($evalue <= $EVALUE ) && ($idval >= $IDENTITY) && ($idval <= $IDMAX) && ($missing <= $unaligned_res));
    }

    
    &debug("...and got a hit!") if $gotahit;
} 
########################################################################################
########################################################################################
#===========================================================================#
# This allows us to check that the query and subject seqs do not	    #
# come from the same species. Here I just define the list of species	    #
# During execution I search for each of them in the query name		    #
# and then pass that as $no_selfs[0] to sub print_hits()		    #
#===========================================================================#
sub avoid_same_species
{
    @species = ("Bruguiera gymnorrhiza","Ctenocephalides felis", "Cryptomeria japonica", "Phytophthora sojae", "Schmidtea mediterranea", "Lolium perenne", "Beta vulgaris","Angiostrongylus cantonensis", "Haliotis asinina", "Isaria japonica", "Schistosoma japonicum", "Suaeda maritima", "Dioscorea nipponica","Humulus lupulus", "Fusarium oxysporum", "Mandarina ponderosa", "Xiphinema index", "Salicornia bigelovii", "Puccinellia tenuiflora","Julodis onopordi", "Oidium neolycopersici", "Meleagris gallopavo", "Heterodera avenae", "Oesophagostomum dentatum", "Clonorchis sinensis","Pseudotsuga menziesii", "Tricholepisma aurea", "Sorghum halepense", "Ambystoma tigrinum", "Agrostis capillaris", "Salix viminalis","Ambystoma mexicanum", "Centaurea maculosa", "Lutzomyia longipalpis", "Theromyzon tessulatum", "Spodoptera frugiperda", "Codonopsis lanceolata","Armillariella tabescens", "Ailuropoda melanoleuca", "Toxoptera citricida", "Moneuplotes crassus", "Lycoris longituba", "Chaetomium cupreum","Leucosporidium scottii", "Eragrostis pilosa", "Phaeostigma major", "Lycopersicon pimpinellifolium", "Apodemus sylvaticus", "Brugia pahangi","Leishmania infantum", "Gibberella circinata", "Malawimonas jakobiformis", "Periophthalmus modestus", "Montastraea faveolata", "Teladorsagia circumcincta","Phytophthora capsici", "Hordeum vulgare", "Hippocampus comes", "Cryptosporidium parvum", "Citrus reshni", "Poa secunda","Vigna unguiculata", "Vibrio vulnificus", "Anoplophora glabripennis", "Neospora hughesi", "Schedonorus arundinaceus", "Necator americanus","Crassostrea gigas", "Dictyocaulus viviparus", "Gekko japonicus", "Diploplastron affine", "Lytechinus variegatus", "Nematostella vectensis","Marsilea vestita", "Brassica juncea", "Scarabaeus laticollis", "Amorphotheca resinae", "Diospyros kaki", "Oncorhynchus mykiss","Micromalthus debilis", "Epidinium ecaudatum", "Coprinus cinereus", "Aeluropus lagopoides", "Haematococcus pluvialis", "Macrotyloma uniflorum","Anopheles gambiae", "Acropora tenuis", "Fusarium culmorum", "Datura innoxia", "Aegilops umbellulata", "Populus alba","Meloidogyne paranaensis", "Oryzias latipes", "Ascidia sydneiensis", "Opsanus beta", "Cycas rumphii", "Cicindela campestris","Babesia bovis", "Daphnia magna", "Cricetulus griseus", "Bacteroides fragilis", "Brachypodium distachyon", "Cyclamen persicum","Pinus pinaster", "Carica papaya", "Paracoccidioides brasiliensis", "Asterina pectinifera", "Pinus radiata", "Bombyx mori","Selaginella lepidophylla", "Nicotiana sylvestris", "Sclerotinia sclerotiorum", "Nilaparvata lugens", "Uromyces viciae-fabae", "Prunus armeniaca","Meloidogyne arenaria", "Mayetiola destructor", "Heterodera schachtii", "Guillardia theta", "Spartina alterniflora", "Prunus persica","Porteresia coarctata", "Glomus mosseae", "Paxillus involutus", "Penicillium chrysogenum", "Acetabularia acetabulum", "Pristionchus pacificus","Festuca arundinacea", "Mesocestoides corti", "Capsicum annuum", "Cladosporium fulvum", "Medicago sativa", "Ostertagia ostertagi","Drosophila erecta", "Metarhizium anisopliae", "Phaseolus coccineus", "Catla catla", "Antrodia cinnamomea", "Pyrococcus furiosus","Viola baoshanensis", "Phytophthora megakarya", "Pecten maximus", "Dicentrarchus labrax", "Blastocladiella emersonii", "Bufo gargarizans","Chaetomium globosum", "Oncometopia nigricans", "Dascillus cervinus", "Senecio aethnensis", "Ipomopsis aggregata", "Plasmodium yoelii","Xiphophorus helleri", "Geotria australis", "Ceratopteris richardii", "Ananas comosus", "Dirofilaria immitis", "Saccharomyces cerevisiae","Porphyra haitanensis", "Tribolium confusum", "Pinus patula", "Lemna gibba", "Ophiostoma piliferum", "Sogatella furcifera","Hydra vulgaris", "Taxus cuspidata", "Fragilariopsis cylindrus", "Sorghum propinquum", "Chiloscyllium plagiosum", "Perca fluviatilis","Hoplodactylus maculatus", "Elaphe quadrivirgata", "Pimephales promelas", "Crinipellis perniciosa", "Avena sativa", "Malus sieboldii","Crocus sativus", "Spermophilus lateralis", "Drosophila willistoni", "Pinctada fucata", "Olea europaea", "Amblyomma variegatum","Fenneropenaeus merguiensis", "Drosophila mojavensis", "Citrus temple", "Gigaspora margarita", "Vitis pseudoreticulata", "Hypocrea jecorina","Bubalus bubalis", "Polyplastron multivesiculatum", "Gigaspora gigantea", "Pratylenchus penetrans", "Oryza punctata", "Ajellomyces capsulatus","Zea mays", "Eschscholzia californica", "Ophiostoma piceae", "Lycopersicon pennellii", "Raphanus sativus", "Ictalurus punctatus","Prunus dulcis", "Streptocarpus dunnii", "Acorus americanus", "Mirabilis jalapa", "Macropus eugenii", "Helianthus petiolaris","Zinnia elegans", "Cyclopterus lumpus", "Plasmodium vivax", "Lactuca virosa", "Gossypium arboreum", "Phaseolus vulgaris","Populus tremuloides", "Camponotus festinatus", "Salicornia brachiata", "Guzmania lingulata", "Acropora cervicornis", "Nicotiana glauca","Diuraphis noxia", "Antheraea yamamai", "Poncirus trifoliata", "Clarias macrocephalus", "Trypanosoma cruzi", "Taraxacum kok-saghyz","Mycosphaerella graminicola", "Lycopersicon hirsutum", "Gnetum gnemon", "Persea americana", "Conidiobolus coronatus", "Carassius auratus","Zeldia punctata", "Leishmania donovani", "Pennisetum ciliare", "Aspergillus flavus", "Tuber borchii", "Thellungiella salsuginea","Caenorhabditis remanei", "Osmerus mordax", "Alternaria brassicicola", "Litopenaeus vannamei", "Dimocarpus longan", "Sulfurospirillum barnesii","Litopenaeus setiferus", "Rhagoletis pomonella", "Naegleria gruberi", "Cercospora nicotianae", "Carabus granulatus", "Hedyotis centranthoides","Drosophila yakuba", "Senecio vulgaris", "Anopheles darlingi", "Deinagkistrodon acutus", "Tagetes erecta", "Kandelia candel","Abutilon theophrasti", "Phlebiopsis gigantea", "Euphorbia tirucalli", "Hydra oligactis", "Ips pini", "Pan troglodytes","Sphodromantis centralis", "Lycopersicon esculentum", "Pratylenchus vulnus", "Petromyzon marinus", "Pleuronectes platessa", "Lonomia obliqua","Oryctolagus cuniculus", "Nippostrongylus brasiliensis", "Gallus gallus", "Paralichthys olivaceus", "Gregarina niphandrodes", "Oryza alta","Sarcocystis falcatula", "Griffithsia japonica", "Cunninghamella elegans", "Agaricus bisporus", "Zamia fischeri", "Sparus aurata","Zamia furfuracea", "Populus tremula", "Oncorhynchus tshawytscha", "Manduca sexta", "Meloidogyne chitwoodi", "Ixodes scapularis","Eucalyptus gunnii", "Thermomyces lanuginosus", "Pleurotus ostreatus", "Ornithodoros porcinus", "Suillus luteus", "Phaeodactylum tricornutum","Papio anubis", "Aphis gossypii", "Linum usitatissimum", "Meloidogyne hapla", "Helianthus paradoxus", "Bombus terrestris","Chamaecyparis formosensis", "Boltenia villosa", "Leucoraja erinacea", "Lactuca serriola", "Schistosoma haematobium", "Populus fremontii","Gigaspora rosea", "Globodera rostochiensis", "Saruma henryi", "Vitis riparia", "Vitis labrusca", "Amblyomma americanum","Hevea brasiliensis", "Bicyclus anynana", "Nicotiana sanderae", "Torpedo marmorata", "Glycine max", "Daphnia pulex","Takifugu rubripes", "Neurospora crassa", "Ciona intestinalis", "Timarcha balearica", "Cryphonectria parasitica", "Marchantia polymorpha","Piper longum", "Pyrus communis", "Populus euphratica", "Bruguiera sexangula", "Cynodon dactylon", "Periplaneta americana","Fungia scutaria", "Triticum aestivum/Thinopyrum", "Quercus robur", "Senecio cambrensis", "Antirrhinum majus", "Alexandrium catenella","Tamarix ramosissima", "Mycetophagus quadripustulatus", "Gadus morhua", "Diploptera punctata", "Aureobasidium pullulans", "Taeniopygia guttata","Ciona savignyi", "Malus domestica", "Populus deltoides", "Ricinus communis", "Sarcocystis neurona", "Nasonia vitripennis","Dysdera erythrina", "Trichuris muris", "Eubalaena glacialis", "Strongylocentrotus purpuratus", "Polyandrocarpa misakiensis", "Medicago truncatula","Entamoeba dispar", "Drosophila virilis", "Haliotis diversicolor", "Scherffelia dubia", "Glomus versiforme", "Solenopsis invicta","Pseudopleuronectes americanus", "Bos taurus", "Epichloe festucae", "Vitis shuttleworthii", "Plasmodium berghei", "Cronartium quercuum","Theobroma cacao", "Chaetoceros compressum", "Euprymna scolopes", "Stevia rebaudiana", "Hemicentrotus pulcherrimus", "Pinus banksiana","Eurydice pulchra", "Ulva linza", "Lolium temulentum", "Lachesis muta", "Acacia mangium", "Hydra magnipapillata","Terfezia boudieri", "Eisenia andrei", "Papaver somniferum", "Acropora palmata", "Coffea canephora", "Aquilegia formosa","Schistosoma mansoni", "Meladema coriacea", "Corynascus heterothallicus", "Plutella xylostella", "Lactuca sativa", "Anolis sagrei","Zingiber zerumbet", "Panax ginseng", "Culex pipiens", "Carpodacus mexicanus", "Physcomitrella patens", "Bombyx mandarina","Cordyceps bassiana", "Trichophyton rubrum", "Aphanomyces piscicida", "Globodera pallida", "Lilium longiflorum", "Helianthus annuus","Ipomoea nil", "Chenopodium quinoa", "Trichosurus vulpecula", "Cicindela littoralis", "Hypocrea lixii", "Anas platyrhynchos","Echinostoma paraensei", "Lupinus albus", "Prunus canescens", "Apium graveolens", "Euphorbia lagascae", "Trypanosoma brucei","Prymnesium parvum", "Ictalurus furcatus", "Oncorhynchus nerka", "Liriodendron tulipifera", "Strongyloides ratti", "Xenopus laevis","Leymus chinensis", "labyrinthulid quahog", "Mimulus guttatus", "Emiliania huxleyi", "Hypsibius dujardini", "Halocynthia roretzi","Schistocerca gregaria", "Ustilago maydis", "Hippoglossus hippoglossus", "Pisolithus tinctorius", "Copidosoma floridanum", "Solanum chacoense","Sarcoptes scabiei", "Vigna radiata", "Citrullus lanatus", "Limnephilus flavicornis", "Alexandrium tamarense", "Prorocentrum donghaiense","Prunus cerasus", "Opisthorchis viverrini", "Polygonatum sibiricum", "Capsicum frutescens", "Lentinula edodes", "Tetrahymena pyriformis","Celuca pugilator", "Diaprepes abbreviatus", "Callinectes sapidus", "Musca domestica", "Ascaris suum", "Secale cereale","Dreissena polymorpha", "Spermophilus tridecemlineatus", "Daucus carota", "Antheraea mylitta", "Theileria orientalis", "Tripsacum dactyloides","Litopenaeus stylirostris", "Haliotis discus", "Candida albicans", "Tineola bisselliella", "Senecio squalidus", "Choristoneura fumiferana","Aegilops speltoides", "Prototheca wickerhamii", "Drosophila pseudoobscura", "Xenopus tropicalis", "Vitis cinerea", "Perillus bioculatus","Venturia canescens", "Cervus nippon", "Trypanosoma carassii", "Pseudourostyla cristata", "Miamiensis avidus", "Eptatretus burgeri","Caenorhabditis japonica", "Cryptococcus laurentii", "Cupiennius salei", "Schizosaccharomyces pombe", "Citrus sinensis", "Vitis aestivalis","Phyllostachys edulis", "Heterocapsa triquetra", "Laminaria japonica", "Eragrostis tef", "Circulifer tenellus", "Mycobacterium smegmatis","Philodina roseola", "Gloeophyllum trabeum", "Hydractinia echinata", "Mentha piperita", "Volvariella volvacea", "Chenopodium amaranticolor","Kluyveromyces cicerisporus", "Marsilea quadrifolia", "Eurotium rubrum", "Glycine latifolia", "Terminalia arjuna", "Nicotiana benthamiana","Mycoplasma hyorhinis", "Robinia pseudoacacia", "Penaeus monodon", "Tortula ruralis", "Ancylostoma caninum", "Populus angustifolia","Gossypium herbaceum", "Neotyphodium coenophialum", "Capra hircus", "Vicia faba", "Convoluta roscoffensis", "Trichinella spiralis","Meloidogyne javanica", "Citrus jambhiri", "Coccidioides posadasii", "Paramisgurnus dabryanus", "Dendrocalamopsis oldhamii", "Fragaria vesca","Curculio glandium", "Oryza minuta", "Porphyra yezoensis", "Dictyostelium discoideum", "Capreolus capreolus", "Cucumis melo","Agrostis stolonifera", "Camellia sinensis", "Lupinus luteus", "Spizellomyces punctatus", "Citrus unshiu", "Drosophila ananassae","Silpha atrata", "Pfiesteria shumwayae", "Gloriosa superba", "Galleria mellonella", "Theileria annulata", "Felis catus","Pyrenophora graminea", "Paramecium tetraurelia", "Populus euramericana", "Tamarix hispida", "Dactylis glomerata", "Aplysia californica","Aerides japonica", "Festuca mairei", "Dermacentor andersoni", "Apis cerana", "Giardia intestinalis", "Cocos nucifera","Spadella cephaloptera", "Saprolegnia parasitica", "Alexandrium fundyense", "Alternaria alternata", "Citrus tangerina", "Macaca fascicularis","Pholiota nameko", "Lupinus angustifolius", "Tritrichomonas foetus", "Citrus paradisi", "Podocoryne carnea", "Triticum turgidum","Macaca nemestrina", "Chondrus crispus", "Digitaria sanguinalis", "Citrus aurantium", "Trichoderma hamatum", "Rhipicephalus appendiculatus","Trichuris vulpis", "Cicer arietinum", "Uncinula necator", "Branchiostoma floridae", "Elaeis oleifera", "Parelaphostrongylus tenuis","Canis familiaris", "Loa loa", "Plasmodiophora brassicae", "Solanum tuberosum", "Chironomus tentans", "Cucumis sativus","Iris hollandica", "Isotricha prostoma", "Enchytraeus japonensis", "Antonospora locustae", "Eptatretus cirrhatus", "Glycine clandestina","Lotus japonicus", "Hyaloperonospora parasitica", "Neotyphodium lolii", "Solanum americanum", "Phytomonas serpens", "Taraxacum officinale","Ochlerotatus triseriatus", "Echinococcus granulosus", "Brassica napus", "Vitis rupestris", "Castanea dentata", "Meriones unguiculatus","Quercus petraea", "Chlamydomonas reinhardtii", "Eimeria maxima", "Monosiga ovata", "Acropora millepora", "Triatoma brasiliensis","Heterorhabditis bacteriophora", "Scutellospora castanea", "Ocimum basilicum", "Emericella nidulans", "Xiphophorus maculatus", "Salvia miltiorrhiza","Agriotes lineatus", "Oryza grandiglumis", "Carcinoscorpius rotundicauda", "Athyrium distentifolium", "Picea sitchensis", "Ancylostoma ceylanicum","Gossypium hirsutum", "Heliocidaris erythrogramma", "Theileria parva", "Magnaporthe grisea", "Oikopleura dioica", "Globodera mexicana","Zantedeschia aethiopica", "Arabidopsis lyrata", "Anguilla japonica", "Helianthus argophyllus", "Drosophila simulans", "Cyprinus carpio","Dermacentor variabilis", "Melaleuca alternifolia", "Lumbricus rubellus", "Nicotiana tabacum", "Taiwania cryptomerioides", "Austrofundulus limnaeus","Pediculus humanus", "Boophilus microplus", "Blumeria graminis", "Monosiga brevicollis", "Solanum sogarandinum", "Mnemiopsis leidyi","Toxoplasma gondii", "Heliconius erato", "Pseudosciaena crocea", "Fusarium sporotrichioides", "Solanum brevidens", "Colletotrichum gloeosporioides","Mengenilla chobauti", "Bigelowiella natans", "Malus sieversii", "Oryza officinalis", "Panulirus argus", "Neospora caninum","Ipomoea batatas", "Aiptasia pulchella", "Meloidogyne incognita", "Laccaria bicolor", "Selaginella moellendorffii", "Leptosphaeria maculans","Oncorhynchus keta", "Oryza sativa", "Carya illinoinensis", "Solanum sparsipilum", "Geomyces pannorum", "Rosa chinensis","Illicium parviflorum", "Picea mariana", "Helix aspersa", "Euglena gracilis", "Biomphalaria glabrata", "Xerophyta humilis","Marmota monax", "Euclidia glyphica", "Nasonia giraulti", "Eleusine coracana", "Malva pusilla", "Bos indicus","Pongo pygmaeus", "Pyrocystis lunula", "Trifolium purpureum", "Gracilaria lemaneiformis", "Sclerotium cepivorum", "Andrographis paniculata","Dunaliella viridis", "Eriocheir sinensis", "Zingiber officinale", "Mastigamoeba balamuthi", "Deschampsia antarctica", "Pennisetum glaucum","Mentha arvensis", "Nicotiana attenuata", "Tursiops truncatus", "Acyrthosiphon pisum", "Boea crassifolia", "Urechis caupo","Coregonus clupeaformis", "Betula pendula", "Amborella trichopoda", "Seculamonas ecuadoriensis", "Pelodiscus sinensis", "Thlaspi caerulescens","Macaca mulatta", "Vigna angularis", "Strongyloides stercoralis", "Acanthoscurria gomesiana", "Gibberella moniliformis", "Aedes aegypti","Monascus pilosus", "Salmo salar", "Mesostigma viride", "Rhopalosiphum padi", "Tupaia belangeri", "Diaphorina citri","Diabrotica undecimpunctata", "Palaemonetes pugio", "Echinococcus multilocularis", "Papilio dardanus", "Gryllus bimaculatus", "Hordeum brevisubulatum","Lingulodinium polyedrum", "Hibiscus syriacus", "Botrytis cinerea", "Plasmopara halstedii", "Laodelphax striatellus", "Ovis aries","Pomatomus saltatrix", "Caenorhabditis elegans", "Anopheles funestus", "Laminaria digitata", "Agrocybe cylindracea", "Eleutherococcus senticosus","Puccinia striiformis", "Tribolium castaneum", "Perinereis aibuhitensis", "Caenorhabditis briggsae", "Panagrolaimus davidi", "Populus simonii","Pisum sativum", "Biphyllus lunatus", "Panicum miliaceum", "Plasmodium falciparum", "Sesbania rostrata", "Coprinopsis cinerea","Gillichthys mirabilis", "Paralichthys lethostigma", "Callosobruchus maculatus", "Manihot esculenta", "Myzus persicae", "Stomoxys calcitrans","Sesamum indicum", "Vitis vinifera", "Cryptococcus neoformans", "Simmondsia chinensis", "Streptocarpus rexii", "Cicindela litorea","Aspergillus niger", "Haemonchus contortus", "Oidiodendron maius", "Picea abies", "Larix kaempferi", "Argopecten irradians","Fundulus heteroclitus", "Hypocrea virens", "Mesobuthus gibbosus", "Poecilia reticulata", "Araneus ventricosus", "Ascaris lumbricoides","Anopheles stephensi", "Pneumocystis carinii", "Platystomus albinus", "Gossypium barbadense", "Eudiplodinium maggii", "Capsicum chinense","Elaeis guineensis", "Limonium bicolor", "Phaseolus acutifolius", "Hebeloma cylindrosporum", "Fagus sylvatica", "Trichoderma viride","Platichthys flesus", "Phytophthora parasitica", "Anopheles albimanus", "Tamarix androssowii", "Triticum timopheevii", "Lethenteron japonicum","Haementeria depressa", "Closterium peracerosum-strigosum-littorale", "Penicillium expansum", "Mangifera indica", "Pinus taeda", "Steinernema feltiae","Heterobasidion annosum", "Cichorium intybus", "Citrus clementina", "Botryllus schlosseri", "Arabidopsis halleri", "Planococcus lilacinus","Pyrus pyrifolia", "Labeo rohita", "Gracilaria gracilis", "Craterostigma plantagineum", "Musa acuminata", "Loligo pealei","Leishmania major", "Entodinium caudatum", "Galega orientalis", "Artemisia apiacea", "Nicotiana langsdorffii", "Glossina morsitans","Tetraodon fluviatilis", "Colletotrichum trifolii", "Phaeosphaeria nodorum", "Homo sapiens", "Lepeophtheirus salmonis", "Phytophthora infestans","Sesuvium portulacastrum", "Toxocara canis", "Cervus elaphus", "Dasytricha ruminantium", "Procambarus clarkii", "Juglans regia","Bambusa edulis", "Descurainia sophia", "Steinernema carpocapsae", "Yucca filamentosa", "Fucus vesiculosus", "Amphidinium carterae","Populus tomentiglandulosa", "Pyrocoelia rufa", "Pinctada maxima", "Ginkgo biloba", "Entamoeba histolytica", "Prosopis juliflora","Festuca rubra", "Parastrongyloides trichosuri", "Heterodera glycines", "Drosophila grimshawi", "Trifolium pratense", "Bothrops insularis","Prunus avium", "Welwitschia mirabilis", "Thellungiella halophila", "Catharanthus roseus", "Dunaliella salina", "Rana catesbeiana","Saccharum officinarum", "Dianthus caryophyllus", "Holothuria mexicana", "Puccinia graminis", "Solanum phureja", "Trichoderma harzianum","Solanum melongena", "Nannochloropsis oculata", "Diabrotica virgifera", "Panicum virgatum", "Jakoba bahamiensis", "Pacifastacus leniusculus","Branchiostoma belcheri", "Onchocerca ochengi", "Eucalyptus globulus", "Aegilops tauschii", "Acipenser transmontanus", "Citrus medica","Seriola quinqueradiata", "Siniperca chuatsi", "Dermatophagoides farinae", "Verticillium dahliae", "Plasmopara viticola", "Equus caballus","Citrus macrophylla", "Arabidopsis thaliana", "Glycine soja", "Pyrus ussuriensis", "Tetrahymena thermophila", "Vespertilio superans","Citrus reticulata", "Triticum monococcum", "Atriplex canescens", "Phanerochaete chrysosporium", "Leishmania mexicana", "Reticulitermes flavipes","Psoroptes ovis", "Thalassophryne nattereri", "Fragaria ananassa", "Homarus americanus", "Brugia malayi", "Melipona quadrifasciata","Paragonimus westermani", "Phytolacca americana", "Oreochromis niloticus", "Monascus aurantiacus", "Amanita muscaria", "Musa paradisiaca","Spilocaea oleaginea", "Chamaecyparis obtusa", "Mustela putorius", "Eimeria tenella", "Haplochromis chilotes", "Dugesia ryukyuensis","Mesembryanthemum crystallinum", "Moniezia benedeni", "Morus alba", "Triticum aestivum", "Coturnix japonica", "Coffea arabica","Fusarium graminearum", "Blomia tropicalis", "Adiantum capillus-veneris", "Apis mellifera", "Ichthyophthirius multifiliis", "Trichostrongylus vitrinus","Phytophthora nicotianae", "Cucurbita pepo", "Picea engelmannii", "Ageratum conyzoides", "Narcissus pseudonarcissus", "Karenia brevis","Acanthopanax sessiliflorus", "Gasterosteus aculeatus", "Melipona scutellaris", "Triphysaria versicolor", "Ixodes ricinus", "Mycobacterium tuberculosis","Lapemis hardwickii", "Artemia parthenogenetica", "Mesocricetus auratus", "Aspergillus fumigatus", "Puccinia coronata", "Heliconius melpomene","Glomus intraradices", "Ceratodon purpureus", "Lactuca perennis", "Homalodisca coagulata", "Alstroemeria peruviana", "Trichomonas vaginalis","Papilio xuthus", "Avicennia marina", "Sparus auratus", "Pinus elliottii", "Artemia franciscana", "Fucus distichus","Eucalyptus tereticornis", "Melittobia digitata", "Locusta migratoria", "Picea glauca", "Populus trichocarpa", "Pyrenophora teres","Ipomoea trifida", "Flammulina velutipes", "Acanthamoeba healyi", "Camponotus vittatus", "Ctenopharyngodon idella", "Gammarus pulex","Allium cepa", "Crassostrea virginica", "Lactuca saligna", "Rhynchosciara americana", "Fusarium solani", "Brassica rapa","Brassica oleracea", "Arabidopsis arenosa", "Rhododendron catawbiense", "Cavia porcellus", "Diplosoma listerianum", "Radopholus similis","Paracentrotus lividus", "Pinus sylvestris", "Allium sativum", "Populus nigra", "Piper colubrinum", "Euphorbia esula","Misgurnus anguillicaudatus", "Lymnaea stagnalis", "Cajanus cajan", "Plumbago zeylanica", "Wuchereria bancrofti", "Squalus acanthias","Lolium multiflorum", "Pimpinella brachycarpa", "Silene latifolia", "Euphorbia pulcherrima", "Phoenix dactylifera", "Marsupenaeus japonicus","Schizophyllum commune", "Dugesia japonica", "Pisolithus microcarpus", "Cirsium arvense", "Helicoverpa armigera", "Mus musculus","Metadinium medium", "Astatotilapia burtoni", "Monacrosporium haptotylum", "Hydrilla verticillata", "Gossypium raimondii", "Eucalyptus grandis","Eoxenoc laboulbenei", "Streptocarpus saxorum", "Eucalyptus nitens", "Gymnema sylvestre", "Pyronema omphalodes", "Gillichthys seta","Molgula tectiformis", "Limonium sinense", "Columba livia", "Ribes americanum", "Ornithorhynchus anatinus", "Bufo viridis","Gibberella zeae", "Dipylidium caninum", "Sphenodon punctatus", "Nicotiana megalosiphon", "Polistes canadensis", "Danio rerio","Hyalomma anatolicum", "Drosophila melanogaster", "Litomosoides sigmodontis", "Mytilus galloprovincialis", "Asparagus officinalis", "Trametes versicolor","Phakopsora pachyrhizi", "Sus scrofa", "Senecio chrysanthemifolius", "Culicoides sonorensis", "Aspergillus oryzae", "Arachis hypogaea","Sitophilus zeamais", "Anagasta kuehniella", "Haematobia irritans", "Epinephelus coioides", "Nuphar advena", "Carcinus maenas","Onchocerca volvulus", "Sorghum bicolor", "Leptinotarsa decemlineata", "Curcuma longa", "Salvelinus alpinus", "Vaccinium corymbosum","Hedyotis terminalis", "Monodelphis domestica" );

}
########################################################################################
########################################################################################
sub no_plants
{
    @plants = ("Abutilon theophrasti", "Acacia mangium", "Acanthopanax sessiliflorus", "Acetabularia acetabulum", "Acorus americanus", "Adiantum capillus-veneris", "Aegilops speltoides", "Aegilops tauschii", "Aegilops umbellulata", "Aeluropus lagopoides", "Aerides japonica", "Ageratum conyzoides", "Agrostis capillaris", "Agrostis stolonifera", "Allium cepa", "Allium sativum", "Alstroemeria peruviana", "Amborella trichopoda", "Ananas comosus", "Andrographis paniculata", "Anoplophora glabripennis", "Antirrhinum majus", "Apium graveolens", "Aquilegia formosa", "Aquilegia pubescens", "Arabidopsis arenosa", "Arabidopsis halleri", "Arabidopsis lyrata", "Arabidopsis thaliana", "Arachis hypogaea", "Artemisia apiacea", "Asparagus officinalis", "Athyrium distentifolium", "Atriplex canescens", "Avena sativa", "Avicennia marina", "Bambusa edulis", "Beta vulgaris", "Betula pendula", "Boea crassifolia", "Brachypodium distachyon", "Brassica juncea", "Brassica napus", "Brassica oleracea", "Brassica rapa", "Bruguiera gymnorrhiza", "Bruguiera sexangula", "Cajanus cajan", "Camellia sinensis", "Capsicum annuum", "Capsicum chinense", "Capsicum frutescens", "Carica papaya", "Carya illinoinensis", "Castanea dentata", "Catharanthus roseus", "Centaurea maculosa", "Ceratodon purpureus", "Ceratopteris richardii", "Chamaecyparis formosensis", "Chamaecyparis obtusa", "Chenopodium amaranticolor", "Chenopodium quinoa", "Cicer arietinum", "Cichorium intybus", "Cirsium arvense", "Citrullus lanatus", "Citrus aurantium", "Citrus clementina", "Citrus jambhiri", "Citrus macrophylla", "Citrus medica", "Citrus paradisi", "Citrus reshni", "Citrus reticulata", "Citrus reticulata ", "Citrus sinensis", "Citrus tangerina", "Citrus temple", "Citrus unshiu", "Closterium peracerosum-strigosum-littorale", "Cocos nucifera", "Codonopsis lanceolata", "Coffea arabica", "Coffea canephora", "Craterostigma plantagineum", "Crocus sativus", "Cryptomeria japonica", "Cucumis melo", "Cucumis sativus", "Cucurbita pepo", "Curcuma longa", "Cycas rumphii", "Cyclamen persicum", "Cynodon dactylon", "Dactylis glomerata", "Datura innoxia", "Daucus carota", "Dendrocalamopsis oldhamii", "Dermatophagoides farinae", "Deschampsia antarctica", "Descurainia sophia", "Dianthus caryophyllus", "Digitaria sanguinalis", "Dimocarpus longan", "Dioscorea nipponica", "Diospyros kaki", "Dunaliella salina", "Dunaliella viridis", "Elaeis guineensis", "Elaeis oleifera", "Eleusine coracana", "Eleutherococcus senticosus", "Eragrostis pilosa", "Eragrostis tef", "Eschscholzia californica", "Eucalyptus globulus", "Eucalyptus grandis", "Eucalyptus gunnii     ", "Eucalyptus tereticornis", "Euphorbia esula", "Euphorbia lagascae", "Euphorbia pulcherrima", "Euphorbia tirucalli", "Fagus sylvatica", "Festuca arundinacea", "Festuca maire", "Festuca rubra", "Fragaria ananassa", "Fragaria vesca", "Galega orientalis", "Gerbera", "Gerbera hybrid", "Ginkgo biloba", "Gloriosa superba", "Glycine clandestina ", "Glycine latifolia", "Glycine max", "Glycine soja", "Gnetum gnemon", "Gossypium arboreum", "Gossypium barbadense", "Gossypium herbaceum", "Gossypium hirsutum", "Gossypium raimondii", "Guzmania lingulata", "Gymnema sylvestre", "Haematococcus pluvialis", "Hedyotis centranthoides", "Hedyotis terminalis", "Helianthus annuus", "Helianthus argophyllus", "Helianthus paradoxus", "Helianthus petiolaris", "Helicosporidium sp.", "Hevea brasiliensis", "Hibiscus syriacus", "Hordeum brevisubulatum", "Hordeum vulgare", "Humulus lupulus", "Hydrilla verticillata", "Illicium parviflorum", "Ipomoea batatas", "Ipomoea nil", "Ipomoea trifida", "Ipomopsis aggregata", "Iris hollandica", "Juglans regia", "Kandelia candel", "Lactuca perennis", "Lactuca saligna", "Lactuca sativa", "Lactuca serriola", "Lactuca virosa", "Larix kaempferi", "Lemna gibba", "Leymus chinensis", "Lilium hybrid", "Lilium longiflorum", "Limonium bicolor", "Limonium sinense", "Linum usitatissimum", "Liriodendron tulipifera", "Lolium multiflorum", "Lolium perenne", "Lolium temulentum", "Lotus japonicus", "Lupinus albus", "Lupinus angustifolius", "Lupinus luteus ", "Lycopersicon esculentum", "Lycopersicon esculentum ", "Lycopersicon hirsutum", "Lycopersicon pennellii", "Lycopersicon pimpinellifolium", "Lycoris longituba", "Macrotyloma uniflorum", "Malus domestica", "Malus hybrid", "Malus sieboldii", "Malus sieversii", "Malva pusilla", "Mangifera indica", "Marsilea quadrifolia", "Marsilea vestita", "Medicago sativa", "Medicago truncatula", "Melaleuca alternifolia", "Mentha arvensis", "Mentha piperita", "Mesembryanthemum crystallinum", "Mesostigma viride", "Mimulus guttatus", "Mirabilis jalapa", "Morus alba", "Musa acuminata", "Musa paradisiaca", "Narcissus pseudonarcissus", "Nicotiana attenuata", "Nicotiana benthamiana", "Nicotiana glauca", "Nicotiana langsdorffii", "Nicotiana megalosiphon", "Nicotiana sanderae", "Nicotiana sp.", "Nicotiana sylvestris", "Nicotiana tabacum", "Nuphar advena", "Ocimum basilicum", "Olea europaea", "Oryza alta", "Oryza grandiglumis", "Oryza minuta", "Oryza officinalis", "Oryza punctata", "Oryza sativa", "Panax ginseng", "Panicum miliaceum", "Panicum virgatum ", "Papaver somniferum", "Pennisetum ciliare", "Pennisetum glaucum", "Persea americana", "Petunia", "Phaseolus acutifolius", "Phaseolus coccineus", "Phaseolus vulgaris", "Phoenix dactylifera", "Phyllostachys edulis", "Physcomitrella patens", "Phytolacca americana", "Picea abies", "Picea engelmannii", "Picea glauca", "Picea mariana ", "Picea sitchensis", "Pimpinella brachycarpa", "Pinus banksiana", "Pinus elliottii", "Pinus patula", "Pinus pinaster", "Pinus radiata", "Pinus sylvestris", "Pinus taeda", "Piper colubrinum", "Piper longum ", "Pisum sativum", "Plumbago zeylanica", "Poa secunda", "Polygonatum sibiricum", "Poncirus trifoliata", "Populus", "Populus alba", "Populus alba ", "Populus angustifolia", "Populus deltoides", "Populus euphratica", "Populus euramericana", "Populus fremontii", "Populus nigra", "Populus simonii", "Populus tomentiglandulosa", "Populus tremula", "Populus tremuloides", "Populus trichocarpa", "Porteresia coarctata", "Prosopis juliflora", "Prototheca wickerhamii", "Prunus armeniaca", "Prunus avium", "Prunus canescens", "Prunus cerasus", "Prunus dulcis", "Prunus persica", "Pseudotsuga menziesii", "Puccinellia tenuiflora", "Pyrus communis", "Pyrus pyrifolia", "Pyrus ussuriensis", "Quercus petraea", "Quercus robur", "Raphanus sativus", "Rhododendron catawbiense", "Ribes americanum", "Ricinus communis", "Robinia pseudoacacia", "Rosa", "Rosa chinensis", "Saccharum", "Saccharum officinarum", "Salicornia bigelovii", "Salicornia brachiata", "Salix viminalis", "Salvia miltiorrhiza", "Saruma henryi", "Schedonorus arundinaceus", "Scherffelia dubia", "Secale cereale", "Selaginella lepidophylla", "Selaginella moellendorffii", "Senecio aethnensis", "Senecio cambrensis", "Senecio chrysanthemifolius", "Senecio squalidus", "Senecio vulgaris", "Sesamum indicum", "Sesbania rostrata", "Sesuvium portulacastrum", "Silene latifolia", "Simmondsia chinensis", "Solanum americanum", "Solanum brevidens", "Solanum chacoense", "Solanum melongena", "Solanum phureja", "Solanum sogarandinum", "Solanum sparsipilum", "Solanum tuberosum", "Sorghum bicolor", "Sorghum halepense", "Sorghum propinquum", "Spartina alterniflora", "Stevia rebaudiana", "Streptocarpus dunnii", "Streptocarpus rexii", "Streptocarpus saxorum", "Suaeda maritima", "Tagetes erecta ", "Taiwania cryptomerioides", "Tamarix androssowii", "Tamarix hispida", "Tamarix ramosissima", "Tamarix sp.", "Taraxacum kok-saghyz", "Taraxacum officinale", "Taxus cuspidata", "Thellungiella halophila ", "Thellungiella salsuginea", "Theobroma cacao", "Thinopyrum intermedium", "Thlaspi caerulescens", "Tortula ruralis", "Trifolium pratense", "Trifolium purpureum", "Triphysaria versicolor", "Tripsacum dactyloides", "Triticum aestivum", "Triticum monococcum", "Triticum timopheevii", "Triticum turgidum", "Ulva linza", "Vaccinium corymbosum", "Vicia faba", "Vigna angularis", "Vigna radiata", "Vigna unguiculata", "Viola baoshanensis", "Vitis", "Vitis aestivalis", "Vitis cinerea", "Vitis labrusca", "Vitis pseudoreticulata", "Vitis riparia", "Vitis rupestris", "Vitis shuttleworthii", "Vitis vinifera", "Welwitschia mirabilis", "Xerophyta humilis", "Yucca filamentosa", "Zamia fischeri", "Zamia furfuracea ", "Zantedeschia aethiopica", "Zea mays", "Zingiber officinale", "Zingiber zerumbet", "Zinnia elegans");

    foreach my $pp (@plants)
    {
	push @no_selfs, $pp;
    }

}
########################################################################################
########################################################################################

sub usage
{   
    open(HELP, "| more") ;
    print HELP <<EndOfHelp;
PROGRAM: 
        alignthingie.pl

USAGE:  
	alignthingie.pl [options] blast.out 

alignthingie will parse a blast outfile and return those HSPs
where a specific residue in the query line is found aligned to
another (or the same) specific residue in the subject line if
that hsp meets certain evalue, identity, score and conservation 
criteria.

COMMAND-LINE OPTIONS:

    -c : Minimum number of conserved residues (or '+') allowed around the matched residue (integer, def : 6)
    -C : Also check for Cs in the subject sequence which align to a '*' in the query
    -e : Maximum e-value allowed (integer, def : 10)
    -i : Minimum (i)dentity percentage allowed (integer, def : 0)
    -I : Maximum (I)dentity percentage allowed (integer, def : 100)
    -l : Minimum number of conserved residues (or '+') allowed on the (l)eft side of the matched residue
         (def : 3)
    -r : Minimum number of conserved residues (or '+') allowed on the (r)ight side of the matched residue
         (def : 3)	 
    -M : (M)aximum score value allowed (integer, def : 10000)
    -m : (m)inimum score value allowed (integer, def : 0)
    -q : String to match in query (def : '\*' )
    -s : String to match in subject     (def : '\*')
    -R : Length of (R)egion around matched residue to check for conservation (integer, def : 6)
         (R-matched residue-R)	 
    -S : Use strict evalue, id, score (0.01,65,50 respectively) and conservation cuttoffs.
         Writing the cons cuttoffs too long, just run "alignthingie.pl -Sd" to see...
    -k : Do not check conservation
    -u : How many (u)naligned residues are allowed with respect to query length.  
         (<query length> - <hsp length> <= <value passed>)

OUTPUT OPTIONS:
    -A :  Return (A)ll HSPs which pass thresholds without looking for any specific aligned residues
    -b :  Print only the (b)est (lowest e-value) hit for each query. If the smallest evalue
          is shared by more than one HSP, all such HSPs will be printed.
    -B :  Print only the (B)est (lowest e-value) hit for each query. If the smallest evalue
          is shared by more than one HSP, only the first such HSP will be printed.
    -d : (d)ebugging mode, very very verbose...
    -f :  No sel(f) : Skips subjects whose name matches (case-INsensitive) the value passed. (string)
    -F :  Generalised no sel(F), takes first chars (till first space) of query and subj names and 
          skips the hit if the 2 are identical.
    -g : (g)ff output. Use "-g q" for guery position gff and "-g s" for subject gff.
    -L :  Print most (L)ikely hits. ie, those with no stop codon before the matched residue and whose
          conservation on the right side of the match is no more than 2 less than that of the left side.
	  (this option is really only usefull for selenoprotein searches)
    -n :  Print only the names of those queries which did NOT return any HSP.
    -p :  Do not return hits against (p)lant species.
    -Q : (Q)uery name or list of names (text file, one name per line) to return HSPs for. Only those HSPs whose 
          query is specified will be printed
    -T : (T)arget (subject) name or list of names (text file, one name per line) to return HSPs for. 
          Only those HSPs whose subject is specified will be printed
    -v : (v)erbose output, prints a . for each query processed.
    -V :  more (V)erbose output, prints a . for each query processed and a ! for each hit found.
    -x : Only return those hits with a redox box CXXU/*
    -X : Read a list of species names and ignore hits between the same species.
    -U : Query name or list of (U)nwanted queries. Quoted list of query names (or text file, one name per line)
         for which NOT to return HSPs.
    
EndOfHelp
close(HELP);

    print STDERR "\n***** @_ *****\n\n" if @_;
    exit(1);

}

########################################################################################
########################################################################################


sub debug
{
    if ($debug)
    {
	print STDERR "@_\n";
    }
}
