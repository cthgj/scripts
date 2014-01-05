#!/usr/bin/perl -w
#################################################################
# Attempt to remove redundancy. Retrieve FASTA seqs for         #
# each of the sequences, blast against each other and 	        #
# remove those those sequences with 100%id, where 	        #
# query!=sbjct and where the entire sequence of the subject     #
# is aligned perfectly against the query for its ENTIRE length. #
#################################################################
use strict;
use Getopt::Std;
use feature qw(switch say);
our $verbose;
our $debug;
our $force;
my (%trembls,%skip_frags,%files);
sub v_say;
require "MY_SUBS.pl";
my %opts;
getopts('vdhrn:f:F:N:s:b:M:T:t:S:',\%opts);
usage() if $opts{h};
usage() if $opts{h};
my $net_file=$ARGV[0] || usage();
my $flat=$opts{f}||undef;
my $tmp_dir=$opts{T}||"/tmp";
my $names_file=$opts{N}||"$tmp_dir/$$.names";
my $blast_out=$opts{b}||"$tmp_dir/$$.cdh";
my $map_file=$opts{M}||undef;
$verbose=$opts{v}||undef;
$debug=$opts{d}||undef;
my $swissprot_fasta=$opts{S}||"$tmp_dir/$$.sp.pep";
my $sim_threshold=$opts{t}||0.95;
my $discard_fragments=$opts{r}||undef;
my $fasta_file=$opts{F}||"$$.pep";

##########################################################
# If (a blast outfile has been given, skip to parsing it #
##########################################################
unless (-e $blast_out) {
    #########################
    # Collect protein names #
    #########################
    my $nf=check_file($net_file,"r") ;
    my %prots;
    while (<$nf>) {
	chomp;
	my @a=split(/\s+/);
	$prots{$a[0]}=$prots{$a[1]}=1;
    }
    open(A, ">$names_file");
    map{print A "$_\n"}keys %prots;
    close(A);

    ##################
    # Get FASTA seqs #
    ##################
    if (-e $fasta_file) {
	-z $fasta_file && do {die("FASTA seq file $fasta_file exists but is empty!")};
	v_say("Using existing FASTA seq file $fasta_file\n");
    } elsif ($flat) {
	
        ##############################################################
        # I am explicitly allowing fragments (-r) to be able 	     #
        # to map frag IDs to SWISSPROT	    			     #
        ##############################################################       
	v_say("Getting FASTA sequences from $flat");
	my $cmd;
	$verbose ?
	    ($cmd="uniprot_parse.pl -vf ID -rFi $names_file $flat > $fasta_file"): 
		($cmd="uniprot_parse.pl -f ID -rFi $names_file $flat > $fasta_file") ;
	debug($cmd);
	`$cmd`;
    } else {
	v_say("Retreieving seqs from UniProt...");
	my $c="uniprot_fetch.pl $names_file  >$fasta_file";
	debug($c);
	`$c`;
    }
#    unlink($names_file);

}


###############################
# Use CD-HIT instead of blast #
###############################
if (-e $blast_out) {
    -z $blast_out && do {die("CD-HIT outfile $blast_out exists but is empty!")};
    v_say("Using existing cd-hit outfile $blast_out\n");

    
} else {
    ########################################
    # Separate trEMBL prots from SWISSPROT #
    ########################################
    open(F, "FastaToTbl $fasta_file |")||die("Could not open FASTA file $fasta_file:$!\n");
 
    ###########################################################
    # Do not overwrite $swissprot_fasta if one has been given #
    ###########################################################
    unless ($opts{S}) {
	open(A, "| TblToFasta > $swissprot_fasta ")||
	die "Could not open $swissprot_fasta for writing\n"; 
    }
    open(B, "| TblToFasta > $tmp_dir/$$.tr.pep")||
	die "Could not open $tmp_dir/$$.tr.pep for writing\n"; 
    while (<F>) {
	/(([a-z0-9]+)_.+?)\s/i || die("Non-uniprot ID in the FASTA file $fasta_file: $_");
	my $nn=$2;
	my $name=$1;
	#######################################################
        # Check if it is a TrEMBL or SWISSPROT entry. 	      #
	# TrEMBLs are always <6 chars>_SPECIES.		      #
        #######################################################
	if (length($nn)==6) {
	    print B;
	    $trembls{$name}++;
	} else {
	    print A unless $opts{S};
	}
    }
    close(A); close(B); close(F);
    ##############
    # Run CD-HIT #
    ##############
#    my $cmd="cd-hit-2d -i $tmp_dir/$$.sp.pep  -i2 $tmp_dir/$$.tr.pep -o $blast_out  -c 0.98 -n 5 -d 100 -s 0.8; mv $blast_out.clstr $blast_out";
    my $cmd="cd-hit-2d -i $swissprot_fasta  -i2 $tmp_dir/$$.tr.pep -o $blast_out -c $sim_threshold -n 5 -d 100 ; mv $blast_out.clstr $blast_out";
    debug($cmd);
    `$cmd`;
    #perl -e 'my $c; my @a; while(<>){chomp; if(/^\>//){$#a>0 && do {print "$c : @a\n";}; $c=$_; @a=();} else{push @a, $_; }}' $blast_out 
}
#########
# Parse #
#########
open(A, "$blast_out")||die "Could not open cluster file $blast_out: $!\n"; 
debug("Opened $blast_out");
my $c; 
my @a=(); 
my %MAP;
while (<A>) {
    chomp;
    if (s/^\>//) {
	$#a>0 && do {
	    ## Get the SWISSPROT
	    map{
#		/sp.+\|(.+_.+)\./;
		if (defined($MAP{$_}) && $MAP{$_} != $a[0]) {
		    die("\$MAP{$_} is already set to $MAP{$_}, when attempting to set to $a[0]") 
		}
		$MAP{$_}=$a[0]; 
	    }@a;
	}; 
	$c=$_; 
	@a=();
    }
    else {
	## If the IDs were downloaded from UniProt
	if (/>..\|/) {
	    while (/>.+\|(.+_.+?)\./g) {
		push @a, $1;
	    }
	}
	## If the IDs were parsed from a local Flat file
	else {
	    />([^\s]+?)\./;
	    push @a, $1;
	}
    }
}
#####################
# Get the last line #
#####################
map{
    $MAP{$_}=$a[0]; 
}@a;
#######################
# Check for fragments #
#######################
if ($discard_fragments) {
    my $trfile="$tmp_dir/$$.tr";
    my $fraglist="$tmp_dir/$$.fr";
    my $trfh=check_file($trfile,"w");
    my $trnum=0;
    foreach my $tr (keys(%trembls)) {
	print $trfh "$tr\n";
	$trnum++;
    }
    ###########################################
    # If we have at least one trEMBL sequence #
    ###########################################
    if ($trnum>0) {
	my $cmd;
	$verbose ?
	    ($cmd="uniprot_parse.pl -f ID -Rvi $trfile $flat | gawk '\$3~/Fragment/'> $fraglist"): 
		($cmd="uniprot_parse.pl -f ID -Ri $trfile $flat | gawk '\$3~/Fragment/'> $fraglist") ;
	debug($cmd);
	`$cmd`;
    
	my $frfh=check_file($fraglist,"r");
	while (<$frfh>) {
	    chomp;
	    my @a=split(/\t/);
	    ########################################################
	    # Has the fragment already been mapped to a SWISSPROT? #
	    ########################################################
	    next if defined($MAP{$a[0]});
	    ###########################
	    # If it hasn't, remove it #
	    ###########################
	    $skip_frags{$a[0]}++;    
	}
    }
}


###############################################
# Now go through the network file and replace #
# any occurances of a mapped TrEMBL prot with #
# its equivalent SWISSPROT prot		      #
###############################################
v_say("Writing network\n");
my $nf1=check_file($net_file,"r") ;
my %seen;
my $skipped_frags=0;
while (<$nf1>) {
    chomp;
    my @a=split(/\t/);
    $MAP{$a[0]}=$a[0] unless defined($MAP{$a[0]});
    $MAP{$a[1]}=$a[1] unless defined($MAP{$a[1]});
    ###########################
    # Skip unmapped fragments #
    ###########################
    if ($discard_fragments && (defined($skip_frags{$a[0]}) || defined($skip_frags{$a[1]}))) {
	# defined($skip_frags{$a[0]}) ?
	#     debug("Skipping $a[0]") : 
	# 	debug("Skipping $a[1]") ;
	$skipped_frags++;
	next;
    }
    s/$a[0]/$MAP{$a[0]}/g;
    s/$a[1]/$MAP{$a[1]}/g;
    ## Now resplit, to sort pair and therefore 
    ## be able to remove duplicate interactions
    @a=split(/\t/);
    my $pair=join("\t", sort($a[0],$a[1]));
    $seen{$pair}++;
    print "$pair\n" unless $seen{$pair}>1;
}
debug("Skipped $skipped_frags fragments");
##################
# Print the %MAP #
##################
my @kk=keys(%MAP);
$#kk > -1 && do {
    my $fh=check_file($map_file,"w");
    map{
	unless ($_ eq $MAP{$_}) {
	    my $a=$_;
	    my $b=$MAP{$_};
	    ########################################
	    # Only print if a change has been made #
	    ########################################
	    print $fh "$a\t$b\n" unless $_ eq $MAP{$_};
	}
    }@kk;
};
###########
# Cleanup #
###########
#unlink()
#system("rm $$*");


sub usage{
    my $us="[dv] -f FLAT_FILE <NET_FILE> ";
    my $desc="This script will take a network file, retreive FASTA seqs for its proteins from a local flat file, and then blast the seequences against each other. Any sequence pairs that share 95% sequence identity, where query!=sbjct and where the entire sequence of the subject is aligned perfectly against the query for its ENTIRE length.";
    my %opts=(
	      "usage" => $us,
	      "desc" => $desc,
	      "n" => "Network file to analyze.",
	      "f" => "Flat file from which FASTA seqs will be extracted.",
	      "F" => "Fasta file containing the sequences of interest. If none is given, the sequences will be extracted from the flat file.",
	      "b" => "Blast output file name (def: <PID>.out)",
	      "N" => "File containing the protein names to be retrived. If none is given, the names are extracted from the network file.",
	      "M" => "Map file linking the original network name to that it was changed to in the non redundand network.",
	      "T" => "Temp. files will be created here (def: ./tmp)",
	      "r" => "Discard fragments. Any interactions involving a protein that is annotated as a fragment will be discarded unless that fragment has been mapped to a Swissprot protein, in which case its interactions are transferred to the Swissprot.",
	      "v" => "Verbose output.",
	      "d" => "Debugging output",
	      "t" => "Similarity threshold for cd-hit-3d. Sequences with less similarity will not be merged.",
	      "S" => "Swissprot FASTA file. This is the file that will be passed as db1 (-i) to cd-hit-2d, the one against which the sequences will be clustered. This is useful when comparing the network against the entire swissprot database, not only the other sequences of that are present in the network",
	      "h" => "Print this message and exit."
	     );
    print_help_exit(\%opts,0);
}


