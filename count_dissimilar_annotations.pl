#!/usr/bin/perl -w
#################################################
# This script will parse a .gaf file and print  #
# the number of dissimmilar annotations (PoGo)  #
# found for each unique protein		        #
#################################################
use strict;
use Getopt::Std;
use DBI();
use autodie;
require "MY_SUBS.pl";
my %opts;
getopts('dvhmieM:T:H:G:o:s:t:',\%opts);
usage() unless $ARGV[0];
usage() if $opts{h};
my $gaf_annotations_file=$ARGV[0];
my $subonto=$opts{o}||'P';
my $genealogy_file=$opts{G}||"$ENV{HOME}/research/new_moon/data/data/all_three.genealogy";
our $debug=$opts{d}||undef;
my $verbose=$opts{v}||undef;
my $exp_codes=$opts{e}||undef;
my $threshold=$opts{t}||0.05;
my $no_multi=$opts{m}||undef;
my $species=$opts{s}||'human';
my $table=$species;
my $db_port=3306;
my $db_host=$opts{H}||"10.1.3.30";
$db_host="127.0.0.1" if $db_host =~ /^l/i;
$db_port=3307 if $db_host eq "127.0.0.1";
my $map_file=$opts{M}||undef;
my (%proteins,%prob,%offspring,%ancestors);
my $TESTS=0;
 my ($dbh,$dsn);
if ($opts{i}) {
    $table="inter_" . $table;
}
my %map;
if ($map_file) {
    open(my $fh, '<', $map_file);
    while (<$fh>) {
	chomp;
	my @a=split(/\s+/);
	$map{$a[0]}=$a[1];
	$map{$a[1]}=$a[0];
    }
}
#################################
# Get the "good" evidence codes #
#################################
my %good_codes=
(
 "IDA"=>1,
 "IEP"=>1,
 "IGI"=>1,
 "IMP"=>1,
 "IPI"=>1,
 "TAS"=>1
);


parse_gaf_file();
my %prot_probs;
my $count=0;
my $tt=scalar(keys(%proteins));
foreach my $prot (keys(%proteins)) {
    $count++;
    print STDERR "$count of $tt\r" if $verbose;
    my @GOs=keys(%{$proteins{$prot}{ALL}});
    for (my $i=0; $i<=$#GOs; $i++) {
	for (my $k=$i+1; $k<=$#GOs; $k++) {
	    $TESTS++;
	    ## Colect the probabilities
	    push @{$prot_probs{$prot}},gopair_probability($GOs[$i],$GOs[$k],$species,$prot);
	    # my $gopair=join("_", sort {$b lt $a} ($GOs[$i],$GOs[$k]));
	    # $prot_probs{$prot}{$gopair}=gopair_probability($GOs[$i],$GOs[$k],$species,$prot);
	}
    }
}
print "\n" if $verbose;
## Do not correct for multiple testing
if ($opts{T}) {
    $TESTS=$opts{T};
}
$TESTS=1 if $no_multi;
## Now, see the dissimilar ones
my $prot_tot=0;
foreach my $prot (keys(%prot_probs)) {
    my $tot=0;
    # foreach my $gopair (keys(%{$prot_probs{$prot}})) {
    # 	my $corr= $prot_probs{$prot}{$gopair} * $TESTS;
    # 	$tot++ if $corr<=$threshold;
    # 	print "$gopair\t$corr\n" if $corr<=$threshold;
    # }
    foreach my $p (@{$prot_probs{$prot}}) {
    	my $corr= $p * $TESTS;
    	$tot++ if $corr<=$threshold;
    }
    print "$prot\t$tot\n";
    $prot_tot++ if $tot>0;
}
print STDERR "$prot_tot proteins had at least 1 pair of dissimilar $subonto annotations\n" if $verbose;
############################################################
sub parse_gaf_file{    
    print STDERR "Parsing GAF file...\n" if $verbose;
    open(A,"$gaf_annotations_file")|| die("Cannot open $gaf_annotations_file:$!\n");
    while(<A>){
	next if /^!/;
	chomp;
	my @tmparray=split(/\t/,$_);
	next unless $tmparray[8] eq $subonto;
	if ($map_file) {next unless defined($map{$tmparray[1]})};
	next if $tmparray[3] !~/^$/;
	##############################
        # Skip ontology roots	     #
        ##############################
	if ($tmparray[4] =~/0008150/ ||
	    $tmparray[4] =~/0005575/ ||
	    $tmparray[4] =~/0003674/){
	    next;
	}
	#########################################
	# If we only want "good" evidence codes #
	#########################################
	if ($exp_codes) {
	    next unless defined($good_codes{$tmparray[6]});
	}
	my $name=$tmparray[1];
	if ($map_file) {
	    $name=$map{$name};
    	}
	$proteins{$name}{ALL}{$tmparray[4]}++;
    }
    close(A);
}
############################################################
sub gopair_probability{
    my $go1=shift;
    my $go2=shift;
    my $table=shift;
    my $prot=shift;
    my $database='pogo';
    my @b = sort {$b lt $a} ($go1,$go2);
    my @c = sort {$a lt $b} ($go1,$go2);
    my $bpair=join("_",@c);
    my $gopair=join("_",@b);
    die("$prot:$go1,$go2\n") if $gopair =~/0008150/;

    ############################
    # Declare database details #
    ############################
    my $host=$db_host;
    my $user="root";
    my $pw="yenapas";
    my $port=$db_port;
   
    ###############################################
    # If this is the first run, connect to the DB #
    ###############################################
    if (scalar keys (%prob)==0) {
	$dsn = "DBI:mysql:database=$database;host=$host;port=$port;"; 
	##Select database 
	$dbh=DBI->connect($dsn,$user,$pw) or die $DBI::errstr;
    }
    if (defined($prob{$gopair})){
	return($prob{$gopair}); 
    }
    if ($go1 eq $go2){  
	$prob{$gopair}=1;
	return($prob{$gopair});
    }
    die("bad gopair $gopair: $bpair\n" ) if defined($prob{$bpair});
    ## If one term is an ancestor of the other prob=1;
    if((defined($offspring{$go1}{$go2}))||
       (defined($offspring{$go2}{$go1}))){
	$prob{$gopair}=1;
	return($prob{$gopair});	
    }
    unless(defined($prob{$gopair})){
	## Define mysql query
	my $query="SELECT P_low from $table where gopair='$gopair'";  
	debug("QQ :$query;\n");
	## Run
	my $result=$dbh->prepare("$query");                        
	## Run query
	$result->execute;                                          
	my $ref = $result->fetchrow_hashref();
	## If a gopair has no prob value, most likely one of its
	## constituent terms does not annotate any proteins DIRECTLY.
	## We can therefore safely ignore the pair by setting P=1.
	$prob{$gopair}=$ref->{'P_low'}||1;
    }
    return($prob{$gopair});
}
############################################################
sub load_ancestors{
    open (ANC,"$genealogy_file")|| die("cannot open $genealogy_file:$!\n");
    while(<ANC>){
	chomp;
	my @terms=split(/\t/);
	my %seen;
	@terms = grep { ! $seen{$_} ++ } @terms;
	my $child=shift(@terms);
	$ancestors{$child}=[@terms];
	map{$offspring{$_}{$child}++}@terms;
    }
    close(ANC);
}
############################################################
sub usage{
    my $us="[options] <gaf file>";
    my $desc="This script will parse a .gaf file and print the number of dissimmilar annotations (PoGo) found for each unique protein.";
	my %opts=
	    (
	     "usage" => $us,
	     "desc" => $desc,
	     "h" => "Print this help and exit",
	     "e" => "Only accept 'good' (identified by experiments) GO terms",
	     "o" => "Subontology, def:P",
	     "G" => "Genealogy file. def:$ENV{HOME}/research/new_moon/data/all_three.genealogy.",
	     "s" => "Species, def: human",
	     "i" => "Check interactome probabilities (def: off).",
	     "m" => "Do not correct for multiple testing",
	     "t" => "Significance threshold, def:0.05",
	     "H" => "pogo host, use 'l' for local tunnel",
	     "T" => "Number of tests for multiple testing correction (if none is given, the real number of tests performed is used.",
	     "v" => "Verbose output"
	    );
    
    print_help_exit(\%opts,0);
}


 
