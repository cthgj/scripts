#!/usr/bin/env perl 
use strict;
use Getopt::Std;
use Math::Round;
use DBI();
our $debug=0;
sub v_say;
sub debug;

#require "MY_SUBS.pl";
my %opts;
getopts('bihs:P:H:',\%opts) || do { print "Invalid option, try '$0 -h' for more information\n"; exit(1); };
usage() if $opts{h};
usage() unless $ARGV[0];
my ($dsn,$dbh);
my $db_port=$opts{P}||3306;
my $db_host=$opts{H}||"10.1.3.30";
my $get_inter=$opts{i}||undef;
my $both=$opts{b}||undef; ## get both high and low probs
## Useful when running via an ssh tunnel
$db_host="127.0.0.1" if $db_host =~ /^l/i;
$db_port=3307 if $db_host eq "127.0.0.1";
my $species=$opts{s}||'human';
$species=&guess_species($species);
my %prob;

## Do we want interaction or annotation probs?
my $mode='annot';
$mode='inter' if $opts{i};


#my $mode=$ARGV[1]||die "Need at least two arguments, a GO term (or list of GOs) and the type of prob we want (inter,annot)\n";
## Read list of GO terms
my @pairs;
my $p;
my @pp;

## If we have been given a list of GOs or GO pairs
if (-f $ARGV[0]) {
    open(G, $ARGV[0])||die "Need a list of GO terms";
    print "GO1    \t   GO2\t        P_low\tP_high\n";
    while (<G>) {
	chomp;
	if (/^(GO:\d+)_(GO:\d+)/ || /^(GO:\d+)\s+(GO:\d+)/) {
	    @pp=($1,$2);
	} 
	else {
	    die "Bad GO file format!\n";
	}
	$p=gopair_probability(@pp,$mode);
	my $kj=join "\t", @pp;
	print "$kj\t$p\n"; 
    }
}
else {
    print "GO1    \t   GO2\t        P_low\tP_high\n";
    my @a=@ARGV;
    s/(GO:\d+)_(GO:\d+)/$1 $2/g foreach @a;
    @pp=map{/(GO:\d+)/g }@a;
    die "No valid GO pairs" unless $#pp>=0;
    $p= &gopair_probability(@pp,$mode);
    my $kj=join "\t", @pp;
    print "$kj\t$p\n"; 
    
}



sub gopair_probability{
    my $go1=shift;
    my $go2=shift;
    my $mode=shift;
    my $low;
    my @b = sort {$b lt $a} ($go1,$go2);
    my @c = sort {$a lt $b} ($go1,$go2);
    my $bpair=join("_",@c);
    my $gopair=join("_",@b);
    #die("$go1,$go2\n") if $gopair =~/0008150/;

    ############################
    # Declare database details #
    ############################
#    my $host="10.1.3.30";
    my $host=$db_host;
    my $port = $db_port;
    my $database="pronto";
    my $user="root";
    my $pw="yenapas";
    my $table="$species";    
    $table="inter_" . $species if $mode ne 'annot';

    ###############################################
    # If this is the first run, connect to the DB #
    ###############################################
    if (scalar keys (%prob)==0) {
	$dsn = "DBI:mysql:database=$database;host=$host;port=$port;"; 
	##Select database 
	$dbh=DBI->connect($dsn,$user,$pw);
    }

    if (defined($prob{$gopair})){
	return($prob{$gopair}); 
    }
    if ($go1 eq $go2){  
	$prob{$gopair}=1;
	return($prob{$gopair});
    }

   ## if this prob has not been queried
    unless(defined($prob{$gopair})){
	## Define mysql query
	my $query="SELECT P_low,P_high from $table where gopair='$gopair'";  
	debug("QQ :$query\n");

	## Run
	my $result=$dbh->prepare("$query");                        
	## Run query
	$result->execute;                                          
	my $ref = $result->fetchrow_hashref();
	$prob{$gopair}=$ref->{'P_low'};
	if ($both) {
	    $prob{$gopair}.="\t" . $ref->{'P_high'};
	}
    }
    unless(defined($prob{$gopair})){
	## Since I am only giving it the under represented ones,
	## any probs missing will be over represented so we can 
	## safely ignore them.
	$prob{$gopair}='unk';
	if ($mode eq 'inter') {
	    $prob{$gopair}=0;
	}
	return($prob{$gopair});
    }

    return($prob{$gopair});
}
############################################################
sub debug{ 
    my $msg=shift;
    chomp($msg);
    print STDERR "\tDEBUG: $msg \n" if $debug; 
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
    elsif($fname=~/([^\/]+).flat/){
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

sub usage{
    print<<EoF;
DESCRIPTION:
  This script will connect to the mysql database at 10.1.3.30 and return the 
  probabilities for each of the GO terms given as input.

  Terms can be given as arguments or in a file. Term pairs, whether in a file or 
  as arguments, should be separated either with as space or an underscore:

    $0 -s human GO:0000377 GO:0006355
    $0 -s human GO:0000377_GO:0006355


USAGE:
 
  $0 [i] -s SPECIES <LIST OF TERMS>

OPTIONS:

  -h: Print this help and exit.
  -H: Database host, default: 10.1.3.30
  -i: Return INTERACTION probabilities.
  -P: Database port, default: 3306.
  -s: Species, one of human, mouse, fly, worm, yeast. Species names can be 
        given in full or as abbreviations or taxids:
         "human" || "hum" || "hs" || "9606"
	 "fly"   || "fb"  || "dm" || "dro" || "dmel" || "7227"
	 "worm"  || "wb"  || "ce" || "ele" || "6239"
	 "mouse" || "mg"  || "mm" || "mus" || "10090"
	 "yeast" || "sgd" || "scc"|| "4932" 
EoF
exit(1);
}
