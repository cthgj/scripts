#! /usr/bin/perl -w

use strict;
use Getopt::Std;
use Term::ANSIColor; 
use Text::Table;
our $verbose;
our $debug;
sub v_say;
sub debug;
require "MY_SUBS.pl";

#Takes a species name and a results directory. For a given list of modes, return all cands from species that were also cands in another species

#die "Need at least 2 arguments, the moonGO outfile (1st arg) and an outfile of metaphors.pl (2nd arg)\n" unless $ARGV[1];

my (%synonyms,%cands, %opts, %orths, %mapped, %has_orth, %moon);
getopts('hvcgm:M:r:s:o:',\%opts) || do { print "Try '" .  this_script_name() . " -h' for more information\n"; exit(1); };

&usage() if $opts{h};
my $count=$opts{c}||undef;
my $no_map=undef;
my $species=$opts{s}||'human';
my $synfile=$opts{m}||undef;
my $results_dir=$opts{r}||"./results";
my $MODES=$opts{M}||"i,B";
my $o_MODES=$opts{o}||"i,B";
my $orths_file=$ARGV[0];
$species=guess_species($species);

$no_map=1 unless $synfile;
$synonyms{LOADED}=0;
my @SPECIES=qw(HUMAN MOUSE DROME CAEEL YEAST);
my %sp_map=(
	    'human' => 'HUMAN',
	    'mouse' => 'MOUSE',
	    'fly'   => 'DROME',
	    'worm'  => 'CAEEL',
	    'yeast' => 'YEAST',
	    'HUMAN' => 'human',
	    'MOUSE' => 'mouse',
	    'DROME' => 'fly',
	    'CAEEL' => 'worm',
	    'YEAST' => 'yeast',
);
######################################
# Get the modes we are interested in #
######################################
my @mm=split(/,/, $MODES);
my %modes;
map{$modes{$_}++}@mm;
my @omm=split(/,/, $o_MODES);
my %o_modes;
map{$o_modes{$_}++}@omm;

## Initialize counts to 0	
my %pair_modes;
foreach my $c_mod (@mm) {
    foreach my $o_mod (@omm) {
	my $mpair=join("<=>",$c_mod,$o_mod);
	$pair_modes{$mpair}=0;
    }
}


##########################################################
# Collect moonGO out files for $species. 		 #
# Read ONLY files, and only moonGO out files of $species #
##########################################################
opendir(RES, $results_dir) || die "Could not open results directory $results_dir for reading: $!\n";
my @species_res=grep { (/^(.+?)\.(\w)\.moon$/) && -f "$results_dir/$_" && defined($modes{$2}) && $1 eq $species} readdir(RES);
close(RES);
#################################################################
# The keys of %{$moon{$species}} are the cands for that species #
#################################################################
%{$moon{$species}}=%{get_cands_from_moon_file(@species_res)}; 

##########################################################
# Collect moonGO out files for ALL other species. 	 #
# Read ONLY moonGO out files                             #
##########################################################
opendir(RES, $results_dir) || die "Could not open results directory $results_dir for reading: $!\n";
my @other_res=grep { (/^(.+?)\.(\w)\.moon$/) && -f "$results_dir/$_" && defined($o_modes{$2}) && $1 ne $species} readdir(RES);
close(RES);
for my $ss (@SPECIES) {
    my $s=$sp_map{$ss};
    next if $s eq $species;
    %{$moon{$s}}=%{get_cands_from_moon_file(@other_res)}; 
}

###########################################
# Collect the orthologs of each candidate #
###########################################
my $orth_fh=check_file($orths_file, 'r');
while (<$orth_fh>) {
    chomp;
    my @a=split(/\s+/);
    ## Get the cand name
    my $cand=get_name(shift(@a),1);
    ## Collect the orthologs (ACs) for this cand (ID)
    map{$orths{$cand}{get_name($_,1)}++;}@a;
}


##############################################################
# Now, for each candidate from $species, check if its orths  #
# are also cands in another species 		             #
##############################################################
my (%tot,%num_found);
foreach my $cand (keys(%{$moon{$species}})){
    my $i=0;
    my $num=0;
    my %k;
    my %cols;
    foreach my $o (keys(%{$orths{$cand}})) {
	my (%terms1,%terms2)=();
	$o=~/^.+_(.+)/ || next; ## skip if we don't have a UniProt ID
	next unless defined($sp_map{$1});
	my $o_species=$sp_map{$1};
	## Skip unless this orth is a candidate in its species
	next unless defined($moon{$o_species}{$o});
	#if ($count) {
	$num_found{$o_species}{tot}++;
	#print "$cand (@ll), $o \$num_found{$o_species}{tot} $num_found{$o_species}{tot} (@lo)\n";
	$tot{$cand}++;
	#}
	my @ll=keys(%{$moon{$species}{$cand}}); 
	my @lo=keys(%{$moon{$o_species}{$o}}); 
	
	foreach my $c_mod (@ll) {
	    foreach my $o_mod (@lo) {
		my $mpair=join("<=>",$c_mod,$o_mod);
		$num_found{$o_species}{pair}{$mpair}++;
	    }
	}

	map{$num_found{$o_species}{met}{$_}++}keys(%{$moon{$species}{$cand}});
	
	unless ($count){
	    if ($i==0){
		print "$cand\t" ;
		if ($opts{g}) {
		    my @met=keys(%{$moon{$species}{$cand}});
		    foreach my $m (@met) {
			my @gos=map{split(/_/)}keys(%{$moon{$species}{$cand}{$m}});
			%terms1=terms2gos(0,@gos);
			my $kk=0;
			map{
			    if ($kk==0) {
				print "$m, ";
				$kk=1;
			    } else {
				print "\t\t   ";			
			    }
			    my @a=split(/_/);
			    $cols{$terms1{$a[0]}}="bold blue";
			    $cols{$terms1{$a[1]}}="bold yellow";
			    print "\"" . color("bold blue")  . $terms1{$a[0]}.color("reset") . "\" and \"";
			    print color("bold yellow") .  $terms1{$a[1]}.color("reset") . "\" \n ";

			}keys(%{$moon{$species}{$cand}{$m}});
		    }
		}
		$i=1;
	    }
	    print "\n" . &get_name($o,1) ;
	
	    ####################################################
	    # Check what methods this cand was found by	   #
	    ####################################################
	    if ($opts{g}) {
		my @met=keys(%{$moon{$o_species}{$o}});
		foreach my $m (@met) {
		    my @gos=map{split(/_/)}keys(%{$moon{$o_species}{$o}{$m}});
		    %terms2=terms2gos(0,@gos);
		    my $kk=0;
		    map{
			if ($kk==0) {
			    $i==0 ? 
				print "\t$m, aa\"" :
			    print "\t\t$m, bb\"";
			    $kk=1;
			} else {
			    print "\t\t   \"";			
			}
			my @a=split(/_/); 

			defined($terms1{$a[0]}) ? 
			    #print color("bold blue").$terms2{$a[0]}.color("reset") :
			    print color($cols{$terms1{$a[0]}}).$terms2{$a[0]}.color("reset") :
				print $terms2{$a[0]};
			print "\" and \"";
			defined($terms1{$a[1]}) ? 
			    print color($cols{$terms1{$a[1]}}).$terms2{$a[1]}.color("reset") : 
				print $terms2{$a[1]};
			print "\"\n" ;

			# $terms1{$a[0]}=~s/$terms1{$a[0]}/color("bold blue").$terms1{$a[0]}.color("reset")/ if defined($terms1{$a[0]});
			# print ": \"$terms2{$a[0]}\" and \"$terms2{$a[1]}\" | "
		    }keys(%{$moon{$o_species}{$o}{$m}});
		}
	    }
	}
	$num++;
	$k{$o_species}++;
    
	#print "\n";
	unless ($i==0 ||$opts{g}) {
	    print "\t$num\t" . keys(%k) . "\n" unless $i==0;
	}
    }  
	print "\n-----------------------------------------\n" unless $i==0;

}
my $tot=0;



foreach my $s (keys %num_found) {
   print "$s\t$num_found{$s}{tot}\t";
   map{
       $num_found{$s}{met}{$_}//=0;
       print "\t$_: $num_found{$s}{met}{$_}"
   }keys(%pair_modes);#{$num_found{$s}{met}});
   map{
       $num_found{$s}{pair}{$_}//=0;
       print "\t$_: $num_found{$s}{pair}{$_}"
   }keys(%pair_modes);#{$num_found{$s}{pair}});
   print "\n";
}

#map{print "$_\t $num_found{$_}{tot}, $num_found{$_}{\n"; }keys %num_found;
print "Total\t" . scalar keys(%tot) . "\n"; 
############################################################
sub get_cands_from_moon_file{
    my %moon_cands;
    foreach my $file (@_) {
	$file =~ /(.+)\.(.+?)\.moon/;
	my $sp=$sp_map{$1};
	my $mode = $2;
	my $ff=check_file("$results_dir/$file", 'r');
	while (<$ff>) {
	    next if /^\s/;
	    next if /^\#/;
	    chomp;
	    /(GO:\d+)\s.+?(GO:\d+)/;
	    my $go=join("_", sort($1,$2));
	    /^(.+?)\s/;
	    my $a=$1;
	    $moon_cands{$1}{$mode}{$go}++;
	}
    }
return(\%moon_cands);

}

############################################################
sub get_name{
    my $name=$_[0];
    return ($name) if $no_map ;
    my $old_name=$name;
    my $reverse=$_[1]||undef;
    if ($synonyms{LOADED}==0){
	if(-e $synfile){
	    open(S,$synfile);
	    while(<S>){
		chomp;
		if(/^([^\t]+?)\t([^\t]+)$/){
		    $synonyms{NAME}{$2}=$1; 
		    $synonyms{ACC}{$2}=$2; 
		    $synonyms{NAME}{$1}=$1; 
		    $synonyms{ACC}{$1}=$2; 
		    # print "\$synonyms{NAME}{$2} : $synonyms{NAME}{$2}\n\$synonyms{ACC}{$2} : $synonyms{ACC}{$2}";
		}
		else{
		    die("Bad format synonyms file : $synfile\n");
		}
	    }
	    close(S);
	}
	else{die("Could not open $synfile\n" );	}	
 	$synonyms{LOADED}=1;
    }
    unless (defined($synonyms{NAME}{$old_name})){
	$synonyms{NAME}{$old_name}=$name;
    }
    unless (defined($synonyms{ACC}{$old_name})){
	$synonyms{ACC}{$old_name}=$name
    }
    $reverse ? 
	($name=$synonyms{NAME}{$old_name}) : ($name=$synonyms{ACC}{$old_name});
  return $name;
}


sub usage{
    my $us="[options] -m <id2ac map file> <ortholog list>";
    my $desc="Takes a species name and a results directory. For a given list of modes, return all candidates that were also cands in another species";
    my $ss=<<EOF;
one of the following:
                 "human" || "hum" || "hs" || "9606"
                 "fly"   || "fb"  || "dm" || "dro" || "dmel" || "7227"
	         "worm"  || "wb"  || "ce" || "ele" || "6239"
	         "mouse" || "mg"  || "mm" || "mus" || "10090"
	         "yeast" || "sgd" || "scc"|| "4932"
EOF

    my %opts=(
	      "usage" => $us,
	      "desc" => $desc,
	      "s" => "Species, $ss",
	      "m" => "Synonyms map file. One protein name per line: \"UNI_ID\tUNI_AC\" ",
	      "r" => "Results directory (def:./results), will collect <species>.moon files from here. ",
	      "M" => "Comma separated list of MoonGO.pl prediction modes we are interested in (def: i,B), only candidates from these modes will be counted.", 
	      "o" => "Comma separated list of MoonGO.pl prediction modes we are interested in for the orthologs (def: i,B), only candidates whose orthologs are also candidates in one of these modes will be counted.",
	      "g" => "Foreach candidate, print the GO pairs that have made it one. Not very human readable, use parser parse_moon_orths.pl"
	     );
    print_help_exit(\%opts,1);
    
}
