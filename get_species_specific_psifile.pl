#!/usr/bin/perl -w
use strict;
use Getopt::Std;
use feature qw(switch say);
our $verbose;
our $debug;
our $force;
sub v_say;
require "MY_SUBS.pl";
my %opts;
getopts('hvSbp:s:',\%opts);

usage() unless $ARGV[0];
usage() if$opts{h};

my $standard=$opts{S}||undef;
$verbose=$opts{v}||undef;
my $species=$opts{s}||undef;
my $binary=$opts{b}||undef;
my $prefix=$opts{p}||undef;
my $taxid;
my %TAXIDS;
my ($fhH, $fhW, $fhF, $fhY, $fhM);
my @FILEHANDLES=($fhH, $fhW, $fhF, $fhY, $fhM);
my @ss=qw(human worm fly yeast mouse);
#print "FF: @FILEHANDLES\n";
if ($standard) {
    for (my $i=0; $i<=$#ss; $i++) {
	$species=get_species($ss[$i],"-s");
	$taxid=get_species($species,"-s",1);
	#print "tt : \$TAXIDS{$taxid}=$FILEHANDLES[$i];$taxid\n";
	$TAXIDS{$taxid}=$species;
	open($fhH,  ">", "human.$prefix.psi");
	open($fhM,  ">", "mouse.$prefix.psi");
	open($fhY,  ">", "yeast.$prefix.psi");
	open($fhW,  ">", "worm.$prefix.psi");
	open($fhF,  ">", "fly.$prefix.psi");
    }
my @a=keys(%TAXIDS) ;
}
else {
    if($species){
	$species=get_species($species,"-s");
	$taxid=get_species($species,"-s",1);
    }
    else{
	$species=guess_species($ARGV[0]);
	$taxid=get_species($species,"-s");
    }
    if ($species eq '-1'){die("Please specify a species (-s)\n");}
}

my $fh=check_file($ARGV[0],'r' )||die("Could not open $ARGV[0] : $!\n");

my %hq;
if ($binary) {
    %hq=%{get_hq_MIs("b")}
}
while (<$fh>) {
    next if /^Total:/;
    next unless /\w/;
    v_say("Psi file, line : [$.]\r", " ", 1) if $. % 10000==0;
    my ($idA, $idB, $detMethod, $orgA, $intType, $pub, $sourceDb, $conf) = split(/\t/);
    $orgA=~s/taxid:(\d+).*/$1/;
    
    ## This should already have been done by 
    ## get_psicquic_interactome.pl. But, just in case...
    $detMethod=~s/.+(MI:\d+).+?$/$1/;

    if($binary){
	next unless defined($hq{$detMethod});
    }

    if ($standard) {
	next unless defined($TAXIDS{$orgA});
	my $fh=get_filehandle($orgA);
	print $fh "$_";
    }
    else {
	next unless $orgA eq $taxid;
	print;
	
    }
}

sub get_filehandle{
    my $txid=shift();
    
    if ($txid eq "9606") {
	return(*$fhH)
    }
    elsif ($txid eq "7227") {
	return(*$fhF)
    }
    elsif ($txid eq "6239") {
	return(*$fhW)
    }
    elsif ($txid eq "10090") {
	return(*$fhM)
    }
    elsif ($txid eq "4932") {
	return(*$fhY)
    }
    else {
	die("Unkown taxid $txid\n");
    }

   
}

sub usage{
      print STDERR <<EndOfHelp;

  USAGE: get_species_specific_psifile.pl [OPTIONS] <psi_file>
  -b : Retrieve only binary interactions.
  -s : Species.
  -S : Standard use. This will create a species specific file for human, worm, fly, yeast and mouse. All other species will be ignored.
  -p : Prefix, for Standard mode. Output files will be named: <SPECIES>.<PREFIX>.psi.
EndOfHelp

exit;

}
