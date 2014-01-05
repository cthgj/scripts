#!/usr/bin/perl
use strict;
use LWP::Simple;
use Getopt::Std;
require "MY_SUBS.pl";
my %opts;
getopts('w:s:d:l:MbLhc',\%opts) || do { print "Invalid option, try '$0 -h' for more information\n"; exit(1); };
#my @services=("APID=http://cicblade.dep.usal.es/psicquic-ws/webservices/psicquic","BioGrid=http://tyerslab.bio.ed.ac.uk:8080/psicquic-ws/webservices/psicquic","IntAct=http://www.ebi.ac.uk/Tools/webservices/psicquic/intact/webservices/psicquic","DIP=http://imex.mbi.ucla.edu/psicquic-ws/webservices/psicquic","MINT=http://mint.bio.uniroma2.it/mint/psicquic/webservices/psicquic","MatrixDB=http://matrixdb.ibcp.fr:8080/webservices/psicquic","iRefIndex=http://irefindex.uio.no:8080/psicquic-ws/webservices/psicquic","STRING=http://string.uzh.ch/psicquic/webservices/psicquic","Reactome=http://www.ebi.ac.uk/Tools/webservices/psicquic/reactome/webservices/psicquic","InnateDB-IMEx=http://www.ebi.ac.uk/Tools/webservices/psicquic/innatedb/webservices/psicquic","MolCon=http://www.ebi.ac.uk/Tools/webservices/psicquic/molcon/webservices/psicquic","Spike=http://spike.cs.tau.ac.il/psicquic-ws/webservices","TopFind=http://clipserve.clip.ubc.ca/topfind-psicquic-ws/psicquic");
my @services=("APID=http://cicblade.dep.usal.es/psicquic-ws/webservices/psicquic","BioGrid=http://tyerslab.bio.ed.ac.uk:8080/psicquic-ws/webservices/psicquic","IntAct=http://www.ebi.ac.uk/Tools/webservices/psicquic/intact/webservices/psicquic","DIP=http://imex.mbi.ucla.edu/psicquic-ws/webservices/psicquic","MINT=http://mint.bio.uniroma2.it/mint/psicquic/webservices/psicquic","MatrixDB=http://matrixdb.ibcp.fr:8080/webservices/psicquic","Reactome=http://www.ebi.ac.uk/Tools/webservices/psicquic/reactome/webservices/psicquic","InnateDB-IMEx=http://www.ebi.ac.uk/Tools/webservices/psicquic/innatedb/webservices/psicquic","MolCon=http://www.ebi.ac.uk/Tools/webservices/psicquic/molcon/webservices/psicquic","Spike=http://spike.cs.tau.ac.il/psicquic-ws/webservices","TopFind=http://clipserve.clip.ubc.ca/topfind-psicquic-ws/psicquic");
&usage() if $opts{h};
list_dbs() if $opts{L};
$opts{c} && do {
    print STDOUT scalar @services . "\n";
    exit;
};

my $DBindex=$opts{d}||undef;
my $want=$opts{w}||undef;
my $request_lines=$opts{l}||20000;

my $species="all";
my %SPECIES;
defined($opts{s}) && do {

        $species=$opts{s};
	if ($opts{s} eq 'five') {
	    %SPECIES=(
		      '9606' => 1,  # human
		      '7227' => 1,  # fly
		      '6239' => 1,  # worm
		      '10090' => 1, # mouse
		      '4932' => 1,  # yeast
		      '559292' => 1, # yeast
		     );
	}
	else {
	    #####################################
	    # Translate species to default name #
	    #####################################
	    $species=get_species($species, "-s", 1);
    }
};
#print "AA : $species\n"; die();
my $binary=$opts{b}||undef; ## get only binary interactions
my %hq; ## this holds the binary method ids
if($binary){
    my $HA=get_hq_MIs();
    %hq = %{$HA};
}
############################### SUBROUTINES ################################
## List database codes
sub list_dbs{
    for my $i (1..$#services+1){
	my @s=split(/=/,$services[$i-1]);
	print "$i\t:\t $s[0]\n";
    }
    exit;
};

sub getXrefByDbName {
    my ( $xrefs, $dbName ) = @_;
    for my $xref (split (/\|/, $xrefs)){
        my ($db, $id, $txt) = split /[:\(\)]/, $xref;
         if( $db eq $dbName ) {
            return $id;
         }
    }    
    return $xrefs;
}

sub readPsicquic {
    my ( $url, $query, $from, $size ) = @_;
    # Remove the trailing /psicquic at the end of the URL so we can build  a REST URL.
    my $finalUrl = $url =~ s/\/psicquic$//;
 #   my $queryUrl = $url . "/current/search/query/" . $query . "taxidA:$species%20ANDtaxidB:$species?firstResult=" . $from . "&maxResults=" . $size;
    
   # my $queryUrl = $url . "/current/search/query/$query";
    my $queryUrl = $url . "/current/search/query/";

    if ($species eq "five") {
	 $queryUrl.="species:9606%20OR%20species:7227%20OR%20species:6239%20OR%20species:10090%20OR%20species:4932%20OR%20species:559292";
    }
    elsif ($species eq "all") {
	$queryUrl.="$query";
    }
    # elsif ($query ne "*" && $species eq "all") {
    # 	$queryUrl.="$query";
    # }
    # elsif ($query ne "*" && $species eq "five") {
    #     $queryUrl.="$query%20AND%20species:9606%20OR%20species:7227%20OR%20species:6239%20OR%20species:10090%20OR%20species:559292";
    # }
    # elsif ($query eq "*" && $species eq "five") {
    #     $queryUrl.="$query%20AND%2species:9606%20OR%20species:7227%20OR%20species:6239%20OR%20species:10090%20OR%20species:559292";
    # }
    # elsif ($query eq "*" && $species ne "all") {
    #     $queryUrl.="$query?%20AND%2species:$species";
    # }
    else {
	$queryUrl.="species:$species";
    }
    $queryUrl.="?firstResult=" . $from . "&maxResults=" . $size;


#    $query%20AND%20species:$species?firstResult=" . $from . "&maxResults=" . $size;
    print STDERR "QQ : $queryUrl\n";

    my $tries=1;
    my $content = get $queryUrl;
    while($tries<=10){
	if (defined($content)) {
	    $tries=11;
	} 
	else {
	    print STDERR "Could not retrieve $queryUrl, retrying($tries)...\n";
	    $content = get $queryUrl;
	}
	$tries++;
    }
	#print STDERR "BB\n$content\n";
    if($tries==10 && not defined($content)){print STDERR "Could not retrieve $queryUrl after 10 tries\n"; return(0)}
    
    # Now list all interactions
    my @lines = split(/\n/, $content);
    my $LINES= @lines;
    my $count = 0;
    for my $line (@lines) {
      $count++;
      my @flds = split(/\t/, $line); # split tab delimited lines      
      # split fields of a PSIMITAB 2.5 line
      my ($idA, $idB, $altIdA, $altIdB, $aliasA, $aliasB, $detMethod, $author, $pub, $orgA, $orgB, $intType, $sourceDb, $intID, $conf) = @flds;
      
      $detMethod=~s/.+(MI:\d+).+?$/$1/;
      
      ## InnateDB (and maybe others) do not give uniprot as id but as alias
      if($idA !~ /uniprotkb/ && $idA !~ /entrez/){
	  if($altIdA =~ /uniprotkb/){$idA=$altIdA}
	  elsif($aliasA =~ /uniprotkb/){$idA=$aliasA}
	  else{$idA.="|$altIdA|$aliasA";}
      }
      if($idB !~ /uniprotkb/ && $idA !~ /entrez/){
	  if($altIdB =~ /uniprotkb/){$idB=$altIdB}
	  elsif($aliasB =~ /uniprotkb/){$idB=$aliasB}
	  else{$idB.="|$altIdB|$aliasB";}
      }
      ## Ignore interactions between different species
      if ($orgA eq '4932' && ($orgB ne '4932' || $orgB ne '559292')){
	  next;
      }
      elsif ($orgB eq '4932' && ($orgA ne '4932' || $orgA ne '559292')){
	  next;
      }
      else {
	  next unless $orgA eq $orgB;      	  
      }

      ## Choose a species
      # if($species){ 
      # 	  $orgA=~/taxid:(\d+)/||next;
      # 	  my $aa=$1;
      # 	  $orgB=~/taxid:(\d+)/||next;
      # 	  my $bb=$1;
      # 	  unless ($species eq "all"){
      # 	      if ($species eq 'five') {
      # 		next unless defined($SPECIES{$aa});
      # 		next unless defined($SPECIES{$bb});
      # 	      }
      # 	      else {
      # 		  next unless $aa eq $species;
      # 		  next unless $bb eq $species;
      # 	      }
      # 	  }
      # }

    ## Do we only want binary interactions?
      if($binary){
	  next unless defined($hq{$detMethod});
      }
      print getXrefByDbName($idA, "uniprotkb") . "\t" . getXrefByDbName($idB, "uniprotkb") . "\t$detMethod\t$orgA\t$orgB\t$intType\t$pub\t$sourceDb\t$intID\t$conf\n";
      
    }
    return ($LINES,$count);
}

################################################################################
# Registry URL that lists all ACTIVE PSICQUIC services
# Documentation: http://code.google.com/p/psicquic/wiki/Registry
my $registryUrl = 'http://www.ebi.ac.uk/Tools/webservices/psicquic/registry/registry?action=ACTIVE&format=txt';


## We can either get ALL interactions or a specific set
my @wanted;
#############################################
# Collect the proteins we are interested in #
#############################################
my @names;
if ($want) {
    my @wanted=split(/[,\s]/,$want);
    $"=" x ";
    print "W: @wanted\n";
    foreach my $ar (@wanted) {
        if (-e $ar){
	    open(A, $ar);
	    while(<A>){
		chomp;
		push @names,$_;
	    }
	}
	else{push @names,$ar};
    }
}
if($want){
    @wanted=split(/\s+/,$want);
}
## Get ALL interactions
else{
    $wanted[0]="*";
}


my $totalCount = 0;

for my $ii (0..$#services){
    ## If we have asked for a specific DB service
    if(defined($DBindex)){                          
	next unless $ii==$DBindex-1;
    }

    my @flds = split(/=/, $services[$ii]);
    my ($serviceName, $serviceUrl) = @flds;
    print STDERR "$serviceName  --->  $serviceUrl\n";
    
    my $current = 0;
    my ($lineCount, $count);
    miq:foreach my $miql (@wanted){
    	do {	    
    	    ## Download a MITAB page ($request_lines lines max per read)  
    	    ## max allowed size (at least for apid) is 2147483647
            ($lineCount, $count) = readPsicquic $serviceUrl, $miql, $current, $request_lines;
            $current = $current + $count;
            ## Get just the first N interactions
            if ($opts{M}){last miq if $current>=1};## just for debugging
    	} while( $lineCount == $request_lines ); 
    }
    #print STDERR "Total: " . $current . " interaction(s).\n\n";
    $totalCount = $totalCount + $current;
}
$DBindex ? 
    print STDERR "Retreived " . $totalCount . " interaction(s) from " . $services[$DBindex-1] . "\n\n" :
    print STDERR "Retreived " . $totalCount . " interaction(s) across " . @services . " service(s)\n\n";


sub usage{
    my $us="[options] <psicquic file>";
    my $desc="This script will query multiple psicquic-enabled databases and retrieve all interactions for the given species. Currently implemented databases are: APID, BioGrid, IntAct, InnateDB, DIP, MINT, MatrixDB, iRefIndex, BIND, STRING, Reactome, InnateDB-IMEx, MolCon, I2D-IMEx, I2D, Spike, TopFind\n";
    my $ser_num=$#services+1;
    my %opts=(
	      "usage" => $us,
	      "desc" => $desc,
	'w' => "Protein IDs we want. If this is set, only interactions involving these proteins will be retrieved.",
	's' => "Species. If the value is 'five', all interactions for human,fly,yeast,worm and mouse will be downloaded.",
	'd' => "Database. Takes a value between 1 and $ser_num. Use -L to list available databases. 0 will query ALL databases.",
	'b' => "Retrieve only binary interactions",
	'c' => "Count the number of available databases and exit",
	'L' => "List database codes",
        'l' => "Number of lines to fetch from the database per request. This parameter does NOT affect the total number of interactions returned, just how many are downloaded by each http request made to the database. (def: 20000)",

	'M' => "Only download the 1st 20000 (or whatever passed as -l) interactions"
	
);
    print_help_exit(\%opts,0);
}
