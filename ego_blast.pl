#!/usr/bin/perl -w

#  This script takes the output of SECISearch in fasta format, creates the clusters_found 
#  and cluster_morethan2 files, tbl files of all the sequences of each cluster in the morethan2
#  file and then blasts the first two sequences of a different species in a cluster against each
#  other. The individual subprocesses can be run as parse_egosecis.sh, cluster.pl and ego_blast.pl.o

#                         USAGE : ego_blast.pl NAME.secis

#use strict;
use Getopt::Std;


#################################################################
##     Get cmd-line options and declare some variables         ##
#################################################################

my %opts;
getopts('SbpscbnChw',\%opts) || do { print "Invalid option, try 'ego_blast.pl -h' for more information\n"; exit(1); };

my $p       = $opts{p} || 0;  ###  do not parse_egosecis
my $s       = $opts{s} || 0;  ###  do not sels
my $c       = $opts{c} || 0;  ###  do not cluster
my $b       = $opts{b} || 0;  ###  do not blast 
my $noseq   = $opts{n} || 0;  ###  do not get full gene sequence
my $noseqcl = $opts{C} || 0;  ###  do not get @clusters foreach seq
my $w       = $opts{w} || 0;  ###  blast the entire cluster against itself not jst seqs w pred. SECIS
my $blsels  = $opts{S} || 0;  ###  blast sels as well 
&usage() unless ($ARGV[0]);
&usage() if $opts{h} ;



my $egoseq = '/seq/databases/EGO/ego8/ego8_112404.seq';
my %clusters = ();
my %ego  = ();
my %sels = ();
my @files;
my $pat;
my $secisout = $ARGV[0];   # get secisout filename
my @clusnames;
my %names = ();
my %seqcl =(); ### $seqcl{seq name} = @{clusters}
my @snames;
my @sclus;
my %dfspcsec;
my $d;
my %uniqclusfound = ();

my %tbl = ();
my %db = ();
	

my @keysa;
my @keysb;

chomp $secisout;

	print STDERR "\$blsels : $blsels\n";

### Get full gene sequence
&ego() unless($noseq);

### Get clusters for each gene
&seqcl() unless($noseqcl);

### Create clusters_found and cluster_morethan2 files
&parse_egosecis($secisout) unless $p;

### Get names of clusters of known Selenoproteins
&sels() unless $s;

### Create tbl files of secises from each cluster
&cluster(); 

### Identify clusters with SECISes predicted in >1 species
&indy_secis();

### Blast 2 gene sequences from each cluster against each other

&blast() unless $b;

#########################################################################################################   
###                                            Functions                                              ###
#########################################################################################################

sub ego   ### Get full gene sequence
{
    {
	open(SEQ,"$egoseq") || die "0 cannot open seq: $!";
	while (<SEQ>)
	{
	    /^(.*?)\s(.*?)$/;
	    $ego{$1} = $2;     ### $ego{seq name} = sequence
	}
	close (SEQ);
	print STDERR "*** get full seq done\n";
    }
}

###############################################################################
###############################################################################

sub seqcl ### Get @clusters
{
open(SEQCL,"/home/ug/cchapple/research/seq/4ego/sequences_and_clusters") ||  die "cannot open sequences_and_clusters: $!";
    while(<SEQCL>)
    {
	/^(.*?)\s+(.*)/og;
	my @f = split /\s+/o, $2 ;
	$seqcl{$1} = [@f];      ### $seqcl{seq name} = @clusters it belongs to
    }
    close(SEQCL);
    print STDERR "*** SEQCL done\n";

}

###############################################################################
###############################################################################

sub parse_egosecis
{
    
    $_[0] =~ /(pat.*?)\.secis/g;
    my $patname = $1;
    if ($patname =~ 'pat_standard_I_II') { $pat = 'std'}
    elsif ($patname =~ 'pat_non_standard_I_II'){$pat = 'non_std'}
    elsif ($patname =~ 'pat_twilight_I_II'){$pat = 'twil'}
    else {$pat = $patname}
    
### Get secis containing sequences' names
    open(SECIS, "$_[0]") || die "cannot open $_[0] : $!";
    while(<SECIS>)
    {
	next unless /^>/;
	/>(.*?)\s+/og;
	$names{$1} = '1';       ### $names{seq name} = defined if seq has predicted SECIS 
    }
    @snames = keys(%names);
    close(SECIS);
   
    foreach my $name(@snames){
	push @clusters, $seqcl{$name};   ### Push @clusters, @listofclusters foreach $secisname
	
    }
    foreach my $clu (@clusters){   ### Cycle through @clusters to find repeated ones
	foreach my    $cluster    (@$clu){
	    if (defined($clusters{$cluster})){
		$clusters{$cluster}++;
	    }
	    else{
		$clusters{$cluster} = 1;
		
	    }
	}
    }
    
    my @keys = keys %clusters;
    
    foreach my $key(@keys){
	if ($clusters{$key} > 1){
	    push @mclus, $key;
	}
    }
    
    open (CLUS2,">cluster_morethan2_$pat")|| die "cannot open cluster_morethan2_$pat :$1" ;
    foreach my $cl(@mclus){ print CLUS2 "$clusters{$cl}\t$cl\n";}
    close(CLUS2);
    
    foreach my $sl(@clusters){
	foreach my $o (@$sl){
	    $uniqclusfound{$o}++;
	}
    }
    open(CLUS, ">clusters_found_$pat")|| die "cannot open clusters_found_$pat : $!";
    my @uniqclus = keys(%uniqclusfound);
    map{print CLUS "$_\n"} @uniqclus;
    
    close(CLUS);
    print STDERR "*** parse_egosecis done\n";  
}

###############################################################################
### Separate known sel clusters (compare with results of previously done blast)
###############################################################################

sub sels
{
    open (HITS, "/home/ug/cchapple/research/seq/4ego/hits") || die "cannot open hits: $!";
    
    my @sel = <HITS>;
    close(HITS);
        
    foreach my $t (@sel)
    {
	chomp $t;
	map {$sels{$_} = '1'} @{$seqcl{$t}}   ### $sels{cluster} = defined if sel cluster
    }
    
    print STDERR "*** sel done\n";
} 
###############################################################################
###############################################################################

sub cluster
{
    open (CLUSTERS, "cluster_morethan2_$pat")|| die "cannot open cluster_morethan2_$pat: $!";
    
    my @clusses= <CLUSTERS>;                  ### get all cluster > 2 names

    foreach my $a (@clusses)
    {
	my @e = split /\s+/o, $a;            ### @e = list of INDIVIDUAL clusters
	map { push @clusnames, $1 if $_ =~ /(TOG\d{6})/go;} @e;   ### @clusnames = clusters_found > 2
    }
    
    
	mkdir("TOG_files", 0777) unless (-d "TOG_files");
	mkdir("TOG_files/sel_clusters") unless (-d "TOG_files/sel_clusters");
		
	system("FastaToTbl $secisout > tmptbl");
	
	open (CL, "/home/ug/cchapple/research/seq/4ego/clusters")|| die "cannot open /home/ug/cchapple/research/seq/4ego/clusters: $!";
	while(<CL>)
	{
	    /^(TOG.*?)\t(.*?)$/o;
	    my @d = split /\s+/o, $2;
	    $clusters{$1} = [@d];           ### $clusters{cluster} = member names foreach cluster  
	}

    unless($c)
    { 	
	open (FILE, "tmptbl")|| die "cannot open tmptbl: $!";
		
	foreach my $c (@clusnames)
	{   	
	    if (defined  $sels{$c})
	    { 	
		print STDERR "#";
		unlink ("TOG_files/sel_clusters/$c\_$pat\.tbl") if (-r "TOG_files/sel_clusters/$c\_$pat\.tbl");
		open (CLUS, ">>TOG_files/sel_clusters/$c\_$pat\.tbl")|| die "cannot open TOG_files/sel_clusters/$c\.tbl: $!";
		map {print CLUS "$_ $ego{$_}\n" if $names{$_};} @{$clusters{$c}}; ### Print only those sequences with predicted SECIS
		close (CLUS);
		next;
	    }
	    else 
	    {	
		print STDERR "*";
		unlink ("TOG_files/$c\_$pat\.tbl") if (-r "TOG_files/$c\_$pat\.tbl");
		open (CLUS, ">>TOG_files/$c\_$pat\.tbl")|| die "cannot open TOG_files/$c\.tbl: $!";
		map {print CLUS "$_ $ego{$_}\n" if $names{$_};} @{$clusters{$c}}; ### Print only those sequences with predicted SECIS
		close (CLUS);
		
	    }
	} 
	
	close(FILE);
	unlink ("tmptbl");
	print STDERR "\n*** cluster done ***\n";
    }
    
}
###############################################################################
###          Separate clusters w SECIS predicted in >1 species.             ###
###############################################################################

sub indy_secis{
    CLUS:    foreach my $clus (@clusnames)
    {
	next if defined  $sels{$clus};
	open(FILE,"TOG_files/$clus\_$pat\.tbl")|| die "cannot open TOG_files/$clus\_$pat\.tbl: $!";
	my @lines = <FILE>;
	$lines[0] =~ /^(.*?)_/;
	my $prevspec = $1;
	my @nms;
	close(FILE);
	
	open(FILE,"TOG_files/$clus\_$pat\.tbl")|| die "cannot open TOGfile: $!";
	
#	print STDOUT "\$prevspec is now $prevspec\n";
	
	while(<FILE>){
	next if $_ =~ /$prevspec/;

	$dfspcsec{$clus}++;    ### DiFferent SPeCies SECis
	push @nms, $clus;
	next CLUS;
    }
	

	close(FILE);
	
	

    }
      open(MCLUS, ">msecis");  ### Clusters w SECIS predic. in different organisms
      open(SCLUS, ">ssecis");  ### Clusters w SECIS predic. in one organism
      foreach my $clus (@clusnames)
      {
	  if ($dfspcsec{$clus}){
	      print MCLUS "$clus\n";
	  }
	  else{
	      print SCLUS "$clus\n"
	  }
      }
      close(MCLUS);
      close(SCLUS);
      print STDERR "*** indy_secis done ***\n";

  }
###############################################################################
###############################################################################

sub blast
{
    unlink("query\.tbl");  ### Delete files possibly present from previous instance 
    unlink("dbase\.tbl"); 

    mkdir("blast", 0777) unless (-d "blast");
    mkdir("blast/single_spec", 0777) unless (-d "blast/single_spec");
    foreach my $clus (@clusnames)
    {
	
	
	if (defined $sels{$clus}) ### Blast only non-sl clusters, unless $blsels is set.
	{
	    if ($blsels == 0){next;}
	    else {open (FILE, "TOG_files/sel_clusters/$clus\_$pat\.tbl") || die "1 cannot open TOG_files/sel_clusters/$clus\_$pat\.tbl : $!";}
	}

unless($blsels){
	    next if defined  $sels{$clus}; print "**********\n";
	};   ### Blast only non-sl clusters, unless $blsels is set.
	%tbl = ();
	%db = ();
	
	print STDERR "################################################\n###\t\$clus : $clus\t###\n################################################\n";
	open (TBL, ">>query\.tbl")|| die "2 cannot open query: $!";
	open(DB, ">>dbase\.tbl");
	my @secises;
	while(<FILE>)  
	{
	    $_ =~ /^(.*?)\s+/og;
	    push @secises, $1 if $names{$1};
	}
	close(FILE);	    
	print STDERR "\@secises : @secises\n";
#	open (FILE, "TOG_files/$clus\_$pat\.tbl") || die "1 cannot open TOG_files/$clus\_$pat\.tbl : $!";
#		my $num = scalar(@{$clusters{$clus}});
	$d = scalar(@secises) / 2;
	print STDERR "****************\n\$d :$d\n******************************\n";

	$p = 0;
	$w = 1 if $d < 1;
	
	@keysa = ();
	@keysb = ();
	
	print STDERR "START : keysa:  @keysa\nkeysb : @keysb\n";
	
	foreach my $secis(@secises)
	{
	    if ($w)
	    {
		&blast_whole_clus($secis);
	    }
	    else
	    {
		&blast_secis_clus($secis);
	    }
	}
	close (TBL);
	close(DB);
	#%tbl = ();
	#%db = ();
	
	system ("TblToFasta query\.tbl > query\.fa") && die "shit1 : $!\n";
	system ("TblToFasta dbase\.tbl > dbase\.fa") && die "shit2 : $!\n";
	
	unlink("query\.tbl"); 
	unlink("dbase\.tbl"); 
	system("pressdb dbase\.fa");
	
	if (defined $dfspcsec{$clus}){     ### If this cluster has SECIS predicted in >1 species
	    print STDOUT "***defined***\n";
	    
	    system ("wu-tblastx dbase\.fa query\.fa > blast/$clus\.out"); 
	}
	else {
	    print STDERR "###not defined###\n";
	    
	    system ("wu-tblastx dbase\.fa query\.fa  -hspmax 3 -S 600 > blast/single_spec/$clus\.out"); 
	}
	
	unlink("dbase\.fa\.ntb");
	unlink("dbase\.fa\.csq");
	unlink("dbase\.fa\.nhd");
	unlink("dbase\.fa");
	unlink("query\.fa");
	print STDERR "****** One Down!!!!! ******\n"; 
    }
    print STDERR "*** blast done ***\n";
}
#####################################################
#####################################################
sub blast_secis_clus
{
    print STDERR "\n***secis blast***\n";
    if ($p < $d)
    {
	print TBL "$_[0]\t$ego{$_[0]}\n" unless defined $tbl{$_[0]};
	$tbl{$_[0]}++;
	$p++;
	 @keysa = keys(%tbl);
	print STDERR "keysa : @keysa\nkeysb : @keysb\n";
    }
    else
    {
	print DB "$_[0]\t$ego{$_[0]}\n" unless defined $db{$_[0]};
	$db{$_[0]}++;
	@keysb = keys(%db);
	print STDERR "keysa : @keysa\nkeysb : @keysb\n";
    }
}

#####################################################
#####################################################
sub blast_whole_clus
{
    print STDERR "\n***whole blast***\n";
    if($p < $d)
    {
	print TBL "$_[0]\t$ego{$_[0]}\n"; 
	$p++;
    }
    
    else
    {
	print DB "$_[0]\t$ego{$_[0]}\n";
	$p++;
    }
}

###############################################################################
###############################################################################

sub usage
{
    print STDERR "\n*** ego_blast.pl : This script takes the output of SECISearch in fasta format and the EGO.seq file from the EGO database, creates the clusters_found and cluster_morethan2 files, tbl files of all the sequences of each cluster in the morethan2 file and then blasts each cluster with a SECIS predited in >1 different species against itself.\n\n############################################\nUSAGE : ego_blast.pl EGO.seq secisout.secis\n############################################\n\n"; 

    print STDERR <<"EOF";

OPTIONS : 
    -b  :  Do not run the blast routine
    -c  :  Do not run clusters
    -C  :  Do not get \@clusters foreach seq
    -m  :  Blast all SECIS clusters, not only those with SECIS predicted in >1 organism.
    -n  :  Do not get \%ego
    -p  :  Do not parse secis
    -s  :  Do not do sels routine
    -w  :  Blast whole clusters against themselves instead of only those sequences with a predicted SECIS in each cluster.
EOF
exit();
}


###############################################
#  Fossil record...                           #
#                                             #         
###############################################

####### Old sels #######
#qw(TOG253979-SelT.tbl TOG262411-TR3sel.tbl TOG286938-glut_perox_3_put.tbl TOG256592-presenilin.tbl TOG264519-CG9313.tbl TOG288351-sperm_assoc_antgn.tbl TOG257832-plants.tbl TOG264955-SelK.tbl TOG292074-SHC_transf.tbl TOG258585-glut_perox.tbl TOG265010-AD-015.tbl TOG298761-ras_rel_gtp_bind.tbl TOG260439-lys_peps-ins_proase.tbl  TOG265639-SelM.tbl TOG299534-rhotekin.tbl TOG262345-SelW.tbl TOG269344-type2_iodothyr_diod.tbl TOG304731-glut_perox_rel.tbl TOG262355-SelX.fa TOG272723-SelP.tbl TOG262355-SelX.tbl TOG279450-KIAA0869.tbl		 );
########################    
  
################## Old blast ##########################

#	my $query = 'la';
	
#	next unless (-r "TOG_files/$clus\_$pat\.tbl");
#	print STDERR "\$clus is $clus\n";
#	open (FILE, "TOG_files/$clus\_$pat\.tbl") || die "1 cannot open $clus\.$pat\.tbl : $!";
#	my @lines = <FILE>;
#	$lines[0] =~ /^(.*?)-(.*?)(:|\s)/;
#	my $sq = $2;
#	my $species = $1;
	
#	print "\$species is $species\n";
#	print STDERR "\$name is $name\n";    
#	my $n = 1;
#	foreach my $line (@lines)
#	{
	    
#	    next if $line =~ /$species/;
	    
#	    $line =~ /^(.*?)-(.*?)\W/;
#	    $query = $2;
	    
#	}
#	if ($query eq 'la')
#	{
#	    $lines[1] =~ /^(.*?)-(.*?)\W/;
#	    $query = $2;	
#	}

#	open (TBL, ">$name\.$pat\.tbl")|| die "2 cannot open name: $!";
#	print TBL "$name $ego{$name}";
#	close (TBL);
#	system ("TblToFasta $name\.$pat\.tbl > $name\.fa");
#	unlink("$name\.$pat\.tbl");
#	
#	open (QUERY, ">$query\.$pat\.tbl") || die "3 cannot open query: $!";
#	print QUERY "$query $ego{$query}";
#	close (QUERY);
#	system ("TblToFasta $query\.$pat\.tbl > $query\.fa");
#	unlink("$query\.$pat\.tbl");
#       print STDERR "\n#################\nTblToFasta $query\.$pat\.tbl > $query\.fa\n#################\n";
#	
#	close (FILE);
#	
#	system ("pressdb $name\.fa");
#   print STDERR "\n#################\npressdb $name\.fa\n#################\n";
#	system ("wu-tblastx $name\.fa $query\.fa > blast/$cluster\.out");
#    print STDERR "\n#################\nwu-tblastx $name\.fa $query\.fa > $cluster\.out\n#################\n";
#    unlink("$name\.fa\.ntb $name\.fa\.csq $name\.fa\.nhd $name\.fa $query\.fa");
	
	
#    



####
#   Problem with setting hash back to 0, so after a few blasts the $db {$_} is always defined from before I think...
#
#





