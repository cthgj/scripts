#!/usr/bin/perl -w
use strict;
use Getopt::Std;
use Cwd; 

my (%all_terms,%ancestors,%stats,%gos,%go_count,%opts,%counted,%go_counter);
getopts('d:S:s:v',\%opts);
my $verbose=$opts{v}||undef;
my $split=$opts{s}||undef; ## split output into 8 files for petitbonum
my $species=$opts{S} || ""; ## suffix for hyper files (species), eg hyper_mgi_0
my $subonto="P";
my $gaf_file=$ARGV[0]||"gene_association.goa_human";
my $geneology_file=$ARGV[1]||"biological_process.geneology";
my $prot_name="Bob";
my $total=0;
my $data_dir=$opts{d}||cwd();


## Read geneology file
open(A,"$geneology_file")||die("cannot open $geneology_file: $!\n");
while(<A>){
    chomp;
    my @terms=split(/\t/);
    my $child=shift(@terms);
    @{$ancestors{$child}}=@terms;
}
close(A);

## Read annotations, GAF file
open(A,"$gaf_file")||die("cannot open $gaf_file: $!\n");
while(<A>){
    next if /^!/;
    chomp;
    my @tmparray=split(/\t/);
    next unless $tmparray[8] eq $subonto;
    next unless $tmparray[3]=~/^$/;  ## skip the NOT annotations
    my $nn=$tmparray[2]; ## $nn == protein name
    $gos{$nn}{$tmparray[4]}++; ## $tmparray[4] == GO
}
close(A);

my @prots=keys(%gos);
my $prot_number=scalar(@prots);
my $c=0;

## Inherit parent terms
foreach my $prot (@prots) {
    $c++;
    printf(STDERR "[Proteins: $c of $prot_number]\r") if $c % 100 == 0;

    ## Count ONLY proteins with >1 direct annotation
    if (scalar(keys(%{$gos{$prot}})) >1){ 
	$total++ ;
    }
    foreach my $goterm (keys(%{$gos{$prot}})){
	$go_count{$goterm}{$prot}++;
	map{$gos{$prot}{$_}++; $go_count{$_}{$prot}++;}@{$ancestors{$goterm}};
    }
    
    
## @GOs is now all the terms of $prot
    my @GOs=keys(%{$gos{$prot}});
    for (my $n=0; $n<=$#GOs;$n++){
	for (my $k=$n+1; $k<scalar(@GOs);$k++){
	    my $go1=$GOs[$n];
	    my $go2=$GOs[$k];
	    $go1=~/GO:(\d+)/||die("no match $go1\n");
	    my $l=$1;
	    $go2=~/GO:(\d+)/||die("no match $go2\n");
	    my $k=$1;
	    my $GOpair;
	    $k>$l ? ($GOpair=$go1 . "xxx" . $go2) : ($GOpair=$go2 . "xxx" . $go1);
	    $stats{$GOpair}++;
	}##end for my $k
    }##end for my $n
}##end foreach prot

printf(STDERR "[Proteins: $c of $prot_number]\r");

my @GOs=keys(%go_count);
my $tt=scalar(@GOs);
print STDERR "\n-------$tt GOs----------\n";

my $file=-1; ## counter to print into different files
my $file_name;
my $limit=$tt/8;
for (my $n=0; $n<=$#GOs;$n++){
    printf(STDERR "[GOs: $n of $tt]\r") if $n % 100 == 0;
#    unless (($GOs[$n] eq 'GO:0015031') || ($GOs[$n] eq 'GO:2000112')){next}

## Print to 8 different files to use the 8 processors on
## petitbonum
    if (($split) &&  ($n % $limit == 0)){
	$file++;
	$file_name=$split . "_" . $file ;	
	if ($n == 0){print STDERR $n+1 . "($n) st opening...\n";open(FILE,">$file_name")}
	else{
	    close(FILE);
	    print STDERR $n+1 . "($n)nd opening...\n";
	    open(FILE,">$file_name");
	};
	
    }

    for (my $k=$n+1; $k<=$#GOs;$k++){
#	unless (($GOs[$k] eq 'GO:0015031') || ($GOs[$k] eq 'GO:2000112')){next}
	    my $go1=$GOs[$n];
	    my $go2=$GOs[$k];
	    # unless ((($go1=~/0007411/) && ($go2=~/0085020/)) ||
	    # 	  (($go1=~/0085020/) && ($go2=~/0007411/))){next}
	    $go1=~/GO:(\d+)/;
	    my $l=$1;
	    $go2=~/GO:(\d+)/;
	    my $k=$1;
	    my $GOpair;
	    $k>$l ? ($GOpair=$go1 . "xxx" . $go2) : ($GOpair=$go2 . "xxx" . $go1);
	   # print "$GOpair\n";
	    my $go1_count=scalar(keys(%{$go_count{$go1}}));
	    my $go2_count=scalar(keys(%{$go_count{$go2}}));
	    if($go1_count==0 || $go2_count==0){
	    	print STDERR  "PROBLEM : $GOpair\n$stats{$GOpair},$go_count{$go1},$go_count{$go1},$total\n";
	    	exit();
	    }
	    my $gopair_count=0;
	    defined($stats{$GOpair}) ? ( $gopair_count=$stats{$GOpair}) : ($gopair_count=0);
            ## OK, now deal with all go pairs
            # R stuff
            ## $gocount : number of times these 2 GOs coincide
            ## $go_count{$go1} : number of prots with go1 
            ## $total-$go_count{$go1} : number of prots without go1 
            ## $go_count{$go2} : number of prots with go2

	    if ($split){
	    	print FILE "hhyper=phyper($gopair_count,$go1_count," . ($total-$go1_count) .",$go2_count,lower.tail = FALSE)\n";
	    	print FILE "lhyper=phyper($gopair_count,$go1_count," . ($total-$go1_count) .",$go2_count,lower.tail = TRUE)\n";
	    	print FILE "hypers=c(\"$go1\",\"$go2\",lhyper,hhyper)\n";
	    	print FILE "write(hypers,file=\"$data_dir/hyper_" . $species . "_" . "$file\",ncolumns=4,append=TRUE,sep=\"\t\")\n";
	    }
            else{
            #                                 #go1go2     #go1                 #!go1              #go2
	    	print STDOUT "hhyper=phyper($gopair_count,$go1_count," . ($total-$go1_count) .",$go2_count,lower.tail = FALSE)\n";
	    	print STDOUT "lhyper=phyper($gopair_count,$go1_count," . ($total-$go1_count) .",$go2_count,lower.tail = TRUE)\n";
	    	print STDOUT "hypers=c(\"$go1\",\"$go2\",lhyper,hhyper)\n";
	    	print STDOUT "write(hypers,file=\"$data_dir/hyper\",ncolumns=4,append=TRUE,sep=\"\t\")\n";
	    }
    }
printf(STDERR "[GOs: $n of $tt]\n") if $n == $#GOs;
}
exit();




# sub load_network{
#     print STDERR "Loading Network...\n" if $verbose;
#     my $A = IO::File->new("< $network_file")|| die("Cannot open $network_file : $!\n");
#     while(<$A>){
# 	next if /^\d+$/;
# 	my ($bait,$target)=split(/\s+/,$_);
#     }
