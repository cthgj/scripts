#!/usr/bin/perl -w
use strict;
use Getopt::Std;
use Cwd; 

my (%all_terms,%synonyms,%ancestors,%stats,%gos,%go_count,%opts,%counted,%go_counter);
getopts('n:g:d:S:D:a:m:s:v',\%opts);
my $verbose=$opts{v}||undef;
my $split=$opts{s}||undef; ## split output into 8 files for petitbonum
my $species=$opts{S} || ""; ## suffix for hyper files (species), eg hyper_mgi_0
my $subonto="P";
my $gaf_file=$opts{a}||"gene_association.goa_human";
my $geneology_file=$opts{g}||"biological_process.geneology";
my $synfile=$opts{m}||die("need a map file (-m)\n");
my $prot_name="Bob";
my $total=0;
my $data_dir=$opts{D}||cwd();
my $network_file=$ARGV[0];
$synonyms{LOADED}=0;
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


my %proteins=&load_network($network_file);
open(A,"$gaf_file")||die("cannot open $gaf_file: $!\n");
while(<A>){
    next if /^!/;
    chomp;
    my @tmparray=split(/\t/);
    next unless $tmparray[8] eq $subonto;
    next unless $tmparray[3]=~/^$/;  ## skip the NOT annotations
    my $nn=&get_name($tmparray[1],"GAF"); ## $nn == protein name
    next unless defined($proteins{$nn});
    $gos{$nn}{$tmparray[4]}++; ## $tmparray[4] == GO
}
close(A);

my $prot_number=(keys(%proteins));
my $c=0;
## Inherit parent terms
foreach my $prot (keys(%proteins)) {
    $c++;
    printf(STDERR "[Proteins: $c of $prot_number]\r") if $c % 100 == 0;

    foreach my $goterm (keys(%{$gos{$prot}})){
	$go_count{$goterm}{$prot}++;
	map{$gos{$prot}{$_}++; $go_count{$_}{$prot}++;}@{$ancestors{$goterm}};
    }
    

}
print STDERR "\n";
$c=0;
my $interactions=0;
my %seen;
foreach my $prot (keys(%proteins)){
    printf(STDERR "[Proteins: $c of $prot_number]\r");
    $c++;
    target:foreach my $target (keys(%{$proteins{$prot}{INTERACTORS}})){
	next target if defined($seen{$prot . "_" . $target});
	$seen{$prot . "_" . $target}=	$seen{$target . "_" . $prot} = 1;
	$interactions++;
	foreach my $bgo (keys(%{$gos{$prot}})){
	    foreach my $tgo (keys(%{$gos{$target}})){
		$bgo=~/GO:(\d+)/||die("no match $bgo\n");
		my $l=$1;
		$tgo=~/GO:(\d+)/||die("no match $tgo\n");
		my $k=$1;
		my $GOpair;
		$k>$l ? ($GOpair=$bgo . "_" . $tgo) : ($GOpair=$tgo . "_" . $bgo);
		$stats{$GOpair}++;
	    
	    }
	}
    }
}

printf(STDERR "[Proteins: $c of $prot_number]\r");

my @GOs=keys(%go_count);
my $tt=scalar(@GOs);
print STDERR "\n-------$tt GOs----------\n";

my $file=-1; ## counter to print into different files
my $file_name;
my $limit=$tt/8;
for (my $n=0; $n<=$#GOs;$n++){
    printf(STDERR "[GOs: $n of $tt]\r") if $n % 100 == 0;

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
	    $k>$l ? ($GOpair=$go1 . "_" . $go2) : ($GOpair=$go2 . "_" . $go1);
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

	  
            #                                 #go1go2     #go1                 #!go1              #go2
	    	print STDOUT "hhyper=phyper($gopair_count,$go1_count," . ($interactions-$go1_count) .",$go2_count,lower.tail = FALSE)\n";
	    	print STDOUT "lhyper=phyper($gopair_count,$go1_count," . ($interactions-$go1_count) .",$go2_count,lower.tail = TRUE)\n";
	    	print STDOUT "hypers=c(\"$go1\",\"$go2\",lhyper,hhyper)\n";
	    	print STDOUT "write(hypers,file=\"$data_dir/hyper\",ncolumns=4,append=TRUE,sep=\"\t\")\n";
	    
    }
printf(STDERR "[GOs: $n of $tt]\n") if $n == $#GOs;
}
exit();


###############################################################
sub load_network {
    my $network_file=shift;
        my %proteins;
    print "nn : $network_file\n";
    print STDERR "Loading Network...\n" if $verbose;
    open(N,"$network_file")|| die("Cannot open $network_file : $!\n");
    while(<N>){
	next if /^\d+$/;
	next if /^!/;
	next if /^\s*$/;
	my ($bait,$target)=split(/\s+/,$_);
	$bait=&get_name($bait,"load_network");
	$target=&get_name($target,"load_network");
	## %proteins will hold all the info on all the 
	## network's proteins
	$proteins{$bait}{INTERACTORS}{$target}++;
	$proteins{$target}{INTERACTORS}{$bait}++;
    }
    close(A);
    return(%proteins);
}

############################################################
sub get_name{
    my $name=$_[0];
    my $old_name=$name;
    my $suffix='_HUMAN';
    my $called_by=$_[1]||"null";
    my $reverse=$_[2]||undef;
    if ($synonyms{LOADED}==0){
	if(-e $synfile){
	    open(S,$synfile);
	    while(<S>){
		## If this is a type 1 synonyms file, ie, 
		## network_name\tcurrent_name\tAccession
		if(/^(.+)\t(.+)\t(.+)/){
		    $synonyms{ACC}{$3}=$3;
		    $synonyms{ACC}{$2}=$3;
		    $synonyms{ACC}{$1}=$3;
		    $synonyms{NAME}{$1}=$2;
		    $synonyms{NAME}{$2}=$2;
		    $synonyms{NAME}{$3}=$2;
		}	
                ## If this is a type 2 synonyms file, ie, 
		## gene name\tlist of synonyms
		elsif(/^(.+)\t(.+)/){
		       my $nn=$1;
		       ## Discard protein if there are naming problems
		       my @syns=split(/\s/,$2);
		       $synonyms{GOOD}{$nn}++;
		       for (my $n=0;$n<scalar(@syns); $n++){
			   $synonyms{$syns[$n]}=$nn;
		       }
		   }
		   else{
		       die("Bad format synonyms file : $synfile\n");
		   }
	    }
	    close(S);
	}
	else{die("Could not open $synfile\n" );
	    # print STDERR "No synonyms file found, parsing annotations... \n NEED TO MODIFY THIS BECAUSE OF THE CHANGE IN SYNONYMS FORMAT" if $verbose;
	#     open(A,"$ARGV[0]")|| die("Cannot open synonyms $ARGV[0]:$!\n");
	#     while(<A>){
	# 	my $silivar="!";
	# 	next if /^$silivar/; ## silly thing cause emacs screws up the syntax highlihting when /!/;
	# 	my @tmparray=split(/\t/);
	# 	map{$synonyms{$_}=$tmparray[2]}split(/\|/,$tmparray[10]);
	#     }
	#     close(A);
	 }	
 	$synonyms{LOADED}=1;
    }
    my $hname;
    $name =~ /_$suffix/i ? ($hname=$name) : ($hname=$name . "_" . $suffix);
    # if(defined($synonyms{ACC}{$name})){
    # 	$name=$synonyms{ACC}{$name};
    # }
    # elsif(defined($synonyms{ACC}{$hname})){
    # 	$name=$synonyms{ACC}{$hname};
    # }
    # elsif(defined($synonyms{GOOD}{$hname})){die("111 : $name\n");
    # 	$name=$synonyms{GOOD}{$hname};
    # }
    if(defined($synonyms{NAME}{$name})){
	$name=$synonyms{NAME}{$name} ;
    }
    elsif(defined($synonyms{NAME}{$hname})){
	$name=$synonyms{NAME}{$hname};
    }
    elsif ($name=~/(.+)_$species/i && defined($synonyms{$1})){die("333 : $name\n");
	$name=$synonyms{$1};
    }
    #else{die("Shit, name problem: $name\n")}
	    

#    print "xxx $old_name:$name\n";
    
    return $name;
   
}
