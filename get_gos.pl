#!/usr/bin/perl -w

use strict;
use Getopt::Std;

my %opts;
my (%synonyms, %proteins,%ancestors);
getopts('hvo:a',\%opts)|| do { print "Invalid option, try '$0 -h' for more information\n"; exit(1); };
unless($ARGV[0]){
    print STDERR "\nUSAGE:\n\t$0 <prot_name/prot_name_file> <GAF_file> <synonyms_file>\n\n"; exit;
}
my $subonto=$opts{o}||".";
my $ancestors=$opts{a}||undef;
$synonyms{LOADED}=0;

my $synfile=$ARGV[2]||die("Need a synfile as \$ARGV[2]\n");
die("Need a file with protein names or a single protein name as ARGV[0]\n") unless $ARGV[0];
if (-e $ARGV[0]){
    open(A,"$ARGV[0]")||die("Could not open $ARGV[0] :$!\n");
	while(<A>){
	    chomp;
	    my $name=&get_name($_,"main");
	    $proteins{$name}{NAME}++
    }
}
else{
    my $name=&get_name($ARGV[0],"main");
    $proteins{$name}{NAME}++;  
}

&parse_gaf_file();
foreach my $prot (keys(%proteins)){
    print "$prot\t $proteins{$prot}{NAME} $proteins{$prot}{GOs}\n";
    my %hash=$proteins{$prot}{GOs};
    my @a=keys(%hash);
}

############################################################
sub parse_gaf_file{    
    open(A,"$ARGV[1]")|| die("Cannot open $ARGV[1]:$!\n");
    while(<A>){
	next if /^!/;
	chomp;
	my @tmparray=split(/\t/,$_);
	next unless $tmparray[8] =~/$subonto/;
	next if $tmparray[3] !~/^$/;
	my $name=&get_name($tmparray[1],"parse_gaf_file");
	next unless defined($proteins{$name}{NAME});
	$proteins{$name}{GOs}{$tmparray[4]}="DIR";
	## Inherit ancestor annotations
	map{$proteins{$name}{GOs}{$_}="ANC"}@{$ancestors{$tmparray[4]}};
	
    }
    close(A);
}

############################################################
sub get_name{
    my $name=$_[0];
    my $old_name=$name;
    my $called_by=$_[1]||"null";
    my $reverse=$_[2]||undef;
    if ($synonyms{LOADED}==0){
	if(-e $synfile){
	    open(S,$synfile);
	    while(<S>){
		#if this is a GAF specific synonyms file, i.e.
		# network_name\tGAF_name
		if(/^(.+)\t(.+)/){
		    $synonyms{NAME}{$2}=$1; 
		    $synonyms{NAME}{$1}=$1; 
		}
		## If this is a type 1 synonyms file, ie, 
		## network_name\tGAF_name\tAccession
		elsif(/^(.+)\t(.+)\t(.+)/){
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
	else{die("Could not open $synfile\n" );	}	
 	$synonyms{LOADED}=1;
    }
   if(defined($synonyms{NAME}{$name})){
	$name=$synonyms{NAME}{$name} ;
    }
   
    return $name;
   
}
