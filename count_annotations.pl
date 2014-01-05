#!/usr/bin/perl -w

use strict;
use Getopt::Std;
use Cwd; 

my (%all_terms,%stats,%go_count,%opts,%counted,%go_counter,%ancestors,%PP,%synonyms);
$synonyms{'LOADED'}=0;
getopts('hps:w:',\%opts);
&usage() if $opts{h};
my $subonto=$opts{s} || "P";
my $gaf_file=$ARGV[0]||"gene_association.goa_human";
my $geneology_file=$ARGV[1]||"biological_process.geneology";
my $prot_name="Bob";
my $total_gt1=0; ## num of prots with >1 direct annotations.
my $print_annots_per_prot=$opts{p}||undef;
my $synfile=$ARGV[1]||die("Need a synfile as \$ARGV[1]\n");
my $wanted=$opts{w}||undef;
## If we just want to count how many DIRECT GOs
## each protein has
if ($print_annots_per_prot){
    my %want;
    if($wanted){
	open(W,"$wanted")||die("Could not open $wanted:$!\n");
	while(<W>){
	    chomp;
	    my $n=&get_name($_);
	    $want{$n}++;
	}
    }
    open(A,"$gaf_file")||die("cannot open $gaf_file: $!\n");
    my %gos;
    while(<A>){
	next if /^!/;
	chomp;
	my @tmparray=split(/\t/);
	next unless $tmparray[8] eq $subonto;
	next unless $tmparray[3]=~/^$/;  ## skip the NOT annotations

	my $prot=&get_name($tmparray[1]); ## $prot == protein name
	$want{$prot}=1 unless $wanted;
	$gos{$prot}{$tmparray[4]}++ if $want{$prot}; ## $tmparray[4] == GO
    }
    foreach my $prot (keys(%want)){
	my @k=keys(%{$gos{$prot}});
	print "$prot\t" . scalar(@k) . "\n";
	
    }
}
else{
    %ancestors=&read_geneology($geneology_file);
#my (%gos,%dgos)=&load_annotations($gaf_file);
    my ($a,$b)=&load_annotations($gaf_file);
    my $prot_number=&get_stats($a,$b);
    print "TOTAL : $prot_number\n";
    print "TOGT1 : $total_gt1\n";
#&calculate_probs($prot_number);
    
#map{print "$_ : $gos{$_}\n"}keys(%gos);
    
}

sub read_geneology{
## Read geneology file
    my $file=shift;
    open(A,"$file")||die("cannot open $file: $!\n");
    while(<A>){
	chomp;
	my @terms=split(/\t/);
	my $child=shift(@terms);
	@{$ancestors{$child}}=@terms;
    }
    close(A);
    return(%ancestors)
}
sub load_annotations{
## Read annotations, GAF file
    my $file=shift;
    my(%gos,%dgos);
    open(A,"$file")||die("cannot open $file: $!\n");
    while(<A>){
	next if /^!/;
	chomp;
	my @tmparray=split(/\t/);
	next unless $tmparray[8] eq $subonto;
	next unless $tmparray[3]=~/^$/;  ## skip the NOT annotations
	my $prot=$tmparray[2]; ## $prot == protein name
	$gos{$prot}{$tmparray[4]}++; ## $tmparray[4] == GO
	$go_count{$tmparray[4]}{$prot}++; ## $tmparray[4] == GO
	
	$dgos{$prot}{$tmparray[4]}++;
#	print 	"\$dgos{$prot} : $dgos{$prot}\n";
	
	map{$gos{$prot}{$_}++; $go_count{$_}{$prot}++;}@{$ancestors{$tmparray[4]}};
    }
    close(A);	
    return(\%gos,\%dgos);
}

sub get_stats{
    my %gos=%{$_[0]};
    my %dgos=%{$_[1]};
    my @prots=keys(%gos);
    my $prot_number=scalar(@prots);
    my $c=0;
    foreach my $prot (@prots) {
	$c++;
	## Count proteins with >1 direct annotation
	if (scalar(keys(%{$dgos{$prot}})) > 1){ 
	    my $is_anc=0;
	    my $start_num=scalar(keys(%{$dgos{$prot}}));
	  go1:foreach my $go1 (keys(%{$dgos{$prot}})){
	    go2:foreach my $go2 (keys(%{$dgos{$prot}})){
		next if $go1 eq $go2;
		foreach my $anc (@{$ancestors{$go1}}){
		    if ($anc eq $go2){
			$is_anc=1;
			$start_num--;
			#print "$prot : $go1 $go2\n";
			#next go1;
		    }
		}
	    }
	      
	  }
	    $total_gt1++ if $start_num>0;
	}
	
    }##end foreach prot
    return($prot_number);
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
		if(/^([^\t]+)\t([^\t]+)$/){die();
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

 
sub usage{
    $0=~s/.+?([^\/]+)$/$1/;
    print <<EndOfHelp;
    
    $0 <GAF> <geneology>

  USAGE:  
	$0 parses a gene_association file and a geneology file and
        returns the total number of proteins, and the number of proteins
	with >1 DIRECT annotation (it will not count a GO if the protein
	it annotates is also annotated to a child term).
   
  OPTIONS:
	-h : Print this help and exit.
	-s : Subontology to use (def:P)
    

EndOfHelp
exit();
}
