#!/usr/bin/perl -w
## Get the numbers necessary to calculate the probabilities of association of each GO pair

use strict;
use Getopt::Std;
use Cwd; 

my (%stats,%go_count,%opts,%counted,%go_counter);
getopts('Rt:g:d:S:s:v',\%opts);
my $verbose=$opts{v}||undef;
my $split=$opts{s}||undef; ## split output into 8 files for petitbonum
my $species=$opts{S} || ""; ## suffix for hyper files (species), eg hyper_mgi_0
my $subonto="P";
my $gaf_file=$ARGV[0]||"gene_association.goa_human";
my $geneology_file=$ARGV[1]||"biological_process.geneology";
my $prot_name="Bob";
my $total_gt1=0; ## num of prots with >1 direct annotations.
my $data_dir=$opts{d}||cwd();
my $GOs_file=$opts{g}||undef; ## List of Go terms we are interested in.
if($GOs_file){print STDERR "WARNING -g option is buggy CHECK THE RESULTS!!\n\n";}
my $temp_out_file=$opts{t}||"./calculate_probabilities.tmp";
my %ancestors=&read_geneology($geneology_file);
#my (%gos,%dgos)=&load_annotations($gaf_file);
my ($a,$b)=&load_annotations($gaf_file);
#$GOs_file ? ($c=&load_wanted_GOs($GOs_file)) : ($c=0);   
my %wanted=%{&load_wanted_GOs($GOs_file)} if $GOs_file;
my ($prot_number)=&get_stats($a,$b);
&calculate_probs($b);

#map{print "$_ : $gos{$_}\n"}keys(%gos);


sub read_geneology{
## Read geneology file
    my $file=shift;
    open(A,"$file")||die("cannot open $file: $!\n");
    while(<A>){
	chomp;
	my @terms=split(/\t/);
	my $child=shift(@terms);
	@{$ancestors{ANC}{$child}}=@terms;
	map{push @{$ancestors{KIDs}{$_}},$child}@terms;
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
	$gos{$tmparray[1]}{$tmparray[4]}++; ## $tmparray[1]== name, $tmparray[4] == GO
	$go_count{$tmparray[4]}{$tmparray[1]}++; ## $tmparray[4] == GO
	### careful about dup GOs
	
	$dgos{$tmparray[1]}{$tmparray[4]}++;
#	print 	"\$dgos{$prot} : $dgos{$prot}\n";
	
	map{$gos{$tmparray[1]}{$_}++; $go_count{$_}{$tmparray[1]}++;}@{$ancestors{ANC}{$tmparray[4]}};
    }
    close(A);
    my $total_gt1;
    map{$total_gt1++ if (keys(%{$dgos{$_}})) > 1}keys(%dgos);
  
    return(\%gos,\$total_gt1);
}

sub get_stats{
    my %gos=%{$_[0]};
    my $total_gt1=$$_[1];
    print STDERR "TT : $$_[0] : $$_[1]\n";
    #%wanted=%{$_[2]} if($GOs_file);
    my @prots=keys(%gos);
    my $prot_number=scalar(@prots);
    my $c=0;
    foreach my $prot (@prots) {
	$c++;
	$verbose && do {
	    printf(STDERR "[Proteins: $c of $prot_number]\r") if $c % 100 == 0;
	};
	## Count proteins with >1 direct annotation
	
	my @GOs=keys(%{$gos{$prot}});
	for (my $n=0; $n<=$#GOs;$n++){
	    for (my $k=$n+1; $k<scalar(@GOs);$k++){
		my @l=sort($GOs[$n],$GOs[$k]);
		$stats{$l[0] . "-" . $l[1]}++;
	    }##end for my $k
	}##end for my $n
    }##end foreach prot
    printf(STDERR "[Proteins: $c of $prot_number]\r") if $verbose;
    return($prot_number);
}

sub calculate_probs{
    my $total_gt1=$_[0];
    my @GOs=keys(%go_count);
    my $tt=scalar(@GOs);
    print STDERR "\n-------$tt GOs----------\n" if $verbose;

    my $file=-1; ## counter to print into different files
    my $file_name;
    my $limit=$tt/8;
    # print "\", total_size(), "\n";
    for (my $n=0; $n<=$#GOs;$n++){
	$verbose && do {printf(STDERR "[GOs: $n of $tt]\r") if $n % 100 == 0;};
	for (my $k=$n+1; $k<=$#GOs;$k++){
	    if($GOs_file){
		next unless defined($wanted{$GOs[$n]}{$GOs[$k]});
	    }
	    my @l=sort($GOs[$n],$GOs[$k]);
	    my $go1_count=scalar(keys(%{$go_count{$l[0]}}));
	    my $go2_count=scalar(keys(%{$go_count{$l[1]}}));
	    if(scalar(keys(%{$go_count{$l[0]}}))==0 || scalar(keys(%{$go_count{$l[1]}}))==0){
	    	print STDERR  "PROBLEM : " , $l[0] . "-" . $l[1] , "\n" , $stats{$l[0] . "-" . $l[1]}, ",$go_count{$l[0]},$go_count{$l[1]},$total_gt1\n";
	    	exit();
	    }
	    $stats{$l[0] . "-" . $l[1]} = 0 unless defined($stats{$l[0] . "-" . $l[1]});
	    print STDOUT $l[0] . "-" . $l[1] ,"\t$total_gt1\t", scalar(keys(%{$go_count{$l[0]}})), "\t", scalar(keys(%{$go_count{$l[1]}})), "\t",$stats{$l[0] . "-" . $l[1]},"\n" ;
	}
    }
}

sub load_wanted_GOs{
    my $GOs_file=shift;
    my %gos;
    open(A,"$GOs_file")||die("Cannot open $GOs_file : $!\n");
    while(<A>){
	chomp;
	my @a=split(/\s+/);
	$gos{$a[0]}{$a[1]}=$gos{$a[1]}{$a[0]}=1;
	
    }
    
#    map{
	#map{$gos{$_}++}@{$ancestors{ANC}{$_}};
	#map{$gos{$_}++}@{$ancestors{KIDs}{$_}};
 #  }keys(%gos);
    return \%gos;  
    
}
