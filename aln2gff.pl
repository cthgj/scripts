#!/usr/local/bin/perl -w

use strict;

use Getopt::Std;


my %opts;
getopts('m',\%opts) || do { print "Invalid option\n"; exit(1); };




$ARGV[0] =~ /-(.*?)\.aln/;
my $genome = $1; 





my %genomes = (
	       pseudoobscura => '~/research/seq/genomes/chrUn.fa',
	       yakuba        => '~/research/seq/genomes/D.yakuba_unfinished/3R_X.fa',
	       simulans      => '/seq/genomes/D.simulans_unfinished/contigs.bases'
	       );

my $min = $opts{m} || 0; ## return gff starting at line with * till end of HSP

my %queries;
my ($subject, $count, $start, $end,$frame,$strand,$query, $qstart, $qend);
my @hits;

while(<>)
{

    if((/\#/) && ($start))
    { 
	my ($e, $s);
	if($start>$end){$s=$end; $e=$start}
	else{ $s = $start; $e = $end}

	my $gff = "$subject\t.\t.\t$s\t$e\t1\t$strand\t.\t.\t.\n";
#	$queries{$query} = $gff;
	push @hits, $gff;
	my($qs,$qe);

	if($qstart>$qend){ $qs=$qend;  $qe=$qstart}
	else{ $qs = $qstart; $qe = $qend}

	my $qgff = "$query\t.\t.\t$qs\t$qe\t1\t+\t.\t.\t.\n";

#	open(Q, ">$query.gff");
#	print Q "$qgff";
#	close(Q);

      
    }
    if (/Query=/)
    { 
	/Query=\s(.*?)-/ || die("no match!\n"); 
	$query = $1;
	
    }
    if(/Query:/)
    {
	if(($min) && (/\*/))
	{
	    /:\s(\d+)\s/;
	    $qstart = $1;
	}
	/(\d+)[^\d]*$/;
	$qend = $1;
    }

    
    
    if (/>/) { />(.*)/; $subject = $1; $count=0;}
    if(/Frame/){/=\s?(.)(\d)/; $frame = $2-1; $strand = $1;}
    
    if(/Sbjct/)
    {
#	$query{$subject} = $query;
	if($count==0)
	{
	    /:\s(\d+)\s/;
	    $start = $1;
	    $count++;
	}
	if(($min) && (/\*/))
	{
	    /:\s(\d+)\s/;
	    $start = $1;
	}
	/(\d+)[^\d]*$/;
	$end = $1;
	
   
    }

	       
}
my @lin;
my %ka;
my $prev;


 @hits = sort(@hits);
#print STDOUT "h : \n@hits\n"; die();

foreach my $g (@hits)
{
    
    $g =~ /^(.*?)\t/; 
  #  print STDOUT "g : $g";
  #  print STDOUT "prev : $prev\n";
    if($prev)
    {
	if($1 eq $prev)
	{
	    push @lin, $g;
	   
	}
	else
	{
	    #print STDOUT "p : $prev\n";
	    push @lin, $g;
	    $ka{$prev} = [@lin];
	    @lin = (); 
#	    print STDOUT "i :  @{$ka{$prev}}"; 
	    $prev= $1; 
	# print STDOUT "p1 : $prev\n";
	}
	
    }
    else{$prev = $1; push @lin, $g;}
    $ka{$prev} = [@lin];
   
}


my @keys = keys(%ka);


my @queries = keys(%queries);

#foreach my $q (@queries)
#{
#open(TMP, ">>$$\.gff");
#print TMP "$queries{$q}";

#print STDOUT "$queries{$q}";
foreach my $key (@keys){
   # print STDOUT "k : $key\n";
map{ print STDOUT "$_"}@{$ka{$key}};
}
#system("SSgff.chr $genomes{$genome} $$\.gff >> $q\.fa");
#    print STDOUT "SSgff.chr $genomes{$genome} $$\.gff >> $q\.fa\n";
#	my $gff = "$subject\t.\t.\t$s\t$e\t1\t$strand\t.\t.\t.\n";

# $queries{$q} =~ /(.*?)\t(.*?)\t(.*?)\t(.*?)\t(.*?)\t(.*?)\t(.*?)\t(.*?)\t(.*?)\t(.*?)/;
# my $nam = $1;
# my $st  = $4;
# my $end = $5;
# my $str = $7;
    

#close(TMP);
#unlink("$$\.gff");
#}



# foreach my $key (@keys)
# {
#     open(TMP, ">>$$\.gff");
#     map{ print TMP "$_"}@{$ka{$key}};
# #    print STDOUT "SSgff.chr putative.fa $query{$key}\.gff > $query\.fa\n"; 
#     print STDOUT "SSgff.chr $genomes{$genome} $$\.gff >> $query\.fa\n"; 
#     system("SSgff.chr $genomes{$genome} $$\.gff >> $query\.fa");
#     close(TMP);
#     unlink("$$\.gff");
#     map{ print STDOUT "$_"}@{$ka{$key}};
# }


