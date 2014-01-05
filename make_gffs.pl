#!/usr/bin/perl -w
use strict;
use Getopt::Std;
use Term::ANSIColor;

my %opts;
getopts('q:t:p:k:vhg',\%opts);
my $NAMES = $opts{q} || &usage();
my $PROTEINS = $opts{p} || &usage();
my($nn,$t,$ll,$exo,$s,$length);
my $delete = $opts{d} || 0;
my $TARGET = $opts{t} || &usage();
my $verbose = $opts{v} || undef;
my @names;

&usage() if $opts{h};

if (-r "$NAMES")
{
    open(A,"$NAMES");
    @names=<A>
}
else
{
    @names = split(/\s+/,$NAMES);
}
my $tname;
if($TARGET =~ /([^\/]+)\.fa/){$tname = $1;}
else{$tname = $TARGET}
my @details;
foreach my $prot(@names)
{
    @details = ();

    my $not_found = 0;
    chomp($prot);
    my $best_score = 0;
    # retrieve protein sequence
    print STDERR "Retrieving sequence $prot from $TARGET... " if $verbose;
    system("retrieveseqs.pl -fn $PROTEINS $prot | perl -ne 'if(/>/){print} else{s/U/\*/g; print;}'> /tmp/$prot.pep");
    print STDERR "done\n" if $verbose;
    
    $exo=$prot . "-" . $tname . ".exo"; 
    my $wise = $prot . "-" . $tname . ".wise";
    my $gff = $prot . "-" . $tname . ".gff";
    if(-r "$exo")
    {

	print STDERR color 'bold red';
	print STDERR "CAUTION : Reusing existing exonerate output file $exo\n" if $verbose;
	print STDERR color 'reset';
    }
    else
    {
	print STDERR "Running exonerate for $prot on $TARGET... ";
	system("exonerate -m protein2genome -q /tmp/$prot.pep -t $TARGET > $exo"); 
	print STDERR "done\n";
    }
    
    open(A,"$exo");
    my ($is_minus_strand,$start,$range) = (0,0,0,0);
   
    my $score = 0;
    while(<A>)
    {
	if(/Raw\sscore:\s?(\d+)/)
	{
	    $score=$1;
	}
	elsif(/Target\srange:\s(\d+)\s->\s(\d+)/)
	{
	    if($1 > $2)
		{
		    $is_minus_strand = 1;
		    $start = $2 -1000;
		    $range = $1-$2+2000;
			
		}
	    else
		{
		    $is_minus_strand = 0;  
		    $start = $1 -1000;
		    $range = $2-$1+2000;
		}
	    
	     if ($score>=$best_score)
	     {
		 $best_score = $score ;
		 @details = ($is_minus_strand,$start,$range); 
	     }
	}
	else{next}
    }
    unless ($details[0])
    {
	print STDERR color 'bold blue';
	print STDERR "$prot not found\n";
	print STDERR color 'reset';
	next;
    }
    next if $not_found == 1;
    if($details[0] == 0)
	{	    
	    print STDERR "a fastasubseq -f $TARGET -s $details[1] -l $details[2] > j\n";
	    system("fastasubseq -f $TARGET -s $details[1] -l $details[2] > j;");
	    print STDERR "Running Genewise for $prot on $TARGET (subseq $details[1]-" . ($details[1] + $details[2]) . ", $details[2] nt)... ";
	    system("genewise -quiet -pretty -cdna -pep -gff /tmp/$prot.pep j > $wise");
	    print STDERR "done\n";
	    unless ($opts{g})
	    {
		system("grep cds $wise > ttt");
		open(AA,"ttt");
		open(GFF,">$gff");
		while(<AA>)
		{
		    /seq\((\d+)/; 
		    my $a=$1;
		    /cds\s+(\d+)\s+(\d+)/; 
		    my $s=$1+$a; 
		    my $e=$2+$a;
		    s/cds\s+(\d+)\s+(\d+)/cds\t$s\t$e/; 
		    print GFF;
		}
		close(AA);
		close(GFF);
		system("gawk 'BEGIN{OFS=\"\t\"}{\$1=\"$tname\"; \$NF=\"$prot\"; print}' $gff > ttt; mv ttt $gff");
	    }
	   
	}
	elsif($details[0] == 1)
	{
	    print STDERR "fastasubseq -f $TARGET -s $details[1] -l $details[2] > j\n";
	    system("fastasubseq -f $TARGET -s $details[1] -l $details[2] > j;");
	    print STDERR "Running Genewise for $prot on $TARGET (subseq $details[1]-" . ($details[1] + $details[2]) . ", $details[2] nt)... ";
	    system("genewise -quiet -pretty -cdna -pep -gff -trev /tmp/$prot.pep j > $wise");
	    print STDERR "done\n";
	    unless ($opts{g})
	    {
		system("grep cds $wise > ttt");
		open(AA,"ttt");
		open(GFF,">$gff");
		while(<AA>)
		{
		    /seq\((\d+)/; 
		    my $a=$1;
		    /cds\s+(\d+)\s+(\d+)/; 
		    my $s=$1+$a; 
		    my $e=$2+$a;
		    s/cds\s+(\d+)\s+(\d+)/cds\t$s\t$e/; 
		    print GFF;
		}
		close(AA);
		close(GFF);
		system("gawk 'BEGIN{OFS=\"\t\"}{\$1=\"$tname\"; \$NF=\"$prot\"; print}' $gff > ttt; mv ttt $gff");
	    }
	}
    else{die("aa : @details\n")}
    system("rm /tmp/$prot.pep");

    if ($delete =~ /e/)
    {
	 system("rm $exo");
    }
    elsif ($delete =~ /w/)
    {
	 system("rm $wise");
    }
  
  
}	

    

sub usage
{
    print STDERR <<EndOfHelp;

    This script will take as input a list of protein names (either a file 
    containing one name per line, or a quoted, space separated list of names),
    a multi-fasta file containing their sequences and a nucleotide fasta file. 
    It will sequentially extract each of the protein sequences from the file,
    use exonerate to align them to the nucleotide target and then genewise to fine 
    tune the alignment based on the highest scoring exonerate result and produce gff output.

    USAGE : make_gffs.pl [-vgd] -q <QUERY_LIST> -t <TARGET_FASTA_FILE> -p <FASTA_QUERIES>
    
    COMMAND-LINE OPTIONS:

    -h : Print this help and exit.
    -d : Specify which files you want to delete [e/w]. For example "-k ew" 
         will tell the program to delete both the Exonerate and the geneWise
         output files. By default all files are kept.
    -g : Do NOT return GFF files.	 
    -p : MultiFasta file containing the queries whose names are given by
         the -q option (any other sequences will be ignored). 
    -q : Queries to extract from the -p MultiFasta file.
    -t : Target (nucleotide) to align queries to.
    -v : Verbose output

EndOfHelp

exit(0);

}
