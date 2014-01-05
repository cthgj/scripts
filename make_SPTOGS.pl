#! /usr/bin/perl -w

# This script will take a blast outfile (ARGV[0]) and the ego .orth file (ARGV[1]), it will
# find the best subject for each query, then retrieve the TOGs that subject belongs to.
# Finally it will use indexfasta/getfasta to retrieve the sequences of all members of
# each of the TOGs found  and cat them into <QUERY_NAME>_<TOGNAME>.fa

use strict;

my (%subj, %best, %togs , %btogs, %subframe, %secis, %score, %hit);
my ($query, $tog, $subject, $frame);
my @seqs;
my $c=0;
my $ff = 0; 
my $FASTA =  "/seq/databases/EGO/ego8/ego8_112404.seq";
my $INDEX =  "/home/ug/cchapple/research/selenoproteins/ego8/index";


$query= "aa";
#&load_blast();
&load_aln();
&get_secis();
&get_togs();

# sub ignore{
# hit:foreach my $hit (@_)
#     {  

#         ($SUB, $QUERY, $NUM) = split /\*\*/, $hit;
# 	${$scores{$hit}}[0] =~ /Expect[^=]+=\s+(.*?),/ || die("Cannot match evalue <> $.\n");
# 	$evalue = $1;
# 	if ($1 =~ /^\s?e/)
# 	{
# 	    $evalue = "1" . $evalue;
# 	}
	
# 	&debug("current evalue : $evalue lowest evalue : $lowest_evalue");
# 	## Debugging stuff, never mind
# 	if ($debug){$evalue < $lowest_evalue ? print STDERR "Smaller than smallest ? YES\n" : print STDERR "Smaller than smallest ? NO\n";}
# 	if($evalue <= $lowest_evalue)
# 	{
# 	    push @best_hits, $hit;
# 	    $lowest_evalue = $evalue;
# 	}
	    
# 	else
# 	{
# 	    next;
# 	}
	
# 	if($single_best)
# 	    {
# 		$best_hit{$QUERY} = $hit;
# 		$best_SUB = $SUB;
# 		$best_NUM = $NUM;
# 	    }	
#     }
# }

sub load_aln
{
    open(ALN, "$ARGV[0]");
    while(<ALN>)
    {
    next unless />/ || /Query=/  || /Frame/; 
    if(/Query=/)
    {
	/Query=\s(.*?)\s/; 
	next if $1 eq $query;
	$c=0 unless $1 eq $query;
	$ff=0 unless $1 eq $query;
	$query=$1; 
	
    } 
    elsif(/>/ && $c==0)
    {
	/>(.*)/; 
	$subject = $1;
	$subj{$query}=$subject; 
	$best{$subject} = $query;
	print "$query :: $subject\n";
#	my $a = $query . "**" . $subject;
#	print "a : $a\n";
	$c=1;
    }
    elsif (/Frame/ && $ff==0)
    {
	$a = $query . "**" . $subject;
	/Frame\s=\s([+-]\d)/;
	$frame = $1;
	$subframe{$a}=$frame;
	$ff = 1;
    }
    elsif(/>/ && $c==1)
    {
	/>(.*)/;
	$hit{$1}++;
    }
    else{next;}

#     else
#     {
# 	/Score\s+=\s+(\d+)/;
# 	print "$1 : score\n";
# 	my $a = $query . "**" . $subject;
# 	$score{$a} = $1;
# 	$best{$subject} = $query if $score{$a} <  $1;
#     }
}
    
}


sub load_blast 
{
    open(BLAST, $ARGV[0]) || die("Need a BLAST file!!, $ARGV[0] : $!\n");
    while(<BLAST>)
    {
	next unless />/ || /Query=/ || /Score/; 
	if(/Query=/)
	{
	    /Query=\s(.*?)\s/; 
	    $query=$1; $c=0;
	} 
	# get the best hit for query
	if(/>/ && $c==0)
	{
	    />(.*)/; 
	    $subject = $1;
	    $subj{$query}=$subject; 
	    $best{$subject} = $query;
	    
	    $c=1;
	}

	## get frame, but only for the first
	if (/Frame/ && $ff==0)
	{
	    /Frame\s=\s([+-]\d)/;
	    $frame = $1;
	    $subframe{$query}=$frame;
	    $ff = 1;
	}
	else{next}
    }    
    my @queries=keys(%subj); 
    map{print "$_\t:\t$subj{$_}\n"}@queries;
    close(BLAST);
}

sub get_secis
{
    open(SECIS, "$ARGV[2]");
    while(<SECIS>)
    {
	next unless />/;
	/>(.*?):/;
	$secis{$1}++;	    
    }
}



sub get_togs
{
    # open orth file and get TOGS
    my $jj = 0;
    open(ORTH,$ARGV[1]); 
    while(<ORTH>)
    {
	if(/>/){
	    ## Sequence names of this tog = $togs{tog}, unless this is the first
	    $togs{$tog} = [@seqs] unless $jj == 0;
	    $jj=1;
	    @seqs = ();
	    />(.*)/;
	    $tog = $1;
	    next;
	}
	else
	{
	    /(.*?)\s/;
	    push @seqs, $1 if $secis{$1} && $hit{$1};
	    
	    ## If this sequence is one of the best hits from blast, keep
	    ## and get the togs it belongs to
	    if (exists $best{$1}){
		push @{$btogs{$1}}, $tog;
	    }
	    else{next}
	}
    }
    print "-----------------------------------------------------\n";

    my @keys=keys(%btogs);

    ## Print out the TOGS of each best hit
    map {print "$best{$_} : $_ : @{$btogs{$_}}\n" }@keys; 

    ## Get the sequences of the TOG
    
    ## $key == name of best subject seq
    my @name_lists;
    foreach my $key (@keys) {
	print STDOUT "key : $key\n";
	foreach my $TOG (@{$btogs{$key}}){
	    print STDOUT "TOG : $TOG\n";
	    my $NAME = $best{$key} . "_" . $TOG ;
	    open(LIST, ">$NAME");
	    map {print LIST "$_\n"}@{$togs{$TOG}};
	    close LIST;
	    push @name_lists, $NAME;
	}
    }
  
    
    system ("retrieveseqs.pl -ivo $FASTA @name_lists");
#    system ("")
}



### Almost works, not quite, seem to be fewer files than expected




