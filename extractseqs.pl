#!/usr/bin/perl -w
 #=====================================================================================================#
 # This script takes as input a list of ids (one id per line, but can have more than 1 ID file)	       #
 # and reads a sequence file containing these sequences. It will return				       #
 # one file per id list with the sequences of that list (-o flag)				       #
 # Or print them all to STDOUT  (default)							       #
 # If sure that names in id and FASTA file identical (*sure*) use -i flag			       #
 # USAGE : retrieveseqs.pl SEQUENCE_FILE ID_FILE(S)                                                    #
 #=====================================================================================================#


# USAGE : retrieveseqs.pl SEQUENCE_FILE ID_FILE(S)

use strict;
use Data::Dumper;
use Getopt::Std;
my %opts;
getopts('ois',\%opts);

&usage unless ($ARGV[0] && $ARGV[1]);
my ($ESTs,@qtfile) = @ARGV;

my %SEQS = ();
my @qlist = ();
my @tlist = ();

my $output = $opts{o} || 0;
my $identical = $opts{i} || 0;
my $sfile = $opts{s}|| undef;


foreach my $file (@qtfile) {   # foreach file with list of names
    if ($file =~ /\.gz$/) {
        open(F,"zcat $file |") || &usage; 
    } else {
        open(F,"< $file")|| &usage;
    };
    my $query;
    while (<F>) {
        my $id; 
        chomp;
        ($id = $_) =~ s/\s*$//og;    # remove trailing space
        defined(%{ $SEQS{$id} }) ||                           ## If doesn't exist, create "%{ $SEQS{$id}" foreach id with 2 keys
            (%{ $SEQS{$id} } = ( TYPE => [], SEQ => '' ));    ## one a list of the sequence ids (list) the other the seq
        $. == 1 && do {                                    
            $query = $id;  # $query= FIRST id in file
            push @qlist, $id; # qlist = list of 1st ids
	    # next;
        };
        push @{ $SEQS{$query}{TYPE} }, $id; # @{ $SEQS{$query}{TYPE} } is a list fo all ids in this file
    };
    close(F);
};
print STDERR "****\nGot names\n****\n";
@tlist = keys %SEQS;    #list of ALL ids

# print STDERR Data::Dumper->Dump([ \%SEQS ], [ qw( *SEQS ) ])."\n@qlist\n";



$identical == 1 ? &get_identical_seqs() : &get_seqs();
&output();



#################################################################
#################################################################
sub get_identical_seqs
{
    my ($seq_name, $seq, $isseq, $hit);
    $isseq = 0;
    if ($ESTs =~ /\.gz$/) {
	open(FILE,"zcat $ESTs |");
    } else {
	open(FILE,"< $ESTs");
    };
    
    my $count=1;
    
    while(<FILE>) {
	
	
	my ($l,$t);
	next if /^\s*$/o;
	chomp;
	$l = $_;
	$l =~ /^>(.*?)(\s+.*)?$/o && do {
	    $count++;
	    print STDERR "." if $count % 100 ==0; 
	    print STDERR "[$count]\n" if $count % 10000 ==0; 
	    $t = $1;               # $t = sequence name (defined further down)
	    $isseq && do {	
		
		exists($SEQS{$seq_name}) && do { 
		    $SEQS{$seq_name}{SEQ} = $seq; 
		}
	    };
	    $seq_name = $t; # quotemeta($t);                        ## quotemeta will backslash all nonstd chars, so $t = $seq_name
	    $seq = '';
	    $isseq = 1;
	    next;                                             ## next moves to sequence line 
	};
	$seq .= $l;                                           ## $seq now gets sequence ( because of next above)
    }
    
    close(FILE);
    $isseq && do {                                            ## get last sequence
	exists($SEQS{$seq_name}) && do { $SEQS{$seq_name}{SEQ} = $seq; }
    }; 
}
#################################################################
#################################################################
sub get_seqs
{

    my ($seq_name, $seq, $isseq, $hit);
    $isseq = 0;
    if ($ESTs =~ /\.gz$/) {
	open(F,"zcat $ESTs |");
    } else {
	open(F,"< $ESTs") || die("cannot open $ESTs : $!\n");
    };
    
    my $count=1;
    #map{print STDOUT "$_\n"} @tlist; die();
    while(<F>) {
	
	
	my ($l,$t);
	next if /^\s*$/o;
	chomp;
	$l = $_;
	$l =~ /^>(.*?)(\s+.*)?$/o && do {
	    $count++;
	    print STDERR "." if $count % 100 ==0; 
	    print STDERR "[$count]\n" if $count % 10000 ==0; 
	    $t = $1;               # $t = sequence name (defined further down)
	    $isseq && do {
#		  print STDOUT "seqname : $seq_name\n";
		($hit) = grep { &do_test($_,$seq_name) } @tlist;
#		  print STDOUT "hit : $hit\n";
		(defined($hit) && exists($SEQS{$hit})) && do     {## if curr seq is desired and its hash is defined
								      $SEQS{$hit}{SEQ} = $seq;                 ## get the id's sequence
								  };
	    };
	    $seq_name = $t; # quotemeta($t);                        ## quotemeta will backslash all nonstd chars, so $t = $seq_name
	    $seq = '';
	    $isseq = 1;
	    next;                                             ## next moves to sequence line 
	};
	$seq .= $l;                                           ## $seq now gets sequence ( because of next above)
    };
    
    close(F);
    $isseq && do {                                            ## get last sequence
	($hit) = grep { &do_test($_,$seq_name) } @tlist;
	(defined($hit) && exists($SEQS{$hit})) && ($SEQS{$hit}{SEQ} = $seq);
    };
    
}

######################################################################################
######################################################################################
sub output
{

    if ($output){                                             ## if we want one outfile per id infile
	foreach my $Q (@qlist)
	{
	    open(F,">$Q\.fa");                                 ## create one outfile per id infile
	    foreach my $T (@{ $SEQS{$Q}{TYPE} }) {              ## foreach id in this file
		my ($S,$l,$i);
		$S = $SEQS{$T}{SEQ};                            ## $S = sequence
		$l = length($S);
		$i = 0;
		print F ">$T\n"; 
		
		## convert to FASTA
		while($i<$l)  
		{
		    
		    print F substr($S,$i,60) . "\n"; 
		    $i=$i+60;
		}
	    };
	    close(F);
	}
    }
    elsif($sfile)
    {
	foreach my $Q (@qlist)
	{
	    foreach my $T (@{ $SEQS{$Q}{TYPE} }) {
		my ($S,$l,$i);
		$S = $SEQS{$T}{SEQ};                            ## $S = sequence

		$l = length($S);
		$i = 0;
		open(F, ">$T\.fa");
		print F ">$T\n"; 
		while($i<$l)  
		{
		    print F substr($S,$i,60) . "\n"; 
		    $i=$i+60;
		}
		close(F);
	    }
	}

    }
    else{ ## print all sequences to STDOUT (usefull when only one file with ids)
	foreach my $Q (@qlist)
	{
	    foreach my $T (@{ $SEQS{$Q}{TYPE} }) {             ## foreach id in this file
		my ($S,$l,$i);
		$S = $SEQS{$T}{SEQ};                           ## $S = sequence
		$l = length($S);
		$i = 0;
		print STDOUT ">$T\n"; 
		
		while($i<$l)                                   ## convert to FASTA
		{
		    print STDOUT substr($S,$i,60) . "\n"; 
		    $i=$i+60;
		}
	    };
	}
    }                              
}



sub do_test() { ### A trick of pep's to check whether the current id is longer than the hash key or not and then search for the shortest in the longest
    my ($q,$t,$ql,$tl,$o);
    ($q,$t) = @_;
    $ql = length($q);
    $tl = length($t);
    $o = ($ql > $tl) ? $tl : $ql;
    $o = ((lc(substr($q,0,$o)) eq lc(substr($t,0,$o))) ? $q : undef);
#    print "$q ~ $t -> ".(defined($o) ? $o : "undef")." <-\n";
    return $o
}

sub usage
{
    print STDERR "retrieveseqs.pl will take one or more lists of ids and extract thier sequences from a FASTA file\n\n";
    
    print STDERR "USAGE : retrieveseqs.pl [-ios] <FASTA sequence file> <desired IDs, one per line>\n\n";
    print STDERR "\t-i : use when the ids in the id file are EXACTLY identical to those in the FASTA file\n\t\t-o : will create one fasta file for each of the id files\n\t\t-s : will create one fasta file per id\n\n";
    die("*** A minimum of two input files is required : $! ***\n");
}




#==========================Known Bugs=========================================#
# When 2 different seqs have v. similar names, script gets confused	      #
# eg MMSELP-I and MMSELP-II, will only retunr MMSELP-1 use -i to avoid	      #
#=============================================================================#
