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
getopts('hVvuoisndf',\%opts);

&usage unless ($ARGV[0] && $ARGV[1]);

my %seqname = ();
my %SEQS    = ();
my %seqlist = ();
my @tlist   = ();
my @notfound;
my $output = $opts{o} || undef;
my $identical = $opts{i} || 0;
my $sfile = $opts{s}|| undef;
my $names = $opts{n}||undef;
my $verbose = $opts{v}|| undef;
my $uniprot = $opts{u}|| undef;
my $V = $opts{V}|| undef;
my %found;
my $filename;
my $got=0;
my $count;
my ($SEQFILE,@qtfile) = @ARGV unless $names;
my $debug = $opts{d} || undef;
my $fast = $opts{f} || undef;
&usage() if $opts{h};


$verbose = 1 if $V;

if($names) ## if names are given as cmdline arguments
{
    $SEQFILE = $ARGV[0];
    @qtfile = split(/\s+/, $ARGV[1]);
    my $id; 
    foreach my $name(@qtfile)
    {
	if($fast)
	{
	    $name =~  /^([^\s]*)/;
	    $id = $1;
	}
	else
	{
	    $id = $name;
	}
	## If it doesn't exist, create "%{ $SEQS{$id}" foreach id with 2 keys
	## one a list of the sequence ids (list) the other the seq
#	defined(%{ $SEQS{$id} }) ||                           
	defined(${ $SEQS{$id}{TYPE} }) ||                           
	    (%{ $SEQS{$id} } = ( TYPE => [], SEQ => '' ));   
	# next;
    }
    push @{ $SEQS{$id}{TYPE} }, $id;
    
}
else ## if file of names
{
    if ((my $k = @qtfile) > 1) 
    {
 	foreach my $file (@qtfile) {  # foreach file with list of names
	    if ($file =~ /\.gz$/) 
	    {
		open(F,"zcat $file |") || &usage; 
	    } 
	    else 
	    {
		open(F,"< $file")|| &usage;
	    };
	    my $query;
	    while (<F>) 
	    {
		my $id; 
		chomp;
		if($fast)
		{
		    /^([^\s]*)/;
		    $id = $1;
		}
		else
		{
		    ($id = $_) =~ s/\s*$//og;    # remove trailing space
		}
		
		push @tlist, $id;
		defined(${ $SEQS{$id}{TYPE} }) ||                           ## If doesn't exist, create "%{ $SEQS{$id}" foreach id with 2 keys
		    (%{ $SEQS{$id} } = ( TYPE => [], SEQ => '', FNAME => $file ));    ## one a list of the sequence ids (list) the other the seq
		$. == 1 && do {                                    
		    $query = $id;  # $query= FIRST id in file
		    
		};
		push @{ $SEQS{$file}{TYPE} }, $id;
	    };
	    close(F);
	};
    }
    else{
	$filename = $qtfile[0];
	if ($filename =~ /\.gz$/) {
		open(F,"zcat $filename |") || &usage; 
	    } else {
		open(F,"< $filename")|| &usage;
	    };
	    my $query;
	    while (<F>) { 
		
		my $id; 
		chomp;	
		if($fast)
		{
		    /^([^\s]*)/;
		    $id = $1;
		    
		}
		else
		{
		    ($id = $_) =~ s/\s*$//og;
		}
		push @tlist, $id;
#		print STDERR "$id\n";
		defined(${ $SEQS{$id}{TYPE} }) ||                           ## If doesn't exist, create "%{ $SEQS{$id}" foreach id with 2 keys
		    (%{ $SEQS{$id} } = ( TYPE => [], SEQ => ''));    ## one a list of the sequence ids (list) the other the seq
		$. == 1 && do {                                    
		    $query = $id;  # $query= FIRST id in file
		};
	    }
    }
}

@tlist = @qtfile if $names;
my @list_of_names = @tlist; ## so that I can pop found names from tlist and then retrieve all anmes again
my $number_of_names = @tlist;
$number_of_names = @qtfile if $names;
#&debug("tlist : @tlist");

if($fast)
{
    &fast_seqs(); 
    print STDERR "@notfound were not found\n" if $notfound[0];
    print STDERR "($got of $number_of_names found)\n";
    exit(0);
}
else{
    $identical == 1 ? &get_identical_seqs() : &get_seqs();

}
&output();
print STDERR "@notfound were not found\n" if $notfound[0];


#################################################################
#################################################################
sub fast_seqs
{
    
    my ($seq_name, $seq, $isseq, $hit);
    $isseq = 0;
    if ($SEQFILE =~ /\.gz$/) {
	open(F,"zcat $SEQFILE |");
    } else {
	open(F,"< $SEQFILE") || die("cannot open $SEQFILE : $!\n");
    };
    
    $count=1;
    my $want=0;
    while(<F>) 
    {
	if(/>/ && $got == $number_of_names){
	    print STDERR "($got of $number_of_names found)\n" if $verbose;
	    exit(0);
	    }
	my ($l,$t);
	$t="bib";
	next if /^\s*$/o;
	chomp;
	$l = $_;

	if($l =~ /^>([^\s]*)\s?/o)
	{
	    $t = $1; 
	    $want=0;
	    $count++;
	    print STDERR "." if $count % 100==0 && $verbose;;
 	    print STDERR "[$count ($got found)]\n" if $count % 5000 ==0 && $verbose;	    
	    if ($uniprot){
		$l=~ s/^>..\|//;
		$l=~ /^(.+?)\s/ ||die("crap!\n");
		my @kk=split(/\|/,$1);
	      koko:foreach my $t (@kk){
		  next koko if $want==1;
		  if(defined($SEQS{$t}))
		  {
		      $found{$t}++;
		      $got++ if $found{$t}==1;
		      $want=1;
		      print "$_\n";
		      print "$t : $want\n";
		  }
		  else{
		      $want=0;
		      print "$t : $want\n";
		  }   
	      }
		$t = $1;
	    }
	    else{
		if(defined($SEQS{$t}))
		{
		    $found{$t}++;
		    $got++ if $found{$t}==1;
		    $want=1;
		    print "$_\n";
		    next;
		}
		else{$want=0;}
	    }
	}
	else
	{
	    
	    if($want==1)
	    {
		print "$_\n";
	    }
	    else
	    {
		next;
	    }
	}
	     
    }
}
#################################################################
#################################################################
sub get_seqs
{

    my ($seq_name, $seq, $isseq, $hit);
    $isseq = 0;
    if ($SEQFILE =~ /\.gz$/) {
	open(F,"zcat $SEQFILE |");
    } else {
	open(F,"< $SEQFILE") || die("cannot open $SEQFILE : $!\n");
    };
    
    $count=1;
    while(<F>) 
    {
	
	if ($got ==  $number_of_names){
	    &output();
	    exit(0);
	};
	my ($l,$t);
	next if /^\s*$/o;
	chomp;
	$l = $_;
	$l =~ /^>((.*?)(\s+.*)?)$/o && do 
	{
	    @tlist = ();
	    map{push @tlist, $_ unless $found{$_}} @list_of_names;
	    $count++;
	    print STDERR "." if $count % 100 ==0 && $verbose; 
	    print STDERR "[$count]\n" if $count % 5000 ==0 && $verbose; 
	    $t = $2;               # $t = sequence name (defined further down)
	    $isseq && do 
	    {
		($hit) = ( grep { &do_test($_,$seq_name)} @tlist);
			(defined($hit) && exists($SEQS{$hit})) && do     
		{## if curr seq is desired and its hash is defined
		     $SEQS{$hit}{SEQ} = $seq;                 ## get the id's sequence
		     $found{$hit}++;
		     $got++;
		     $seqname{$hit} = $seq_name;
		     print STDERR "!" if $V;
		     
		 };
	    };
	    $seq_name = $t; 
	    $seq = '';
	    $isseq = 1;
	    next;                                             ## next moves to sequence line 
	};
	$seq .= $l;                                           ## $seq now gets sequence ( because of next above)
    };
    
    close(F);
    $isseq && do                                             ## get last sequence
    {
	($hit) = grep { &do_test($_,$seq_name) } @tlist;
	(defined($hit) && exists($SEQS{$hit})) && do     
		{## if curr seq is desired and its hash is defined
		     $SEQS{$hit}{SEQ} = $seq;                 ## get the id's sequence
		     $found{$hit}++;
		     $got++;
		     $seqname{$hit} = $seq_name;
		     print STDERR "!" if $V;
		     
		 };
    };

}
#################################################################
#################################################################
sub get_identical_seqs
{
    my ($seq_name, $seq, $isseq, $hit);
    $isseq = 0;
    if ($SEQFILE =~ /\.gz$/) {
	open(FILE,"zcat $SEQFILE |");
    } else {
	open(FILE,"< $SEQFILE");
    };
    
     $count=1;
    
    while(<FILE>) {

	&output() if $got ==  $number_of_names;
	my ($l,$t);
	next if /^\s*$/o;
	chomp;
	$l = $_;
	$l =~ /^>((.*?)(\s+.*)?)$/o && do {
	    next unless exists($SEQS{$1});
	    @tlist = ();
	    map{push @tlist, $_ unless $found{$_}} @list_of_names;
	    $count++;
	    print STDERR "." if $count % 100 ==0 && $verbose; 
	    print STDERR "[$count]\n" if $count % 5000 ==0 && $verbose;
	    $t = $1;               # $t = sequence name (defined further down)
	    $isseq && do {	
		exists($SEQS{$seq_name}) && do { 
		    $SEQS{$seq_name}{SEQ} = $seq; 
		    $seqname{$seq_name} = $seq_name;
		    print STDERR "!" if $V;
		    $got++; 
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
    }; print STDERR "[$count]\n"; 
}
######################################################################################
######################################################################################
sub output
{

    print STDERR "[$count ($got/$number_of_names found]\n";   
    if ($output)                                 ## if we want one outfile per id infile
    {              
	foreach my $file (@qtfile)
	{
	    if($names)
	    {
		foreach my $Q (@qtfile)
		{
		    open(F,">$Q.fa");
		    my ($S,$l,$i);
		    $S = $SEQS{$Q}{SEQ};                           ## $S = sequence
		    $l = length($S);
		    $i = 0;
		    print F ">$Q\n";
		    while($i<$l)                                   ## convert to FASTA
		    {
			print F substr($S,$i,60) . "\n"; 
			$i=$i+60;
		    }
		}
	    }
	    else
	    {
		open(F,">$file\.fa");                                  ## create one outfile per id infile
		foreach my $name (@{ $SEQS{$file}{TYPE} }) {              ## foreach id in this file
		    my ($S,$l,$i);
		    $S = $SEQS{$name}{SEQ};                               ## $S = sequence
		    $l = length($S);
		    $i = 0;
		    if($SEQS{$name}{SEQ}) {  print F ">$seqname{$name}\n";}	
		    
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
    }
    elsif($sfile)
    {
	foreach my $name (@list_of_names)
	{
	    $name =~ /^(.*?)\s/;
	    my $q = $1 || $name;
	    open(F, ">$q.fa");
	    my ($S,$l,$i);
	    $S = $SEQS{$name}{SEQ};                            ## $S = sequence
	    $l = length($S);
	    $i = 0;
	    if($SEQS{$name}{SEQ}) {  print F ">$seqname{$name}\n";}
	    while($i<$l)  
	    {
		print F substr($S,$i,60) . "\n"; 
		$i=$i+60;
	    }
	    
	    close(F);
	}
	
    }
    elsif($names)
    { 
	foreach my $Q (@qtfile)
	{
	    my ($S,$l,$i);
	    $S = $SEQS{$Q}{SEQ}; ## $S = sequence
	    
	    $l = length($S);
	    $i = 0;
	    print STDOUT ">$seqname{$Q}\n";
	    while($i<$l)                                   ## convert to FASTA
	    {
		print STDOUT substr($S,$i,60) . "\n"; 
		$i=$i+60;
	    }
	}
    }
    else{ ## print all sequences to STDOUT (usefull when only one file with ids)
		foreach my $name (@list_of_names)
		{
		    
		    my ($S,$l,$i);
		    $S = $SEQS{$name}{SEQ};                           ## $S = sequence
		    $l = length($S);
		    $i = 0;
		    if($l>0) { print STDOUT ">$seqname{$name}\n";}
		    else{push @notfound, $name;  }
		    
		    while($i<$l)                                   ## convert to FASTA
		    {
			print STDOUT substr($S,$i,60) . "\n"; 
			$i=$i+60;
		    }
		}
	}
}



sub do_test() { ### A trick of pep's to check whether the current id is longer than the hash key or not and then search for the shortest in the longest
    my ($id_passed,$seqname,$ll,$seqnamel,$o);

    ($id_passed,$seqname) = @_;
    $ll = length($id_passed);
    $seqnamel = length($seqname);
    my ($target,$query) = ($ll > $seqnamel) ? ($id_passed,$seqname) : ($seqname,$id_passed);
      $o = (($target =~ /$query/i) ? $seqname : undef);
    &debug("id : $id_passed, seq : $seqname, l : $ll, seqn :$seqnamel");
    &debug("o:$o") if $o;
      return $o
}


sub debug
{
    if ($debug)
    {
	print STDERR "@_\n";
    }
}


sub usage
{
    print STDERR <<EndOfHelp;

    retrieveseqs.pl will take one or more lists of ids and extract their sequences from 
    multi FASTA file
    
    USAGE : retrieveseqs.pl [-viofsn] <FASTA sequence file> <desired IDs, one per line>

    -v : verbose output, print a progress indicator (a "." for every 1000 sequences processed)
    -V : as above but a "!" for every desired sequence found.
    -f : fast, takes first characters of name "(/^([^\\s]*)/)" given until the first space as the search string
         make SURE that those chars are UNIQUE.
    -i : use when the ids in the id file are EXACTLY identical
         to those in the FASTA file
    -h : Show this help and exit.
    -o : will create one fasta file for each of the id files
    -s : will create one fasta file per id
    -n : means that the last arguments (after the sequence file)
         passed are a QUOTED list of the names desired.
    -u : assumes uniprot format ids (separated by |)
EndOfHelp

die("\n*** A minimum of two input files is required : $! ***\n\n");
}


#==========================Known Bugs=========================================#
# When 2 different seqs have v. similar names, script gets confused	      #
# eg MMSELP-I and MMSELP-II, will only return MMSELP-I use -i to avoid	      #
#=============================================================================#


