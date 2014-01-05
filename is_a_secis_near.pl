#!/usr/local/bin/perl -w

## this script is very crude, make sure it works


use strict;
use Term::ANSIColor; 
use Getopt::Std;

my %opts;
getopts('ghd:D:o',\%opts) || die("no opts\n");
my $GFF = $ARGV[0] || die();
my $orfs = $opts{o} || undef;
my $geneid_gff = $opts{g} || undef;
my $distance = $opts{d} || 1000;
my $distance1 = $opts{D} || 2500;

my (@seqs, @lines );
my (%genes_m, %genes_p, %gffminus, %gffplus, %secises_p, %secises_m, %hassecis_m, %hassecis_p, %secisfound);



&usage() if $opts{h};

# If input file is a geneid ORFs (-zZ, only ORF lines kept) output file:
if ($orfs)
{
    open(GFF, "sort -nk 4 $GFF |") || die("no gff : $!\n");
    while(<GFF>)
    {
	s/^\s+//;
	my ($NAME, $type, $gff_start, $gff_end, $score, $strand, @rest) =  split(/\s+/, $_); 
	push @{$gffminus{$NAME}}, $_ if $strand eq '-';
	push @{$gffplus{$NAME}}, $_ if $strand eq '+';

    }

}
elsif($geneid_gff)
{
    open(GFF, "sort -nk 4 $GFF |") || die("no gff : $!\n");
    while(<GFF>)
    {
	s/^\s+//;
	my ($type, $gff_start, $gff_end, $score, $strand, @rest) =  split(/\s+/, $_);
	my $NAME1 = $rest[scalar(@rest)-1];
	my $NAME = $NAME1;
	$NAME  =~ s/_\d+$//;
	my $a = $NAME . "\t" . $type . "\t" .  $gff_start . "\t" .  $gff_end . "\t" .  $score . "\t" .  $strand . "\t" .  "$NAME1\n";
	push @{$gffminus{$NAME}}, $a if $strand eq '-';
	push @{$gffplus{$NAME}}, $a if $strand eq '+';
    }
}
else
{
    open(GFF, "sort -nk 4 $GFF |") || die("no gff : $!\n");
    while(<GFF>)
    {
	s/^\s+//;
	next if /^\#/;
	my ($NAME, $source, $type, $gff_start, $gff_end, $score, $strand, @rest) = split(/\s+/, $_); 
	my $NAME1 = $rest[scalar(@rest)-1];
	my $a = $NAME . "\t" . $source . "\t" .$type . "\t" .  $gff_start . "\t" .  $gff_end . "\t" .  $score . "\t" .  $strand . "\t" .  "$NAME1\n";
	push @{$gffminus{$NAME}}, $a if $strand eq '-';
	push @{$gffplus{$NAME}}, $a if $strand eq '+';
# 	push @{$gffminus{$NAME}}, $NAME, $type, $gff_start, $gff_end, $score, $strand, @rest if $strand eq '-';
# 	push @{$gffplus{$NAME}}, $NAME, $type, $gff_start, $gff_end, $score, $strand, @rest  if $strand eq '+';
    }
}

close(GFF);
my @minus = keys(%gffminus);
my @plus = keys (%gffplus);

open(SECISFILE, "<$ARGV[1]") || die();
while(<SECISFILE>)
{
    next unless />/;
    
    />(.*?)\s+\[(\d+)\s+-\s+(\d+)\]/ || die("The second argument must be a SECISearch outfile : $!\n");
    my $ll = $1 . " " . $2 . " " . $3;
    if ($3 < $2) { push @{$secises_m{$1}}, $ll }
    else { push @{$secises_p{$1}}, $ll }
}
close(SECISFILE);


foreach my $seq_name (@plus)
{
  gff:foreach my $gff_line (@{$gffplus{$seq_name}}) 
  {
    
      my ($NAME, $source, $type, $gff_start, $gff_end, $score, $strand, @rest) =  split(/\s+/, $gff_line); 
      secis:foreach my $secis (@{$secises_p{$NAME}})
      {
	  my ($secis_name, $secis_start, $secis_end) = split(/\s+/,$secis);
	  
	  next secis if $secis_start < $gff_end;
	  if(($type eq "Terminal")  || ($type eq "Single"))
	  {
	      
	      if ($secis_start - $gff_end > $distance)
	      {
		  next secis;
	      }
	      else
	      {
		  
		  push @{$hassecis_p{$NAME}}, $gff_line ;
		  my $k = $NAME . " " . $gff_start . " " . $gff_end . " " . $strand;
		  my $distance_from_STOP = $secis_start - $gff_end;
		  $secisfound{$k} = $secis_start . " - " .  $secis_end . "( $distance_from_STOP from STOP )";
		  next gff;
	      }
	  }
	  else
	  {
	      if ($secis_start - $gff_end > $distance1)
	      {#print "aaaaaa $secis_name, $secis_start, $secis_end :: $gff_end -- $distance1\n";
		  next secis;
	      }
	      else
	      {		  
		  push @{$hassecis_p{$NAME}}, $gff_line ;
		  my $k = $NAME . " " . $gff_start . " " . $gff_end . " " . $strand;
		  my $distance_from_STOP = $secis_start - $gff_end;
		  $secisfound{$k} = $secis_start . " - " .  $secis_end . "( $distance_from_STOP from STOP )";
		  next gff;
	      }
	  }  
      }
  }
}
foreach my $seq_name (@minus)
{
  gff:foreach my $gff_line (@{$gffminus{$seq_name}}) 
  {
      my ($NAME, $source, $type, $gff_start, $gff_end, $score, $strand, @rest) =  split(/\s+/, $gff_line); 
      secis:foreach my $secis (@{$secises_m{$NAME}})
      {
	  my ($secis_name, $secis_start, $secis_end) = split(/\s+/,$secis);
	  next secis if $secis_end > $gff_start;
	  if(($type eq "Terminal")  || ($type eq "Single"))
	  {
	      if ($gff_start - $secis_end > $distance)
	      {
		  next secis;
	      }
	      else
	      {
		  
		  push @{$hassecis_p{$NAME}}, $gff_line ;
		  my $k = $NAME . " " . $gff_start . " " . $gff_end . " " . $strand;
		  my $distance_from_STOP = $gff_start - $secis_end;
		  $secisfound{$k} = $secis_start . " - " .  $secis_end . "( $distance_from_STOP from STOP )";
		  next gff;
	      }
	  }
	  else
	  {
	      if ($gff_start - $secis_end > $distance1)
	      {
		  next secis;
	      }
	      else
	      {
		  push @{$hassecis_p{$NAME}}, $gff_line ;
		  my $k = $NAME . " " . $gff_start . " " . $gff_end . " " . $strand;
		  my $distance_from_STOP = $gff_start - $secis_end;
		  $secisfound{$k} = $secis_start . " - " .  $secis_end . "( $distance_from_STOP from STOP )";
		  next gff;
	      }
	  }
      }  
  }
}
my @names = keys(%hassecis_p);
map
{
    foreach my $ss (@{$hassecis_p{$_}})
    {
	if($geneid_gff)
	{
	    my ($NAME, $type, $gff_start, $gff_end, $score, $strand, $NAME1) =  split(/\s+/, $ss);
	    print STDOUT "$NAME1\n"
	}
	elsif($orfs)
	{
	    my ($NAME, $type, $gff_start, $gff_end, $score, $strand, @rest) =  split(/\s+/, $ss); 
	    my $ll = $NAME . " " . $gff_start . " " . $gff_end . " " . $strand;
	    my $tmp = @rest;
	    print ">$NAME" . "_" . $gff_start . "_" . $gff_end . "_" . $strand . "\n";
	    my $l = length($rest[$tmp-1]);
	    my $ii = 0;
	    while($ii<$l)                                   ## convert to FASTA
	    {
		print substr($rest[$tmp-1],$ii,60) . "\n"; 
		$ii=$ii+60;
	    }
	}
	else
	{
	    my ($NAME, $source, $type, $gff_start, $gff_end, $score, $strand, @rest) =  split(/\s+/, $ss); 
	    my $k = $NAME . " " . $gff_start . " " . $gff_end . " " . $strand;
	    print "$NAME $secisfound{$k}\n";	    
	}
	
    } 
}


@names;
@names = keys(%hassecis_m);
map
{
    foreach my $ss (@{$hassecis_m{$_}})
    {
	if($geneid_gff)
	{
	    my ($NAME, $type, $gff_start, $gff_end, $score, $strand, $NAME1) =  split(/\s+/, $ss);
	    print STDOUT "$NAME1\n"
	}
	elsif($orfs)
	{	
	    my ($NAME, $type, $gff_start, $gff_end, $score, $strand, @rest);
	    my $ll = $NAME . " " . $gff_start . " " . $gff_end . " " . $strand;
	    my $tmp = @rest;
	    print ">$NAME" . "_" . $gff_start . "_" . $gff_end . "_" . $strand . "\n";
	    my $l = length($rest[$tmp-1]);
	    my $ii = 0;
	    while($ii<$l)                                   ## convert to FASTA
	    {
		print substr($rest[$tmp-1],$ii,60) . "\n"; 
		$ii=$ii+60;
	    }
	}
	else
	{
	    my ($NAME, $source, $type, $gff_start, $gff_end, $score, $strand, @rest) =  split(/\s+/, $ss); 
	    my $k = $NAME . " " . $gff_start . " " . $gff_end . " " . $strand;
	    print "$NAME $secisfound{$k}\n";	    
	}
    }
}
@names;



sub usage
{
    open(HELP, "| more") ;
    print HELP <<EndOfHelp;
  PROGRAM: 
    is_a_secis_near.pl
	
  USAGE:  
	is_a_secis_near.pl [options] <feature file> <.secis file>

is_a_secis_near.pl will take a set of features (gff, geneid orf, or geneid gff)
and a SECISearch out file and return those features which are within 
a given distance from a SECIS element.

COMMAND-LINE OPTIONS:

    -d : Maximum distance between feature and secis allowed (for terminal/single exons)
         Default : 1000
    -D : Maximum distance between feature and secis allowed (for other exons)
         Default : 2500
    -h : Show this help and exit
    -g : Geneid gff file 
    -o : Geneid orfs file
    

EndOfHelp
close(HELP);

exit(1);
}
