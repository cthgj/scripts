#!/usr/local/bin/perl -w


# This script takes as input a geneid ( in order to have the name of the predicted protein) GFF file and a SECISearch .secis 
# file (refering to the same FASTA seq as the GFF) and returns the gene from the GFF file closest to each of the SECISes.  
# REMEMBER to only give terminal exons in the gff file, otherwise it is not very useful
#
# The -c option gives colorfull output, -o is for geneid orfs file (prob a once off option)
#
#  USAGE : positionsecis.pl <GFF> <SECISearch.secis>
use strict;
use Term::ANSIColor; 
use Getopt::Std;

my (@seqs, @lines, @secises_p, @secises_m);
my (%genes_m, %genes_p, %gffminus, %gffplus);

my $laststart;
my($a, $type, $gff_start, $gff_end, $score, $strand, @rest, $kk,$secis_start, $STOP, $START, $name);
my $NAME = 'j';
my @haha;
## Load gff lines into %gffplus and %gffminus depending on strand
my $lastgffstart = 0;
my $found;

my %opts;
getopts('co',\%opts) || die("no opts\n");
my $GFF = $ARGV[0] || die();
my $color = $opts{c} || undef;
my $orfs = $opts{o} || undef;

if ($orfs)
{
    open(GFF, "sort -nk 4 $GFF |") || die("no gff : $!\n");
    while(<GFF>)
    {

	($NAME, $type, $gff_start, $gff_end, $score, $strand, @rest) =  split(/\s+/, $_); 
	push @{$gffminus{$NAME}}, $_ if $strand eq '-';
	push @{$gffplus{$NAME}}, $_ if $strand eq '+';

    }

}
else
{
    open(GFF, "sort -nk 4 $GFF |") || die("no gff : $!\n");
    while(<GFF>)
    {
	next if /^\#/;
	($NAME, $a, $type, $gff_start, $gff_end, $score, $strand, @rest) = split(/\t/, $_); 
	push @{$gffminus{$NAME}}, $_ if $strand eq '-';
	push @{$gffplus{$NAME}}, $_ if $strand eq '+';

    }
}

my @minus = keys(%gffminus);
my @plus = keys (%gffplus);




## Read SECIS outfile and push SECIS name, start and end into @secises
open(F, "<$ARGV[1]") || die();
while(<F>)
{
    next unless />/;
    />(.*?)\s+\[(\d+)\s+-\s+(\d+)\]/ || die("The second argument must be a SECISearch outfile : $!\n");
    my $ll = $1 . " " . $2 . " " . $3;
    if ($3 < $2) { push @secises_m, $ll }
    else { push @secises_p, $ll }
}
close(F);


## Find the gene closest to each - strand SECIS element
foreach my $scaf_name (@minus)
{ 
#    print STDOUT "new scaffold\n";
#    my @gffs = @{$gffminus{$scaf_name}};
    $laststart = 0;
   
   	secis:foreach my $secis (@secises_m)
	{  
	   ($name, $START, $STOP) = split (/\s+/, $secis);
	   $found = 0;
	   gff:foreach my $gff_line (@{$gffminus{$scaf_name}})
	   {
	      
	       if ($orfs)
	       {
		      
		   ($NAME, $type, $gff_start, $gff_end, $score, $strand, @rest) =  split(/\t/, $gff_line) ;
		   chomp $rest[1]; 
		   if($name eq $NAME)
		   {
		       $color ? &get_minus($name, $NAME, $START, $STOP, $rest[1], $gff_end, $gff_start) : &get_minus_nocolor($name, $NAME, $START, $STOP, $rest[1], $gff_end, $gff_start);
		   }
		   
	       }
	       else
	       {
		   ($NAME, $a, $type, $gff_start, $gff_end, $score, $strand, @rest) = split(/\t/, $gff_line) ;
		   chomp $rest[1]; 
		   if($name eq $NAME)
		   {
		       $color ? &get_minus($name, $NAME, $START, $STOP, $rest[1], $gff_end, $gff_start) : &get_minus_nocolor($name, $NAME, $START, $STOP, $rest[1], $gff_end, $gff_start);
		   }
	       }
	   }
       }
}
# foreach my $gff_line (@{$gffminus{$scaf_name}})
#     {print STDOUT "new gff line\n";
# 	($NAME, $a, $type, $gff_start, $gff_end, $score, $strand, @rest) = split(/\t/, $gff_line) ;
# 	chomp $rest[1]; 
	
# 	foreach my $secis (@secises_m)
# 	{  print STDOUT "new SECIS\n";
# 	    ($name, $START, $STOP) = split (/\s+/, $secis);
# 	   if($name eq $NAME)
# 	   {
# 	       print STDOUT "*********$name == $NAME**********\n";
# 	       &get_minus($name, $NAME, $START, $STOP, $rest[1], $gff_en, $gff_start); 
# 	   }
# 	   else { print STDOUT "did next\n"; next}
# #	    print "did next\n" && next unless $name eq $NAME;
	   
# 	}
#     } 
# }
foreach my $scaf_name (@plus)
{
    $laststart = 0;
  secis:foreach my $secis (@secises_p)
  {
      $found = 0;
      # my $secis_start;
      ($name, $START, $STOP) = split (/\s+/, $secis);
      foreach my $gff_line (@{$gffplus{$scaf_name}})
      {
	   
	  if ($orfs)
	  {	     
	      ($NAME, $type, $gff_start, $gff_end, $score, $strand, @rest) =  split(/\s+/, $gff_line) ;
	      chomp $rest[1]; 

	      if($name eq $NAME)
	      {
		  
		  $color ? &get_plus($name, $NAME, $START, $STOP, $rest[1], $gff_end, $gff_start) : &get_plus_nocolor($name, $NAME, $START, $STOP, $rest[1], $gff_end, $gff_start);
	      }
	      
	  }
	  else
	  {
	      ($NAME, $a, $type, $gff_start, $gff_end, $score, $strand, @rest) = split(/\t/, $gff_line) ;
	      chomp $rest[1]; 
	      
	      if($name eq $NAME)
	      {
		  $color ? &get_plus($name, $NAME, $START, $STOP, $rest[1], $gff_end, $gff_start) : &get_plus_nocolor($name, $NAME, $START, $STOP, $rest[1], $gff_end, $gff_start);
#	      print STDOUT "xx $NAME $found\n";
	      }
#	  else { print STDOUT "did next\n"; next}  
	  }   
      }
  }
}

my @k = keys(%genes_m);
map {print "$_\t:\t$genes_m{$_}\n"} @k;


my @l= keys(%genes_p);
map {print "$_\t:\t$genes_p{$_}\n"} @l;


############################################################################################################
sub get_minus
{ 
    ($name, $NAME, $START, $STOP, $rest[1], $gff_end, $gff_start) = @_;
   #($secis_start, $STOP, $name) = @_;
    $secis_start = $START;
#    print "$name == $NAME && $secis_start < $gff_start && $gff_start < $laststart \n";
    ## Force it to run for the first time when $laststart(==0) < $gff_start by definition 
#    if($name eq $NAME && $secis_start < $gff_start && $laststart == 0)
#     if($laststart == 0)
#     {
# 	print "one $gff_start $laststart ($strand)\n";
# 	$kk = $gff_start - $secis_start;
# 	$genes_m{$NAME} = $rest[1] . "  -  \t\tend : $gff_end\tsecis : $secis_start\t$kk nt downstream";
# 	$laststart = $gff_start;
# 	if($kk < 1500)
# 	{
# 	    $genes_m{$NAME} =~ s/(.*)/color("green").$1.color("reset")/oge;
# 	}
# 	else
# 	{
# 	    $genes_m{$NAME} =~ s/(.*)/color("red").$1.color("reset")/oge;
# 	}
#     }
#     ## Now, run as normal since $laststart will be >0
#     elsif($name eq $NAME && $secis_start < $gff_start && $gff_start < $laststart) # && $gff_start - $secis_start < 1000)
#     {
# 	print "three " . ($laststart < $gff_start ? " YES":" NO") ." $gff_start $laststart ($strand)\n";
# 	$kk = $gff_start - $secis_start;
# 	$genes_m{$NAME} = $rest[1] . "    \t-\tend : $gff_start\tsecis : $secis_start\t$kk nt downstream";
# 	$laststart = $gff_start;
# 	if($gff_start - $secis_start < 1500)
# 	{
# 	    $genes_m{$NAME} =~ s/(.*)/color("green").$1.color("reset")/oge;
# 	}
# 	else
# 	{
# 	    $genes_m{$NAME} =~ s/(.*)/color("red").$1.color("reset")/oge;
# 	}
# 	$found = 1;
# #	print STDOUT "xxx : $found\n";

#     }
#     else {print "four $laststart < $gff_start" . ($laststart < $gff_start ? " YES":" NO") ."($strand)\n"; $found = 0; $laststart = $gff_start;}
    
    if ($secis_start < $gff_start)
    {
	if($laststart == 0)
	{  
	    $kk = $gff_start - $secis_start;
	    $genes_m{$NAME} = $rest[1] . "  -  \t\tend : $gff_start\tsecis : $secis_start\t$kk nt downstream";
	    $laststart = $gff_start;
	    if($kk < 1500)
	    {
		$genes_m{$NAME} =~ s/(.*)/color("green").$1.color("reset")/oge;
	    }
	    else
	    {
		$genes_m{$NAME} =~ s/(.*)/color("red").$1.color("reset")/oge;
	    }
	}
	if($gff_start < $laststart)
	{  # print STDOUT "*********$name == $NAME**********\n";
	   # print "hih1\n";
	   # print "three " . ($laststart < $gff_start ? " YES":" NO") ." $gff_start $laststart ($strand)\n";
	    $kk = $gff_start - $secis_start;
	    $genes_m{$NAME} = $rest[1] . "    \t-\tend : $gff_start\tsecis : $secis_start\t$kk nt downstream";
	    $laststart = $gff_start;
	    if($gff_start - $secis_start < 1500)
	    {
		$genes_m{$NAME} =~ s/(.*)/color("green").$1.color("reset")/oge;
	    }
	    else
	    {
		$genes_m{$NAME} =~ s/(.*)/color("red").$1.color("reset")/oge;
	    }
	    
	}
	
    }


    #return $found;
}


## Find the gene closest to each + strand SECIS element
sub get_plus
{     
    ($name, $NAME, $START, $STOP, $rest[1], $gff_end, $gff_start) = @_;
    $secis_start = $STOP;
     if ($name eq $NAME && $secis_start > $laststart && $secis_start > $gff_end) # && $secis_start - $gff_end < 1500)
     {	 
	 $kk = $secis_start - $gff_end; 		
	 $genes_p{$NAME} = $rest[1] . "  +  \t\tend : $gff_end\tsecis : $secis_start\t$kk nt downstream";
	 if($secis_start - $gff_end < 1500)
	 {
	     $genes_p{$NAME} =~ s/(.*)/color("green").$1.color("reset")/oge;
	 }
	 else
	 {
	     $genes_p{$NAME} =~ s/(.*)/color("red").$1.color("reset")/oge;
	 }
	 $laststart = $gff_start;
     }
}
	


sub get_minus_nocolor
{ 
    ($name, $NAME, $START, $STOP, $rest[1], $gff_end, $gff_start) = @_;
     $secis_start = $START;
    
    if ($secis_start < $gff_start)
    {
	if($laststart == 0)
	{  
	    $kk = $gff_start - $secis_start;
	    $genes_m{$NAME} = $rest[1] . "  -  \t\tend : $gff_start\tsecis : $secis_start\t$kk nt downstream";
	    $laststart = $gff_start;
	    # if($kk < 1500)
# 	    {
# 		$genes_m{$NAME} =~ s/(.*)/color("green").$1.color("reset")/oge;
# 	    }
# 	    else
# 	    {
# 		$genes_m{$NAME} =~ s/(.*)/color("red").$1.color("reset")/oge;
# 	    }
	}
	if($gff_start < $laststart)
	{  
	    $kk = $gff_start - $secis_start;
	    $genes_m{$NAME} = $rest[1] . "    \t-\tend : $gff_start\tsecis : $secis_start\t$kk nt downstream";
	    $laststart = $gff_start;
	    # if($gff_start - $secis_start < 1500)
# 	    {
# 		$genes_m{$NAME} =~ s/(.*)/color("green").$1.color("reset")/oge;
# 	    }
# 	    else
# 	    {
# 		$genes_m{$NAME} =~ s/(.*)/color("red").$1.color("reset")/oge;
# 	    }
	    
	}
	
    }

}

sub get_plus_nocolor
{   
    ($name, $NAME, $START, $STOP, $rest[1], $gff_end, $gff_start) = @_;
    $secis_start = $STOP;
     if ($name eq $NAME && $secis_start > $laststart && $secis_start > $gff_end) # && $secis_start - $gff_end < 1500)
     {	 
	 $kk = $secis_start - $gff_end; 		
	 $genes_p{$NAME} = $rest[1] . "  +  \t\tend : $gff_end\tsecis : $secis_start\t$kk nt downstream";
	 # if($secis_start - $gff_end < 1500)
# 	 {
# 	     $genes_p{$NAME} =~ s/(.*)/color("green").$1.color("reset")/oge;
# 	 }
# 	 else
# 	 {
# 	     $genes_p{$NAME} =~ s/(.*)/color("red").$1.color("reset")/oge;
# 	 }
# 	 $laststart = $gff_start;
      }
}
