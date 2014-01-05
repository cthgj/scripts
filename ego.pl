#!/usr/bin/perl -w

#
#Create tbl file with cluster name, seq id and seq from the EGO sequences.
#

use strict;

my @lines;
my $name;
my $seq;

open (FILE, "/home/ug/cchapple/research/seq/EGO/ego3_020103_seq.tbl")|| die "cannot open file $_";

 line: while (<FILE>)
{ print STDERR "@\n";
  
  m/^(.*?)\s(\w*)\n/o;
  
  $name = $1;
  $seq = $2;

  open (CLUS, "/home/ug/cchapple/research/seq/EGO/ego3_020103_orth_clusters")|| die "cannot open file $_";
  
  @lines = <CLUS>;
  foreach my $line (@lines)
  {
      print STDERR "#";
      
      if ($line =~ m/$name/)
      {
	  print STDERR "\n¿\n";
	  $line =~ /^(.*?)\s/go;
	  my $cluster = $1;
	  
	  print STDOUT "$cluster $name $seq\n"; print STDERR "?\n";
	  
	  $name = 0;
	  $seq = 0;
	  
	  next line; 
      }   
  }
  close (CLUS);
  
}

close (FILE);
