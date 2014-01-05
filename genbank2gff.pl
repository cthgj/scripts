#!/usr/local/bin/perl -w

use strict;

use Bio::SeqIO;

my $usage = "genbank2gffexons.pl\n";

my $in  = Bio::SeqIO->new(-fh => \*STDIN,
                          -format => "genbank");
#my $out = Bio::Tools::GFF->new(-fh => \*STDOUT,
#                               -gff_version => 2);

my $i = 1;
while ( my $seq = $in->next_seq() ) {
  my $seqfullname = "gb|".$seq->accession_number()."|".$seq->display_id()."|".$seq->length()."|bp";
  my $seqname = $seq->display_id();

  print STDERR "$seqname\t$seqfullname\n";
  if ($i > 1) {
    print STDOUT "#\$\n";
  }
  print STDOUT "$seqname\thuman\tSequence\t1\t",$seq->length(),"\t.\t+\t.\t$seqname\n";
  for my $feature ($seq->get_SeqFeatures()) {
    if ($feature->location->isa('Bio::Location::SplitLocationI') &&
        ($feature->primary_tag eq 'mRNA' || $feature->primary_tag eq 'CDS')) {
      my $strand = $feature->strand < 0 ? '-' : '+';
      foreach my $location ($feature->location->sub_Location()) {
        print STDOUT "$seqname\thuman\t",$feature->primary_tag,"\t",
                      $location->start,"\t",$location->end,"\t.\t$strand\t.\t$seqname\n";
      }
    }

  }

  $i = $i + 1;
}
