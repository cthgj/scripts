#!/usr/bin/perl -w
use strict;
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Utils::Exception qw(throw);
use Bio::SimpleAlign;
use Bio::AlignIO;
use Bio::LocatableSeq;
use Getopt::Long;


my $target_species=$ARGV[0]|| die("Need a target species as ARGV[0] (taxon id)\n");
# 4932 

Bio::EnsEMBL::Registry->load_registry_from_db(
    -host => 'ensembldb.ensembl.org',
    -user => 'anonymous',
    -port => 5306);


my %taxa=(
	  '10090'=>'mouse',
	  '10116'=>'rat',
	  '9913'=>'cow',
	  '9544'=>'macaque',
	  '9598'=>'chimp'
);
my @species=sort(keys(%taxa));

############################
# Read desired ENSEMBL IDs #
############################
my @names;
open(A,"$ARGV[1]")||die(" Need a list of desired ENSEMBL IDs as ARGV[1] : $!\n");
while (<A>) {
    chomp;
    push @names, $_;
}
close(A);

my $member_adaptor = Bio::EnsEMBL::Registry->get_adaptor('Multi', 'compara', 'Member');
my $homology_adaptor = Bio::EnsEMBL::Registry->get_adaptor('Multi', 'compara', 'Homology');

my $transcript_adaptor = Bio::EnsEMBL::Registry->get_adaptor( 'Human', 'Core', 'Transcript' );
my $tot=scalar(@names);
my $c=0;

my %has;
foreach my $name (@names) {
    $c++;
    $has{$name}=0;
    print STDERR "$c of $tot\r";
    my $member = $member_adaptor->fetch_by_source_stable_id('ENSEMBLGENE',$name) || next;
  #  next;
  #  next unless $member;
    my $foo=$homology_adaptor->fetch_all_by_Member_paired_species($member,"Saccharomyces cerevisiae", "ENSEMBL_ORTHOLOGUES");
    print $foo->stable_id(), "\n";
}

