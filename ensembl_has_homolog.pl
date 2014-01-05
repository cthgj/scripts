#!/usr/bin/perl -w
use strict;
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Utils::Exception qw(throw);
use Bio::SimpleAlign;
use Bio::AlignIO;
use Bio::LocatableSeq;
use Getopt::Long;

unless ($ARGV[0]) {
    $0=~/.+?([^\/]+)$/;
    my $p_name=$1;

print <<Eof;
DESCRIPTION:
     This script will check whether a given list of proteins have homologs in a target species.

USAGE:
     $p_name  <target taxid> <list of ensembl gene IDs>

The default taxid is 4932, yeast. 

Eof

exit();
}

my $target_species=$ARGV[0]|| die("Need a target species as ARGV[0] (taxon id)\n");
# 4932 

Bio::EnsEMBL::Registry->load_registry_from_db(
    -host => 'ensembldb.ensembl.org',
    -user => 'anonymous',
    -port => 5306);



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
    my $homologies = $homology_adaptor->fetch_all_by_Member($member);
    my %foo;
    my %found;
    my @seqs;
    h:foreach my $homology (@{$homologies}) {
	foreach my $m (@{$homology->get_all_Members}) {
	    if ($m->taxon_id eq $target_species) {
		$has{$name}=1;
		next h;
	    } 
	}
    }
}

map{print "$_\t$has{$_}\n"}keys(%has);
