#!/usr/bin/perl -w
use strict;
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Utils::Exception qw(throw);
use Bio::SimpleAlign;
use Bio::AlignIO;
use Bio::LocatableSeq;
use Getopt::Long;

Bio::EnsEMBL::Registry->load_registry_from_db(
    -host => 'ensembldb.ensembl.org',
    -user => 'anonymous',
    -port => 5306);

my $member_adaptor = Bio::EnsEMBL::Registry->get_adaptor('Multi', 'compara', 'Member');
my $member = $member_adaptor->fetch_by_source_stable_id('ENSEMBLGENE','ENSG00000004059');

my $family_adaptor = Bio::EnsEMBL::Registry->get_adaptor('Multi','compara','Family');
my $families = $family_adaptor->fetch_all_by_Member($member);

my $foo=$family_adaptor->get_Member_by_source_taxon('ENSP00000269305',9606);
print "$foo\n";die();


foreach my $family (@{$families}) {
    print join(" ", map { $family->$_ }  qw(description description_score))."\n";

    foreach my $member (@{$family->get_all_Members}) {
        print $member->stable_id," ",$member->taxon_id,"\n";
    }

    my $simple_align = $family->get_SimpleAlign();
    my $alignIO = Bio::AlignIO->newFh(
        -interleaved => 0,
        -fh          => \*STDOUT,
        -format      => "phylip",
        -idlength    => 20);

    print $alignIO $simple_align;
}
