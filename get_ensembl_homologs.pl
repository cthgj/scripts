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
open(A,"$ARGV[0]")||die(" Need a list of desired ENSEMBL IDs as ARGV[0] : $!\n");
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
foreach my $name (@names) {
    $c++;
    print STDERR "$c of $tot\r";
    my $member = $member_adaptor->fetch_by_source_stable_id('ENSEMBLGENE',$name) || next;
  #  next;
  #  next unless $member;
    my $homologies = $homology_adaptor->fetch_all_by_Member($member);
    my %foo;
    my %found;
    my @seqs;
    foreach my $homology (@{$homologies}) {
	foreach my $m (@{$homology->get_all_Members}) {
	    next unless defined($taxa{$m->taxon_id});
	    ## save each orth by species in %foo
	    $foo{$m->taxon_id}=$m->stable_id;
	    $found{$m->taxon_id}++;
	    my $s=$m->get_canonical_Member;
	    push @seqs, "$taxa{$m->taxon_id} ", $s->stable_id,"\t", $s->sequence, "\n";
	}  
    }
    next unless scalar(keys(%found))==scalar(keys(%taxa));
    open(my $fh, "|-", "TblToFasta > $name.pep");
    print $fh @seqs;
    close($fh);
    # print "$name\t";
    # foreach my $sp (@species) {
    # 	$foo{$sp}||='-';
    # 	print "\t$foo{$sp}";
    # }
    # print "\n";

}


    ##################
    # Get alignments #
    ##################
 #    foreach my $homology (@{$homologies}) {
# 	my $simple_align = $homology->get_SimpleAlign();
# 	my $alignIO = Bio::AlignIO->newFh
# 	    (
# 	     -interleaved => 0,
# 	     -fh => \*STDOUT,
# 	     -format => "clustalw",
# 	     -idlength => 20
# 	    );

# print $alignIO $simple_align;
# die();


    # foreach my $this_homology (@$homologies) {
   #     my $homologue_genes = $this_homology->gene_list();
   #     my ($id1,$id2)=($$homologue_genes[0]->stable_id,$$homologue_genes[1]->stable_id);
   #     print "$id1 : $id2\n";
   #     # print join(" and ", @$homologue_genes), " are ",
   #     # 	   $this_homology->description, "\n";
   # }
