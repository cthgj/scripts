#!/usr/bin/perl -w
use strict;
use lib '/users/rg/mmariotti/libraries/ensembl_api/modules/' ;
#use lib '/home/ug/talioto/myperlmods/ensembl/modules';
use Getopt::Std;
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::DBSQL::DBAdaptor;

my %opts;
getopts('hvt:i:s:g:',\%opts);
&usage() if $opts{'h'};

my $reg_conf="/users/rg/mmariotti/.ensembl_init";
my $temp_folder='/home/mmariotti/temp/';
if ($opts{"t"}) { $temp_folder = $opts{"t"} ;  }

my $cmnd="";
my $species=$opts{'s'};
my $ext_ID=$opts{'i'};

&debug( "using tempfolder:$temp_folder\n");

Bio::EnsEMBL::Registry->load_all($reg_conf);




my $gene_adaptor=    Bio::EnsEMBL::Registry->get_adaptor(  $species, "core", "Gene");
my $transcript_adaptor=    Bio::EnsEMBL::Registry->get_adaptor(  $species, "core", "Transcript");


my $subj_name;
print ("ext_id: $ext_ID\n");

my $transcript;
my $gene; my $gene_ID;
my @transcripts; # =     @{ $transcript_adaptor->fetch_all_by_external_name( $ext_ID) };

$gene_ID = $opts{'g'};


#
my $description;
my $chrom_name;


print "transcripts found: ". scalar @transcripts . "\n";


$gene = $gene_adaptor->fetch_by_stable_id($gene_ID);

print_DBEntries( $gene->get_all_DBEntries() );



foreach my $transcript ( @{ $gene->get_all_Transcripts() } ) {
    print "TRANSCRIPT ", $transcript->stable_id(), "\n";
    print_DBEntries( $transcript->get_all_DBEntries() );

    # Watch out: pseudogenes have no translation
    if ( defined $transcript->translation() ) {
        my $translation = $transcript->translation();

        print "TRANSLATION ", $translation->stable_id(), "\n";
        print_DBEntries( $translation->get_all_DBEntries() );
    }
}







foreach $transcript (@transcripts){
	
	print "ENSEMBL_TRANSCRIPT_ID:	".$transcript->stable_id."\n";
	$gene = $gene_adaptor->fetch_by_transcript_stable_id                 ($transcript->stable_id);


	$description = $gene->description();
	$chrom_name =  "chr".$gene->slice->seq_region_name;

	print ("GENE_DESCRIPTION: $description\n");
}





sub debug
{	my $msg = $_[0];	
	print STDERR ("$msg") if $opts{v}; }

sub usage
{
    open(HELP, "| more") ;
    print HELP <<EndOfHelp;

This program queries the Ensembl database using an external reference as query (ex: swissprot or gi ), and print stuff for the trascripts annotated for that xeternal reference

USAGE:
	get_ensembl_from_gi [options] -i external_id



options:
- t:	folder for temporary files ##USELESS BY NOW

- v:	verbose mode	

	
EndOfHelp
close(HELP);
	

}


sub print_DBEntries
{
    my $db_entries = shift;

    foreach my $dbe ( @{$db_entries} ) {
        printf "\tXREF %s (%s)\n", $dbe->display_id(), $dbe->dbname();
    }
}
