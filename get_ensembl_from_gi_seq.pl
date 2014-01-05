#!/usr/bin/perl -w
use strict;
use lib '/home/ug/mmariotti/libraries/ensembl_api/modules/' ;
#use lib '/home/ug/talioto/myperlmods/ensembl/modules';
use Getopt::Std;
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::DBSQL::DBAdaptor;

my %opts;
getopts('hvt:i:s:g:',\%opts);
&usage() if $opts{'h'};

my $reg_conf="/home/ug/mmariotti/.ensembl_init";
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
print STDERR ("ext_id: $ext_ID\n");

my $transcript;
my $gene; my $gene_ID;
my @transcripts; # =     @{ $transcript_adaptor->fetch_all_by_external_name( $ext_ID) };




#
my $description;
my $chrom_name;


print STDERR "transcripts found: ". scalar @transcripts . "\n";

#$gene = $gene_adaptor->fetch_by_stable_id($gene_ID);

$transcript = $transcript_adaptor->fetch_by_stable_id($ext_ID);






#foreach my $transcript ( @{ $gene->get_all_Transcripts() } ) {
    if(defined($transcript)){#->stable_id())){
#    die("shit\n");
    }
    else{
    print STDERR "trying somethingn else\n";
    @transcripts =@{  $transcript_adaptor->fetch_all_by_external_name($ext_ID)};
    foreach my $transcript ( @transcripts  ) {
    
      print "TRANSCRIPT ", $transcript->stable_id(), "\n";
    }
    
    

#    $transcript= $transcripts[0];
 
    
    }
    
    
    
    print_DBEntries( $transcript->get_all_DBEntries() );

    # Watch out: pseudogenes have no translation
    if ( defined $transcript->translation() ) {

        my $translation = $transcript->translation();
	my $seq= $transcript->seq()->seq();
	print ">".$transcript->stable_id() ."\n$seq\n";


#        print "TRANSLATION ", $translation->stable_id(), "\n";
#        print_DBEntries( $translation->get_all_DBEntries() );
    }
#}







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
