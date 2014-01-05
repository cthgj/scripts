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




my $prot_adaptor=    Bio::EnsEMBL::Registry->get_adaptor(  $species, "core", "Translation");
#my $transcript_adaptor=    Bio::EnsEMBL::Registry->get_adaptor(  $species, "core", "Transcript");

print $prot_adaptor . "\n";
my $a=$prot_adaptor->fetch_by_dbID($ext_ID);

print $a;
print "\n";




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
