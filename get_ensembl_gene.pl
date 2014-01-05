#!/usr/bin/perl -w
use strict;
use lib '/users/rg/mmariotti/libraries/ensembl_52/modules/' ;
#use lib '/home/ug/talioto/myperlmods/ensembl/modules';
use Getopt::Std;

use Bio::EnsEMBL::DBSQL::DBAdaptor;

my %opts;
getopts('gYsTho:i:t:h',\%opts);
if ($opts{'h'}){&usage; die();}
#my $reg_conf="/users/rg/mmariotti/.ensembl_init";


use Bio::EnsEMBL::Registry;
#Bio::EnsEMBL::Registry->load_all($reg_conf);
#Bio::EnsEMBL::Registry->load_registry_from_db(    -host => 'ensembldb.ensembl.org',    -user => 'anonymous',    -port => 5306);

my $reg_conf="/users/rg/mmariotti/Ensembl_genomes/.ensembl_init"; 
Bio::EnsEMBL::Registry->load_all($reg_conf);


my $species;
$species=$opts{'o'};

my $gene;
my $gene_adaptor=    Bio::EnsEMBL::Registry->get_adaptor(  $species, "core", "Gene"); #$query->{species}

if ($opts{'i'}){
	$gene = $gene_adaptor->fetch_by_stable_id($opts{'i'})   ;
}
elsif ($opts{'t'}){
	$gene = $gene_adaptor->fetch_by_transcript_stable_id($opts{'t'}) || die('ERROR: cannot fetch gene from transcript ID'. $opts{'t'})  ;
}

#print "DESCRIPTION " . $gene->get_description() ."\n";

my $transcript;

my $chrom_number=$gene->slice->seq_region_name();
my $coord_system_to_use="chromosome";
my $to_find_is_minus_strand;

if ($gene->strand() == -1 ){ $to_find_is_minus_strand=1}

my $transcripts = $gene->get_all_Transcripts();
while ( $transcript = shift @{$transcripts} ) {

	if (!$opts{'T'} || $transcript->stable_id eq $opts{'t'} ){


				if (! $opts{'g'}){
	        print "# " .$gene->stable_id ."\t" .feature2string($transcript) . "\n";
				}
				if (! $opts{'s'}){				
					print &get_gff($transcript, $gene->description);
				}
	}
}



sub feature2string
{
    my $feature = shift;

    my $stable_id  = $feature->stable_id();
    my $seq_region = $feature->slice->seq_region_name();
    my $start      = $feature->start();
    my $end        = $feature->end();
    my $strand     = $feature->strand();

    my $out= sprintf( "%s : %s:%d-%d (%+d)",
        $stable_id, $seq_region, $start, $end, $strand );

    if ( $feature->biotype() eq 'protein_coding' ) {
	    my $translation =  $feature->translation()->seq();
	    $out .= " " .$translation ;
    }

    return $out
    
}


sub get_gff{

	my $transcript_obj= $_[0];
	my $description=$_[1];

	my $outgff='';
	my $notes;


	my $stable_id=$transcript_obj->stable_id();
	$notes= "$stable_id (EnsEMBL ID)" ; #$exon->seq()->seq();
	my $gene_strand;
	if ($transcript_obj->strand() > 0 ){ $gene_strand='+'}
	else {$gene_strand= '-'}
	if ($opts{'p'}) {print ("TRANSL: ". $transcript_obj->translation()->seq(). "\n"); }


	my $temp; my $new ;

	  my $csa = Bio::EnsEMBL::Registry->get_adaptor(  $species, "core", "coordsystem");

	  #
	  # Get all coord systems in the database:
	  #
	my @coord_systems= @{ $csa->fetch_all() };

	if ($opts{'Y'}) {
		foreach $temp ( @coord_systems) {
			print $temp->name() . "\n";
		};
		die()
	}
	$temp=0;



	foreach $new (@coord_systems){
		if ("chromosome" eq $new->name() ) {
			$temp=1;
		}
	}

	if (!$temp){
		
		my $candidate;

		if ($chrom_number =~ /(\w+)_(\d+)/){
			$candidate=$1;
		}
		elsif ($chrom_number =~ /^Ultra/){
			$candidate="ultracontig";
		}
		elsif ($chrom_number =~ /^Contig/){
			$candidate="contig";
		}
		elsif ($chrom_number =~ /^supercont/){
			$candidate="supercontig";
		}

	
		if ( $candidate) {
			foreach $new (@coord_systems){
				if ($candidate eq $new->name() ) {
					$coord_system_to_use=$candidate;
					$temp=1;
				}
			}
		}
	}
	
	if (!$temp){	$coord_system_to_use = $coord_systems[0]->name()}

	

	$temp=0;
	if( $opts{'y'}){
		$coord_system_to_use=$opts{'y'};
		foreach $new (@coord_systems){
			if ($new->name() eq $coord_system_to_use){
				$temp=1;
			}
		}
		if (!$temp){
			print "ERROR the specified coordinate system is not available\n";
			die()
		}
	}




	my $subj_name=$chrom_number;
	if ( $coord_system_to_use eq "chromosome" ) {
		#$subj_name='chr' . $subj_name;
	}


	$new= $transcript_obj->transform($coord_system_to_use);
	my $index_system=0;
	while (! $new && $index_system <scalar(@coord_systems)){
		
		$coord_system_to_use=$coord_systems[$index_system]->name();
		$new= $transcript_obj->transform($coord_system_to_use);
		$index_system+=1;
	}

	my $cds_start= $new->coding_region_start();
	my $cds_end = $new->coding_region_end();


	if (! $description){$description = $notes}
	if (! $opts{'E'} ){
		$outgff.= "$subj_name\tEnsEMBL_gene:$stable_id\tgene\t" . $new->start() ."\t" . $new->end() ."\t.\t$gene_strand\t.\t$description\n"  ;
	}
	my $exon;
	my $intron;


	my @exons = @{$transcript_obj->get_all_translateable_Exons}	;			#CHECK DIFFERENCES
	my @introns = @{$transcript_obj->get_all_Introns};
	my $ex_index = 0;
	my $ex_end;
	my $ex_start;

		

	$exon = $exons[$ex_index]->transform($coord_system_to_use);

	$outgff.= "$subj_name\tEnsEMBL_gene:$stable_id\tcds\t" . $exon->start() ."\t" . $exon->end() ."\t.\t$gene_strand\t.\t$notes\n"  ;
	$ex_index++;

	while ($ex_index < scalar @exons)	{
		$intron = $introns[$ex_index-1]->transform($coord_system_to_use);
#		$notes = "$stable_id (EnsEMBL ID)";

		if (! $opts{'E'} ){
			$outgff.= "$subj_name\tEnsEMBL_gene:$stable_id\tintron\t" . $intron->start() ."\t" . $intron->end() ."\t.\t$gene_strand\t.\t$notes\n"  ;
		}
		$exon = $exons[$ex_index]->transform($coord_system_to_use);
#		$notes= ''; #$exon->seq()->seq();

		$ex_index++;
		$ex_start=$exon->start();
		$ex_end = $exon->end();


		if ($ex_index == scalar @exons){
			#it was last exon. contains stop codon. Deleting it-

			if ($to_find_is_minus_strand){
				$ex_start+=3;
			}
			else{   $ex_end-=3 }
}

		$outgff.= "$subj_name\tEnsEMBL_gene:$stable_id\tcds\t" .$ex_start ."\t" . $ex_end ."\t.\t$gene_strand\t.\t$notes\n"  ;

}



	return $outgff

}


sub usage
{   


    open(HELP, "| more") ;
    print HELP <<EndOfHelp;
This program utilizes Ensembl API and fetches a gene annotation, printing some information for each transcript.
usage:  get_ensembl_gene.pl  -o species -i gene_ensembl_id
or:	get_ensembl_gene.pl  -o species -t one_transcript_ensembl_id

with no options, full output is provided. This includes summary, a line for each transcript like this:

# ENSDORG00000013088	ENSDORT00000013089 : GeneScaffold_5275:14580-24661 (+1) MASVRSAVGPLLAAASARPVPFGRPLPPPAPRSAWAAAMEPVPRWLAGLRFDNRALRELPVEAPPPGPEGAPSAPRPVPGACFARARPVPLQRPRVVALSGPALALLGLDAPPAAEAEAALFFSGNALLPGAEPAAHCYCGHQFGQFAGQLGDGAAMYLGEVCTAAGGRWELQLKGAGPTPFSRQADGRKVLRSSIREFLCSEAMFHLGIPTTRAGACVTSESTVVRDVFYDGNPRYEKCAVVLRIAPTFIRFGSFEIFKSTDEHTGRAGPSVGRNDIRIQMLNYVISSFFPEIHTAHACESSPVQRNAAFLREVTRRTAQMVAEWQCVGFCHGVLNTDNMSIVGLTIDYGPFGFLDRYDPDHVCNASDNAGRYAYCKQPEVCKWNLHKLAEALEPELPLALGEAIIAEEFDAEFQKHYLHKMRRKLGLVRMEQDEDRALVARLLETMRLTGADFTNTFSVLSTFPVETEAGLSDFLAVLTSQCASLEELKLAYRPQMDPRQLSMMLMLAQSNPQLFALIGTQANVTKELERMEQQSRLEQLSPAELQNRNEEHWTAWLQEYRARLDKDKECAGDTAAWQAERVQVMHMNNPKYVLRNYIAQKAIEAAENGDFSEVRRVLKLLETPYHRDGEAAAGLEATSPEGASGGQCSYSSRPPLWAAELCVTSS

And on the next line, gff output:

GeneScaffold_5275       EnsEMBL_gene:ENSDORT00000013089 gene    14580   24661   .       +       .       RIKEN cDNA 1300018J18 gene Gene [Source:MGI Symbol;Acc:MGI:1919007]
GeneScaffold_5275       EnsEMBL_gene:ENSDORT00000013089 cds     14580   14911   .       +       .       ENSDORT00000013089 (EnsEMBL ID)
GeneScaffold_5275       EnsEMBL_gene:ENSDORT00000013089 intron  14912   14912   .       +       .       ENSDORT00000013089 (EnsEMBL ID)
GeneScaffold_5275       EnsEMBL_gene:ENSDORT00000013089 cds     14913   14922   .       +       .       ENSDORT00000013089 (EnsEMBL ID)



options:
-g	gff output only. 
-s  summary output only.
-T	print output only for the transcript specified with -t (works as a get_ensembl_transcript)

-Y  print available coordinate systems and die
-y [ARG]	use specified coordinate system   (normally, a builtin parser try to find the correct coordinate system)

EndOfHelp
close(HELP);
    exit(1);

}

