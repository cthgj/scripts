#!/usr/bin/perl -w
use strict;
use lib '/users/rg/mmariotti/libraries/ensembl_52/modules/' ;   # --> path to ensembl API software. must be same version of genome and annotation releases
use Getopt::Std;

use Bio::EnsEMBL::DBSQL::DBAdaptor;
my %opts;
getopts('EaYFhCvpS:s:e:c:t:r:o:d:f:I:y:',\%opts);
if ($opts{'h'}){&usage; die();}
if (! $opts{"v"}) {$opts{"v"}=0}
my $db_source='core';
if ($opts{"d"}) { $db_source = $opts{"d"} ;  }

use Bio::EnsEMBL::Registry;

my $reg_conf="/users/rg/mmariotti/Ensembl_genomes/.ensembl_init"; # --> ensembl registry. if you don't want to use it, comment this and the next line and uncomment the other one (my $registry...) 
Bio::EnsEMBL::Registry->load_all($reg_conf);

#my $registry = Bio::EnsEMBL::Registry->load_registry_from_db(    -host => 'ensembldb.ensembl.org',    -user => 'anonymous',    -verbose => "".$opts{"v"}. "",    -port => 5306);  # --> this will load the registry online, so the most recent annotations will be fetched. careful cause the latest release you can download depends on the software version you have.

my $temp_folder='/home/mmariotti/temp/';
if ($opts{"t"}) { $temp_folder = $opts{"t"} ;  }
my $seq_to_find= $opts{"S"};
my $species =  $opts{"o"} ;

my $cmnd="";
my $chrom_number; 
my $pos_start   ; 
my $pos_end  ;
my $to_find_is_minus_strand;
my $t; my @a;
my $coord_system_to_use="chromosome";

if ($opts{"f"}){

	open(GFF_FILE, "< " . $opts{"f"} );
	my $gff_line= <GFF_FILE>; 	 
	 while ($gff_line){
		 chomp $gff_line; 
		#print "$gff_line \n" ;

		 $gff_line =~  /^([^\t]+)\t([^\t]+)\t([^\t]+)\t(\d+)\t(\d+)\t([^\t]+)\t(\+|\-)\t.*$/ || die ("Error, the gff file is badly formatted. $!");
	
#		 $gff_line =~  /^(\w+)\t(\w+).+$/ || die ("Error, the gff file is badly formatted. $!");

		 if ( !$pos_start || $4<$pos_start ){
			 $pos_start   = $4 ; 
		 }
		 if ( !$pos_end || $5>$pos_end ){
			 $pos_end  = $5 ;
		 }

		 $to_find_is_minus_strand= ($7 eq "-"  );
	
		 $chrom_number = $1 ;
		 $gff_line= <GFF_FILE>;

	 } 
	close(GFF_FILE)
}


else{
 $chrom_number = $opts{"c"} ; 
 $pos_start   = $opts{"s"} ; 
 $pos_end  = $opts{"e"} ;
 $to_find_is_minus_strand= $opts{"r"}; if (! $to_find_is_minus_strand) { $to_find_is_minus_strand= 0 ; }
}


if (! $seq_to_find ){&debug("Protein sequence was not defined. use option -S "); $seq_to_find=''} else{ &debug("seq_to_find=$seq_to_find\n");}
&debug( "using tempfolder:$temp_folder\n");
&debug( "to_find_is_minus_strand: " .$to_find_is_minus_strand ."\n");
&debug( "chrom_number: $chrom_number \n");
&debug( "pos start: $pos_start \n");
&debug( "pos_ens: $pos_end \n");

$cmnd='echo ">gene_to_find_seq'."\n$seq_to_find".'" > '."$temp_folder/seq_to_find";
`$cmnd`;





 #BUGGG HERE!

my $slice_adaptor=    Bio::EnsEMBL::Registry->get_adaptor(  $species, $db_source, "Slice"); #$query->{species}





if ($opts{"C"}){		#check mode 
	my $slice = $slice_adaptor->fetch_by_region('toplevel', $chrom_number, 1 , 30 ); 
	die("TEST OK")}

#$slice_adaptor->fetch_all('toplevel');

#my $slice = $slice_adaptor->fetch_by_region('chromosome', $chrom_number, $pos_start , $pos_end );
my $slice = $slice_adaptor->fetch_by_region( 'toplevel', $chrom_number, $pos_start , $pos_end );

my $gene;
my $current_score = 0;
my $max_score = 0;
my $best_transcript;
my $best_gene='pituffo';	
my $gene_desc='';
my $transcript;
my $strand;
my $seq_obj;
my $translation;
my $desc;
my $gene_strand;
my $skip_others=0;
my $temp;

foreach $gene (@{$slice->get_all_Genes}) {

	&debug("GENE : $gene ID: ".$gene->stable_id(). " "); 		#
	&debug(" biotype: " .$gene->biotype . " n_transcripts:" . scalar @{$gene->get_all_Transcripts} . "\n"  );
	next if $skip_others;

	next if $gene->biotype ne "protein_coding";
	$gene_strand= $gene->strand();
	&debug("---- strand:".	$gene_strand . "\n");
	next if ($to_find_is_minus_strand &&  $gene_strand>0 ) || (  ! $to_find_is_minus_strand &&  $gene_strand<0 );



	#
	$gene_desc=$gene->description();
	#&debug ("GENE:$gene_desc\n");
#	my $geneseq = $gene->seq();

	foreach $transcript (@{$gene->get_all_Transcripts}){

		next if $skip_others;

		$desc =  $transcript->description() || &debug( "ERROR: NO TRASCRIPT DESCRIPTION\n");
		&debug("TRANSCRIPT : $transcript ID: ".$transcript->stable_id(). " \n"); 		#
#		&debug("TRANSCRIPT desc: $desc\n");
		
		if ($opts{'a'}) 	{		

			if ($opts{'F'}){ print  &get_fasta($transcript,$gene_desc)}
			else{ print &get_gff($transcript, $desc) }
			$best_gene='something found'
		}




		$seq_obj =  $transcript->seq();		#Bio::seq	object

		$strand =  $transcript->strand();
		$translation =  $transcript->translate()->seq();		#peptide seq!!!

		&debug("strand: $strand\n");
		&debug("translation: $translation\n");


#		my $dna_seq= $seq_obj->seq();			#actual DNA seq
#		print("dna seq:$dna_seq\n");

		if (! $opts{'S'} && !$opts{'a'} )  {			
			$best_transcript = $transcript;
			$best_gene = $gene;  
			$skip_others=1;
			next;
				}


		$cmnd='echo ">current_transcript_seq'."\n$translation".'" > ' . " $temp_folder/current_transcript_seq  ; cat $temp_folder/seq_to_find >> $temp_folder/current_transcript_seq " ;
		`$cmnd`;
		&debug ("launching tcoffee\n");
		$cmnd="cd $temp_folder; t_coffee current_transcript_seq  > t_coffee_out 2>/dev/null";
		`$cmnd`;

		&debug ( "reading $temp_folder/t_coffee_out\n");
		$cmnd="grep -o \"*\" $temp_folder/t_coffee_out  | wc";
		$t= `$cmnd` ;
		@a=split(/ +/, $t);
		$current_score=$a[1];


		&debug ( "score: $current_score\n");
		if ( $current_score> $max_score) {
			$max_score = $current_score;
			$best_transcript = $transcript;
			$best_gene = $gene;
			&debug("setting BEST gene and transcript\n");
}



		
}


}

if ($best_gene eq 'pituffo'){
	die("no protein coding-genes with the specified orientation in that slice")
}
else {

	if (!$opts{'a'}){
	
		if ($opts{'F'}	) {
			print  &get_fasta($best_transcript,$best_gene->description()) 
		}
		else{
			&debug( "BEST ONE: $best_gene $best_transcript\n");
			$temp = &get_gff($best_transcript, $best_gene->description());

		print ("$temp");
		}

	}
}


&debug("end\n");

sub get_fasta{
	my $transcript=$_[0];
	my $title=$_[1];
	my $out='';

	if (! $title){ $title = "Ensembl_transcript:" . $transcript->stable_id() }
	else {$title .= " (Ensembl_transcript:" . $transcript->stable_id() . ")"}
	$out.=">". $title . "\n" . $transcript->translate()->seq() ."\n";
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

	
	#print STDERR $coord_system_to_use . "\n";

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


	my @exons = @{$transcript_obj->get_all_translateable_Exons}	;	#CHECK DIFFERENCES
	my @introns = @{$transcript_obj->get_all_Introns};
	my $ex_index = 0;
	my $ex_end;
	my $ex_start;

		
	&debug("using coordinates system $coord_system_to_use \n");

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


sub debug
{	my $msg = $_[0];	
	print STDERR ("$msg") if $opts{"v"}; }

sub usage
{   


    open(HELP, "| more") ;
    print HELP <<EndOfHelp;
This program fetch the most similar annotated transcript in Ensembl to a gene prediction given as a gff file. 
usage:  get_annotation  -o species -S prot_sequence  -f gff_file  
or: 	get_annotation  -o species -S prot_sequence  -c chromosome_name -s start -e end  [-r] #->is gene is on minus strand 
    
    in case protein sequence is provided, the transcript with translation featuring the highest score when aligned to given protein sequence is fetched.
    in case it is not provided, the first one will be reported.
    
    options:

	-d database source. "core" by default, but you could specify "vega" as well for some genomes.
    	-C run in check mode, to see if any genes could be fecthed. it tests only the organism name. if the test goes ok, it prints TEST OK and dies.
	-y coordinate system name
	-Y print available coordinate systems and die
	-v verbose
	-a get all transcripts in specified range
	-E essential output: only cds lines
	
EndOfHelp
close(HELP);
    exit(1);

}
