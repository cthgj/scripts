#!/usr/bin/perl -w

use strict;
use Getopt::Std;
require "MY_SUBS.pl";

my (%opts,%synonyms,%interactors,%terms_to_GOs,%want);
getopts('hgdT:N:o:s:S:O:D:C:',\%opts);
$synonyms{LOADED}=0;
my $DATADIR=$opts{D}||"$ENV{HOME}/research/new_moon/data";
my $species;
if ($opts{s}) {
     $species=get_species($opts{s}, "-s");
}
else {
    $species='human';
}
my $degree=$opts{d}||undef;
my $synfile=$opts{S}||"$DATADIR/$species.map";
my $net_file=$opts{N}||"$DATADIR/../$species.gr";
#my $cand=$ARGV[0]||die "USAGE:\n\t$0 <cand_prot>";
my $gaf_file=$opts{O}||"$DATADIR/gene_association.$species";
my $subonto=$opts{o}||'P';
my $go_terms_file=$opts{T}||"$DATADIR/GO.terms_alt_ids";
my $have_already_read_terms_file=0;
my $return_go=$opts{g}||undef;
my $class_file=$opts{C}||"$DATADIR/$species.BP.clas";
usage() unless($ARGV[0]);
usage() if $opts{h};

if (-e $ARGV[0]) {
    my $fh=check_file($ARGV[0], 'r');
    while (<$fh>) {
	chomp;
	$want{&get_name($_,1)}++;
    }
}
else {
    my @a=split(/,/, $ARGV[0]);
    map{$want{&get_name($_,1)}++}@a;
}




open(G,"$gaf_file")||die "Could not open $gaf_file :$!\n";
while(<G>){
    next if /^!/;
    chomp;
    my @tmparray=split(/\t/,$_);
    next unless $tmparray[8] eq $subonto;
    next if $tmparray[3] !~/^$/;
    my $name=&get_name($tmparray[1],1);
    if(defined($want{$name}) and $tmparray[8] eq $subonto){
	$interactors{$name}{TERMS}{&terms_to_GOs($tmparray[4],1)}++;
	$interactors{$name}{GOs}{$tmparray[4]}++;
	# push @{$interactors{$name}{TERMS}}, &terms_to_GOs($tmparray[4],1);
	# push @{$interactors{$name}{GOs}}, $tmparray[4];
    }
   
}
# my @k=keys(%interactors);
# for (my $i=0; $i<scalar(@k); $i++){
#     print &get_name($k[$i],0) . "\t@{$interactors{$k[$i]}{TERMS}}\n";
# }
$"=":::";
foreach my $name (keys(%want)){
    my @terms=keys(%{$interactors{$name}{TERMS}});
    $terms[0]//="None";
    map{
	my $a=&terms_to_GOs($_,0);
	s/$_/$a/}@terms if $return_go;
    #s/\s/_/g for @terms;
    print &get_name($name,0) . "\t@terms\n";
}

######################################################################

sub get_name{
    my $name=shift;
    my $mode=shift;
    if ($synonyms{LOADED}==0){
	open(S,$synfile);
	while(<S>){
	    chomp;
	    my @kk=split(/\t/);
	    $synonyms{NAME}{$kk[1]}=$kk[0]; 
	    $synonyms{ACC}{$kk[0]}=$kk[1]; 
	    $synonyms{NAME}{$kk[0]}=$kk[0]; 
	    $synonyms{ACC}{$kk[1]}=$kk[1]; 
	}
	$synonyms{LOADED}=1;
    }
#    print "\$synonyms{$name} : $synonyms{$name}\n";
    return("aa") unless defined($synonyms{NAME}{$name});
    $mode==0 ? return($synonyms{NAME}{$name}) : return($synonyms{ACC}{$name});
    # if(defined($synonyms{NAME}{$name})){return($synonyms{NAME}{$name})}
    # elsif(defined($synonyms{ACC}{$name})){return($synonyms{ACC}{$name})}
}

############################################################
sub terms_to_GOs{
    my $term=shift;
    my $mode=shift; ## 0 will return GO:xxx, 1 will return term name
#    unless($term eq 'biological_process'){die("crapiola : -$term-\n");}
    $term=~s/_/ /g;
    if($have_already_read_terms_file==0){
	open(T,"$go_terms_file")|| die("Cannot open terms file : $!\n");
	while(my $line=<T>){
	    next if $line=~/^\!/; 
	    chomp($line);
            ## For some reason, latest GO terms file had
	    ## biological_process instead of biological process.
	    ## So, deal with that:
	    $line=~s/_/ /g; 

	    my @aa=split(/\t/, $line);
	    my @a=@aa;
	    my @terms=($aa[0]);
	    if($line=~/obs$/){pop(@a);}
	    if($aa[1] =~/GO:/){
		push @terms,split(/\s+/,$aa[1]);
	    }
	    if($a[$#a] ne $subonto){
		$terms_to_GOs{TERMS}{$a[$#a-1]}="BAD";
		map{$terms_to_GOs{GOs}{$_}="BAD";}@terms;
	    }
	    else{
		$terms_to_GOs{TERMS}{$a[$#a-1]}=$a[0];
		map{$terms_to_GOs{GOs}{$_}=$a[$#a-1];}@terms;
	    }
	}
	close(T);
	$have_already_read_terms_file=1;
    }

#    &debug("term : $term, id:$terms_to_GOs{TERMS}{$term}, id:$terms_to_GOs{TERMS}{$term} " );
#    print STDERR "term : $term, id:$terms_to_GOs{TERMS}{$term} \n";
    $mode==0 ? 
	return($terms_to_GOs{TERMS}{$term}) :
	return($terms_to_GOs{GOs}{$term}) ;
    
}


sub usage{
    my $us="[options] <protein|list of proteins>";
    my $desc="This file will take a single protein name, or a file containing a list and return its annotations.";

    my %opts=(
	      "usage" => $us,
	      "desc" => $desc,
	      "s" => "The species we want, can be either a long (human), or short (hs) form, or taxid. Def: human",
	      "S" => "Synnonyms file mapping Uniprot ACs to IDs",
	      "N" => "Network file, used only if the degree of the protein is requested. Def: \$species.gr",
	      "d" => "Degree, print the degree of each query protein in the network.",
	      "O" => "Gene ontology file, gaf format. Def: \$DATADIR/gene_association.\$species ",
	      "D" => "Data directory. Def: $ENV{HOME}/research/moonlight/data",
	      "o" => "Subontology of interest, one of P,C,F. Def: P",
	      "T" => "GO terms file. Def: \$DATADIR/GO.terms_ids_obs.",
	      "g" => "Return GO terms intead of GO term names, i.e. GO:123 instead of \"signalling\"",
	      "C" => "Class file. Def: \$DATADIR/\$species.BP.clas",
	      "h" => "Print this help and exit"
	     );
    print_help_exit(\%opts);


}
