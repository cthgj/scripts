#!/usr/bin/env perl 
use Bio::EnsEMBL::Registry;
use autodie;
use Math::Round;
use Getopt::Std;
getopts('vhtua:',\%opts);
if ($opts{h} || !$ARGV[0]) {
    usage()
} 
my $trim=$opts{t}||undef;
my $codon_use=$opts{u}||undef;
###################################
# Get the amino acids of interest #
###################################
my $acids=$opts{a}||die("Need a list of desired amino acids (-a)\n") ;
my @aas=split(",",$acids);
my %wanted_aas;  ## These are the amino acids we are interested in.
my %wanted_cols; ## This will hold the columns we a re interested in.
my %codon_usage; ## This will hold the codon usage stats for each desired codon
$wanted_aas{$_}++ foreach @aas; 

while (<>) {
    next if /^\s*$/; ## Skip empty lines 
    chomp;
    my $to_be_printed; ## This will be the string printed for each line if run with -t
    my @columns=(split(/\t/));
    ##############################################
    # If this is the first line, get the headers #
    ##############################################
    if ($.==1) {
	$to_be_printed.="$columns[0]" if $trim; ## Print the 1st data column
	$wanted_cols{0}=1;     ## Keep the 1st data column
    	for (1..$#columns){
	    #############################################################
	    # Data columns are of the format "AA_CODON" or "AA_CODON_%" #
	    #############################################################
	    if ($columns[$_]=~/^..._/ || $columns[$_]=~/^Stop_/) {
		my @a=split(/_/,$columns[$_]);
		my $aa=shift(@a); ## This is the current column's amino acid
		if (defined($wanted_aas{$aa})){
		    $wanted_cols{$_}=1;
		    $to_be_printed.= "\t$columns[$_]" if $trim; 
		}
	    }
	    ####################################
            # Get the gene info columns	       #
            ####################################
	    else {
		$wanted_cols{$_}=1;
		$to_be_printed.=  "\t$columns[$_]" if $trim;
	    }
	}
    print "$to_be_printed\n" if $trim;
    }
    ##########################################
    # If this is not the first line, process #
    ##########################################
    else {
	my @want=sort(keys(%wanted_cols));
	
	######################################################
        # If we are just trimming the file, print the	     #
	# desired columns.				     #
        ######################################################
	if ($trim){
	    $to_be_printed.=$columns[shift(@want)];
	    foreach (@want){
		$to_be_printed.= "\t$columns[$_]" if defined($wanted_cols{$_});
	    }
	print "$to_be_printed\n" if $trim;
	}
	#############################################
        # If we are calculating codon usages	    #
        #############################################
	elsif ($codon_use) {
	    
	}
    }

}

sub usage{
    my $script=`basename $0`;
    chomp($script);
    print STDERR <<EndOfHelp;

USAGE: 

    $script [OPIONS] -a <AminoAcid> <SPECIES NAME>

DESCRIPTION:
 
    This script will parse the codon usage file produced by ensembl_get_codon_count.pl.
    It will return the usage statistics for the requested amino acids. Multiple amino
    acids can be specified as a comma separated list.

     Use '>' to redirect the output to a file (it will be lost if you do not do this).

EXAMPLES:

     $script -a Ser,Thr human.csv > human.Ser.stats
     $script -v "Saccharomyces cerevisiae"  > yeast.csv

OPTIONS:
    
    -a : Three letter amino acid code. Only data for this aa will be returned.
         Multiple amino acids can be given as a comma separated list, eg 
         Thr,Ser,Cys.
    -h : Print this help and exit.
    -t : Trim the file, this will print only the data columns (name, ID etc) and 
         the columns of the amino acids you asked for.
    -u : Print the average codon usage for each specified amino acid.
    -v : Verbose, print progress messages.
    
EndOfHelp
exit;


}
