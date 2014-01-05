#!/usr/bin/perl -w

## Wrapper script to launch calculate_probabilities.pl or calculate_interactome_probabilities.pl

use Getopt::Std;

my %opts;
getopts('hipvea:D:G:g:m:n:o:s:',\%opts); ## inter
usage() if $opts{h};

my @bool_opts=('p','v','h','e');
my @string_opts=('a','D','G','g','m','n','o');

#########################################
# If we are running in interactome mode #
#########################################
my ($args,$bargs,$prog_name);
if($opts{i}){
    $prog_name="calculate_interactome_probabilities_c1.pl";
}
else{
    $prog_name="calculate_probabilities_c1.pl";
}

foreach my $option(@bool_opts){
    next unless $opts{$option};
    $bargs.=$option;
}
foreach my $option(@string_opts){
    next unless $opts{$option};
    $args.=" -$option $opts{$option}";
}
my $cmd;
$bargs?
    ($cmd="$prog_name -$bargs$args") : 
    ($cmd="$prog_name $args");

print STDERR "COMMAND: $cmd\n\n";
system("$cmd");



sub usage{
    $0=~/.+\/(.+)/;
    my $name = $1;
    open(HELP, "| more") ;
    print HELP <<EndOfHelp;

$name is a wrapper script that calls either calculate_probabilities_c1.pl (default) 
or calculate_interactome_probabilities_c1.pl (-i).

USAGE:  
    $name [options] 

OPTIONS:
    -a : Gene (a)nnotation file. GAF format, from GeneOntology.org.
    -d : (D)ebug mode
    -g : Ontology (g)enealogy file. Each line has a go:pair followed 
         by a TAB separated list of its ancestors
    -G : Go alt file (from GeneOntology). Contains ontology information for
         each GO, as well as synonyms, deprecated terms etc. 
	 (def: ./data/GO.terms_alt_ids)
    -i : Calculate (i)nteractome probabilities.
    -o : desired sub(o)ntology, one of F,C,P or combinations thereof, e.g. CF.
    -n : (N)etwork file (for mode -i).
    -m : (M)ap file, allows translation between network name and GAF name.
    -v : Print progress information
EndOfHelp
close(HELP);

exit;
    
}
