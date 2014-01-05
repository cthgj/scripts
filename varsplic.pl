#!/usr/bin/perl

# POD documentation - main docs before the code

=head1 NAME - varsplic.pl

=head1 DESCRIPTION

Generates sequence for alternative protein isoforms from  
UniProt Knowledgebase (UniProtKB) entries

=head1 USAGE

perl varsplic.pl -input filename [-error filename] [-stats filename] 
[-pseudo filename] [-fasta filename] [-statsfile filename] [-threshold int] 
[-which [full|allforms]] [-linelength int]  
[-uniqspids] [-showdesc] [-locations] [-dbcode] [note] [-describe filename] 
[-filter] [-count] [-noftids] [-crosscheck] [-check_vsps] [-varseq]
[-variant] [-conflict] [-varsplic DEPRECATED] [-uniqids int DEPRECATED]
[initiator DEPRECATED]

the input option is used to specify a file of entries in UniProt format to be
expanded. If this option is not specified, the program will read from STDIN.

if the error option is specified, error messages are printed to a named file.

if the pseudo option is specified, all new isoforms are printed to file in 
pseudo UniProt format.  If no filename is specified, the program will print
to STDOUT (note that to produce both FASTA and pseudo-UniProt format output,
you must specify a destination file for at least one type of output.  The other
will default to STDOUT if not specified)

if the fasta option is specified, all new isoforms are printed to file in FASTA 
format. If no filename is specified, the program will print to STDOUT (note 
that to produce both FASTA and pseudo-UniProt format output, you must specify a
destination file for at least one type of output.  The other will default to
STDOUT if not specified)

if the stats option is specified, statistics are printed to file

if the statsfile option is specified, program will read in numbers between 0
and 100 (new-line delineated) and use these to derive categories for
statistical analysis, i.e. if the numbers 20 and 40 are specified, program will
count the number of variants showing between  0 and 20% change, 40 and 60%
change and 60 and 100% change.  If option is not specified, default categories
are used.

if the threshold option is specified, new isoforms are only printed to file if 
their sequences show greater than the specified percentage change from the 
primary variant.

if the which option is specified, the user can choose to print new records for
existing isoforms as well as for new isoforms.

if the user specifies '-which full', new records (in the specified format) are
produced for all existing records in the input file as well as for all 
alternative isoforms.

if the user specifies '-which allforms', new records are produced for all 
isoforms of all alternatively-spliced proteins (including the isoforms listed 
in the input file).  New records are not produced from records currently in
the database that have no known splice variants.

the default (no 'which' option specified) is to produce new records only for
isoforms not present in the input file.

if the linelength option is specified, program will take the specified value as
the maximum possible linelength for the first line of a FASTA entry.

the uniqspids option is used to modify the ID lines of entries produced in
pseudo-UniProt format.  It does not affect entries produced in FASTA format.
If specified, any suffix applied to the parental AC to create the unique
isoform_ID of a variant will also be applied to the ID, e.g. if the parental
ID is P12345 the ID should be in the format ABCD_ECOLI-6.

if the showdesc option is specified, the program will print out the full
UniProtKB/TrEMBL description in the FASTA header line (followed by a
description of the variant).  Default behaviour is only to print out the
variant description instead of the full description (for sequences derived
from entries with > 1 variant).

if the note option is specified, the program will add isoform-specific free
text notes to the header lines of FASTA files (or to the DE lines of
pseudo-UniProt format files), after the normal description.   This feature only
works when performing VARSPLIC expansion only (i.e. it will not work if another
type of feature, or any combination of features, are expanded). 

if the locations option is specified, the program will add the co-ordinates of
the first and last changed amino acid in the original sequence to the header
lines of FASTA files (or to the DE lines of pseudo-UniProt format files), after
the normal description.  A semi-colon (and a space) separates the rest of the 
description from the location, which has the following form: 
Change_Location:<first_co_ordinate>-<last_co_ordinate>.

if the describe option is specified, program will print an output file
describing the co-ordinates of the new variants in terms of the parental
variant.  The format used is: 
parental_accession_number : isoform_identifier : 
start position of first region of alternative sequence , 
  stop position of first region of alternative sequence, 
  length of sequence inserted in variant, 
  length of sequence removed in variant; 
start position of second region of alternative sequence etc.

if the filter option is specified, program will filter out any lines beginning
'**' or '++' from input file (used internally by the UniProt production team) 
from input file before attempting processing.  If such lines are present in the
input file, this option should be specified to guarantee correct program
behaviour.

if the count option is specified, program will report to STDERR each time it 
has processed 1000 entries.

if the noftids option is specified, the program will expand VARSPLIC sequences
without requiring the existence of FTIds.

if the crosscheck option is specified, the program will report on cases where
two cross-referencing entries both describe identical sequences.

if the check_vsps option is specified, the program will report on cases
where the VSP lists in the CC and FT sections are inconsistent.

EXTRA OPTIONS

The default behaviour of the program is to produce protein sequences for all
isoforms described using the features identified by the 'VAR_SEQ' key
(features identified by  the 'VARSPLIC' and 'INIT_MET' keys prior to UniProt
release 8.0).  These features describe (in conjunction with additional
information provided in the ALTERNATIVE PRODUCTS comment lines) protein 
isoforms generated through alternative splicing and/or initiation patterns to 
those that produce the sequence displayed in the entry.

It is additionally possible to generate protein sequences for natural variants
and experimental artifacts (described in the feature table of an entry using 
the 'VARIANT' and 'CONFLICT' keys); and to generate protein sequences using
specified combinations of these feature types.  Output is produced in a
slightly altered format, reflecting the fact that a given form may, for example
be created by including a particular conflict in a particular splice variant.

if varseq option is specified, program performs varsplic expansion

if variant option is specified, program performs variant expansion

if conflict option is specified, program performs conflict expansion

with the release of Uniprot 8.0, the init_met option has been discontinued (as 
alternative initiation is now described in the same way as alternative 
splicing).  The 'varsplic' option has been retained for backwards 
compatibility, but it's behaviour is now identical to that of the new 'varseq' 
option, i.e. it will generate sequences for alternatively initiated proteins, 
as well as for alternatively spliced proteins.

NOTE that it is not necessary to specify the 'varseq' or 'varsplic' options to
perform expansion of VAR_SEQ features.  Expansion of these features is 
performed by default.  Use of the 'varseq'/'varsplic' option is necessary only
to perform expansion of 'VAR_SEQ' features in conjunction with expansion of 
'VARIANT' and 'CONFLICT' features.

if more than one of the varseq/variant/conflict options are specified, the
program produces records for all combinations, so if there are 2 additional 
splicing isoforms (in addition to the form shown in the 'parent' record), 3
variant forms and 4 conflict forms, the program will produce 3 x 4 x 5 = 60
records for this entry (unless no 'which' option is specified, in which case 59
records will be produced, i.e. the sequence displayed in the original entry 
will be omitted). See note on 'combinations' below for a further example.

=head1 DEPRECATED

the varsplic option is now deprecated, its functionality is duplicated by the
varseq option since file format changes made in UniProtKB release 8.0.

the uniqids option can still be specified but no longer has any effect.  This 
is because under the new UniProtKB FASTA format (since UniProtKB release 9.0), 
unique IDs are produced automatically using the standard Swissknife toFasta() 
method.

the initiator option can still be specified, but will only have an effect if 
the program is applied to data preceding UniProtKB release 9.4. fs specified, a 
methionine will be added to the first base of all sequences generated from a 
UniProtKB entry (provided this particular variant's sequence starts at position 
1 of the displayed sequence) if the entry contains annotation indicating the 
presence of an (unshown) initiator methionine.  No such annotation exists in
entries from UniProtKB 9.4 onwards (iniatator methionines are always included
in the sequence of an entry).

the dbcode option could previously be used to specify the insertion of a 
database-specific code into the FASTA header line.  The presence of such a code
is standard in the UniProt FASTA header since UniProtKB release 14.0 and this 
option is now deprecated and without function.

=head1 REQUIRES

The Swissknife package, see 
ftp://ftp.ebi.ac.uk/pub/software/swissprot/Swissknife/

=head1  ALGORITHM SKETCH
 
-read in data 1 record at a time

  -for each of varseq, variant and conflict

    -parse the relevant feature tables

    -combine information from the CC blocks and feature tables as appropriate
     to complete the isoform - feature mapping

    -are there any variants?

      - if no

        -print parent record if specified

      - if yes

        -for each isoform (representing 1 possible combination of features of
         each specified key type)

         -for each feature

           -is this feature found in this isoform?

             -if yes

               -is it possible to add this feature to those already
                considered in constructing this (hypothetical) isoform?

                -if yes

                  -process sequence changes defined in this feature

         - print alternative isoform to file as specified            

  

=head1 INFORMATION RELATING TO SPECIFICATION OF VARSEQ, VARIANT, CONFLICT AND
INIT_MET OPTIONS

Each VAR_SEQ, VARIANT, CONFLICT line describes one location at
which the sequence of the protein has been reported to differ from that
displayed in the entry.  Each alternative isoform will consist of a
particular combination of these 'variants'.  For example, a given isoform may
contain 2/3 alternative sequences reported in VAR_SEQ features, but the
conflicts reported in the CONFLICT features may only have been found in
isoforms whose splice pattern matches that displayed in the entry).

This information (if which combination of sequence features have been reported
in actual isoforms) is not necessarily present or parsable in UniProtKB 
records. Therefore, the following approach has been taken:

VAR_SEQ

The 'VAR_SEQ' feature is consistently annotated.  The information of which
(alternatively spliced and initiated) isoforms exist, and how their sequence 
can be derived from the sequence displayed in the entry, is clearly given, and
is reliably parsed by this program.  Isoform names are parsed from the feature
description.

VARIANT

VARIANT features are used to annotate naturally occurring polymorphisms.

A polymorphism annotated in a 'VARIANT' feature may belong to a particular
clone or strain, in which case similar processing can be performed to that
performed for 'VAR_SEQ' features, i.e. a single new record is produced for
each clone or strain, containing all the polymorphisms that characterise it.

However, polymorphisms  of one entry may also have been found in 
patients suffering from a particular disease.  This does not mean that these
features occur concurrently.  When a disease is referred to in an VARIANT
feature in a UniProtKB entry, this can be identified because the name of that
disease will also occur in a corresponding CC DISEASE line. If the polymorphism
is identified in this way as being present in a disease (as opposed to a clone
or a strain), a separate sequence is generated for each polymorphism, such that
each new form has only 1 difference from the parent sequence.

The names of variant forms are parsed from the feature description.  Where no
name is given, where extant the feature ID is used.  Where there is no feature
ID, the form is named 'VARIANT n', where this form is the nth variant form
mentioned in the record.

For forms derived for polymorphisms features that occur in a disease, the
nomenclature is: disease_name-form_identifier.  The form identifier
is, as for unnamed variants (i) if possible, a feature ID (ii) if not, the
number of this variant in order of reference in the parent record.

The annotation of 'VARIANT' features is quite irregular, this limits the
program's performance in parsing them correctly.

CONFLICT

Conflicts are derived from references.  It is assumed that all the conflicts
mentioned in a particular reference are found concurrently, unless there is
specific annotation that states otherwise, i.e. if two 'CONFLICT' features are
described as being 'IN REF 1', a single new record will be generated containing
both conflicts.

Where two conflicting sequences have been obtained from 1 reference, these are
often distinguished by comments located after a semi-colon in the feature
description.  Comments after a semi-colon are thus considered part of the name
of a conflict form (e.g. a conflict may be identified by  'REF. 1' or 'REF. 1;
IN AAB40657'.  By contrast, for 'VAR_SEQ' and 'VARIANT' features, comments
after semi-colons are usually not informative for the purpose of determining
which features are found in which forms, and are therefore ignored.

INIT_MET

Use of the INIT_MET feature key in a UniProtKB record has traditionally
indicated one of two possible things, according to the position of the 
feature.  In position 0, it indicates the presence of an initiator methionine 
cleaved off the primary translation during the production of the functional 
protein product, and not included in the sequence displayed in the entry.  In 
any other position, it indicates an initiator methioine that is an alternative 
to the usual initiator methioine at positions 0 or 1.

With the release of UniProtKB release 8.0, the existence internal initiators is 
now described using the VAR_SEQ key.  Use of INIT_MET to indicate a 
subsequently cleaved initiator is also being discontinued (in future, the 
methionine will be included at position 1 of the sequence, and a CHAIN feature 
beginning at position 2 will be used to represent the mature peptide.  However, 
some 'INIT_MET' features located at position 0 have not yet been removed from 
the database.  The "missing" methionines can be replaced by using the 
'initiator' option (in combination with any of the 'varseq', 'variant' and 
'conflict' options).

With the release of UniProtKB release 9.4, all legacy data has been converted
to the new format and no INIT_MET features located as position 1 remain in the
database

COMBINATIONS

No attempt is made to determine which conflicts are found in which splice 
forms, etc.  Instead, where > 1 option is specified (from 'varseq', 'variant',
'conflict' and 'init_met'), the program will produce new records for each
combination, e.g. if there are 3 'VAR_SEQ' features, two found in isoform B and
one found in isoform C, and one 'CONFLICT' feature, the program will produce 
six records, i.e. a record for the "standard" isoform, a record for isoform B,
a record for isoform C, and an additional record for each splice isoform with
the reported conflict.  However, if a given combination is impossible (for
example, a polymorphism occurs in an exon absent in one isoform), an entry is
not printed to file and an error message is printed instead.

OUTPUT  

If neither the 'varseq', 'variant', or 'conflict' options are specified,
or if only the 'varseq' option is specified, output is produced as specified in
the README.Otherwise, output is produced as specified below.

AC NUMBERS (extra options only)

AC numbers for each new entry are constructed as follows: 

parental_AC_number
-varsplic_number
-variant_number
-conflict_number
-init_met_number

such that P12345-00-00-00-00 would possess the same sequence as the parent
record, and P12345-01-02-04-00 would possess the splicing variations belonging
to the first alternative splice form, the variant features belonging to the
second alternative variant form, and the conflicts belonging to the fourth
alternative conflict form, starting from the default initiator methionine.

in the interests of restricting the size of these surrogate AC numbers, each 
component of the complete identifier is only included if it corresponds to a
specified option.  So if the user specifies any two of the varsplic, variant,
conflict and init_met options, the program will generate surrogate IDs such as
P12345-01-02, whereas if the user specifies any three options, the program will
generate surrogate IDs such as P12345-01-02-03.

DE LINES (extra options only)

If output is produced in pseudo UniProt format, the DE line is constructed
as follows: 
Isoform NAME; Variant NAME; Conflict NAME; Init met NAME; from ENTRY_AC
e.g.
Isoform Long; Variant Displayed; Conflict REF. 3;  Init met 1 from 
P78314

Any annotated note associated with one particular isoform in the parent
entry can be added to the DE line of the appropriate child entry by use of the
note option, as described above.

STATISTICS

These are produced separately for each type of feature; if multiple options
are specified, the statistics for each type of feature are the same that would
have been produced had only that option been specified .

e.g. if there is one alternatively spliced isoform documented in the parental
record, the record will be recorded as describing a protein with 2 
isoforms, even if each alternative splicing pattern has been combined with each
of several polymorphisms to generate many more alternate sequences.

=head1 AUTHOR/CONTACT

B<Paul Kersey>, pkersey@ebi.ac.uk

=head1 VERSION - 2.2.7

=cut

# ** Initialisation

use vars qw($opt_input 
            $opt_error $opt_stats $opt_fasta $opt_pseudo $opt_statsfile 
            $opt_threshold $opt_which $opt_linelength 
            $opt_varsplic $opt_variant $opt_conflict $opt_varseq 
            $opt_initiator
            $opt_uniqids $opt_uniqspids $opt_showdesc $opt_dbcode $opt_describe 
            $opt_note 
            $opt_filter $opt_count $opt_noftids 
            $opt_crosscheck $opt_check_vsps
            $opt_locations);

use strict;
use SWISS::Entry;
use SWISS::CC;
use SWISS::DE;
use Getopt::Long;

use Data::Dumper;

&GetOptions("input=s", 
            "error=s", "stats=s", "fasta:s", "pseudo:s", "statsfile=s", 
            "threshold=i", "which=s",  "linelength=s",
            "varsplic", "variant", "conflict", "varseq",
            "initiator",
            "uniqids:i", "uniqspids", "showdesc", "dbcode", "describe=s",
            "note",
            "filter:s", "count", "noftids",
            "crosscheck", "check_vsps", "locations");
            
# global variables

my @featureTypes; # expand these types of features
my %featureTypes; # hash of same data
my @range; # stores the ranges for which statistics are to be returned
my $records; # counts the number of records parsed

# hashes for cross-reference checking

my %xrefs;
my %checksums;

my $globalHash; # for sorting

my %allCalc;
my %noChange;
my %maxNumberOfIsoforms;
my %numberOfIsoforms;
my %errors;
my $newFormCount;

# global variables for format options

my $i;
my $featureType;
my $percent;

# opening input/output files

if (defined $opt_input)  {
    
  open (IN, "$opt_input") || die "Could not open $opt_input";
  
} else { 
    
   *IN = *STDIN;
} 

if ($opt_error) {
  open (ERROR, ">$opt_error") || die "Could not open $opt_error";
}

if ($opt_stats) {
  open (STATS, ">$opt_stats") || die "Could not open $opt_stats";
}

my $pseudoF = 0;

if (defined $opt_pseudo) {
  
  if ($opt_pseudo =~ /\w/) {
    
    open (PSEUDO, ">$opt_pseudo") || die "Could not open $opt_pseudo";
  
  } else { 
    
     *PSEUDO = *STDOUT;
     $opt_pseudo = "yes";
     $pseudoF = 1; 
  }
}

if (defined $opt_fasta) {
  
  if ($opt_fasta =~ /\w/) {
  
    open (FASTA, ">$opt_fasta") || die "Could not open $opt_fasta";
  
  } elsif ($pseudoF != 1) { 
  
     *FASTA = *STDOUT;
     $opt_fasta = "yes";
  
  } else {
  
    die "To produce both FASTA and pseudo-UniProt format output, you MUST " .
        "specify a destination file for at least one type of output.  The " .
        "other will default to STDOUT if not specified.\n";
  }
}

if ($opt_describe) {
  open (DESCR, ">$opt_describe") || die "Could not open $opt_describe";
}

# statsfile: this option allows the user to configure the statistics they 
# want to see

if ($opt_statsfile) {
  
  open (STATSFILE, "$opt_statsfile") || die "Could not open $opt_statsfile";
  my $q = 0;
  
  while (<STATSFILE>) {
  
    chomp $_;
    
    if (($_ > 0) && ($_ < 100)) {
    
      my $term = $q . "-". $_;
      push @range, $term;
      $q = $_;
      
    } else {
    
      die "Invalid terms in statistics input file";
    }
  }
  
  my $term = $q . "-100";
  push @range, $term;

} else {

  @range = qw(0-1 1-2 2-5 5-10 10-20 20-50 50-100);
}

# wrapping of FASTA file, preferred linelength will usually be >= 60

if ($opt_linelength) {

  if ($opt_linelength < 60) {

    print STDERR "Warning! Short Linelength specified for FASTA output " .
                 "file\n";
  }
}

# ensure compatability with new and old option names

if ($opt_varseq) {

  $opt_varsplic = 'Y';
}

# perform varsplic expansion by default

if ((! $opt_varsplic) && (! $opt_variant) && (! $opt_conflict)) {

  push @featureTypes, 'VAR_SEQ';

} else {

  if ($opt_varsplic) {
  
    push @featureTypes, 'VAR_SEQ';
  } 
  
  if ($opt_variant) {
  
    push @featureTypes, 'VARIANT';
  }
  
  if ($opt_conflict) {
  
    push @featureTypes, 'CONFLICT';
  }
}

my $varsplicOnly = 0;

if ((scalar @featureTypes == 1) && ($featureTypes[0] eq 'VAR_SEQ')) {

  $varsplicOnly = 1;
}

for my $featureType (@featureTypes) {

  $featureTypes{$featureType} = "yes";
}

# comment type

my %commentType;

# can use two alternative techniques for processing VARSPLIC data

if (! $opt_noftids) {

  $commentType{'VAR_SEQ'} = 2;

} else {

  $commentType{'VAR_SEQ'} = 1;
}

$commentType{'VARIANT'} = 1;
$commentType{'CONFLICT'} = 1;

# type order, as hash and list

my %typeOrder;
$typeOrder{'VAR_SEQ'} = 0;
$typeOrder{'VARIANT'} = 1;
$typeOrder{'CONFLICT'} = 2;

my @deName;
$deName[$typeOrder{'VAR_SEQ'}] = "Isoform";
$deName[$typeOrder{'VARIANT'}] = "Variant";
$deName[$typeOrder{'CONFLICT'}] = "Conflict";


my @specifications = sort _byTypeOrder keys %typeOrder;

# set default threshold to zero

unless ($opt_threshold)  {
  
  $opt_threshold = 0;
}

# read an entire record at a time

$/ = "\n\/\/\n";

my $records = 0;

ENTRY: while (<IN>) {

  # re-initialise variables for each new record

  my %isoformNames;
  my %isoform2id;
  my %isoform2note;
  my %count;
  my %crcHash;
  my $count;
  
  my %featId2start = ();
  my %featId2end = ();
  my %featId2desc = ();
  my %featId2featureType = ();
  my %isoform2feature;
  my %isoformNames;
  my %isoform2id;
  
  # intialisation
  
  @{$isoformNames{'VAR_SEQ'}} = ();
  @{$isoformNames{'VARIANT'}} = ();
  @{$isoformNames{'CONFLICT'}} = ();
  
  @{$isoform2feature{'VAR_SEQ'}} = ();
  @{$isoform2feature{'VARIANT'}} = ();
  @{$isoform2feature{'CONFLICT'}} = ();
  
  @{$isoform2id{'VAR_SEQ'}} = ();
  @{$isoform2id{'VARIANT'}} = ();
  @{$isoform2id{'CONFLICT'}} = ();
  
  %{$isoform2note{'VAR_SEQ'}} = ();
  
  my $featId = 0;
  
  $records++;
  
  # read in entry
  
  $_ =~ s/\n\*\*.*//gm if $opt_filter;
  my $entry = SWISS::Entry->fromText($_);
  my $ac = ${$entry -> ACs -> list()}[0];
  my $isHuman;
  my $addMethionine = 0;
  
  if ($opt_initiator) {
  
    my @initMets = $entry -> FTs -> get("INIT_MET");
    
    for my $initMet (@initMets) {
    
      if (($$initMet[1] == 0) && ($$initMet[2] == 0)) {
      
        $addMethionine = 1;
      }
    }
  }
  
  # code to be activated when CC DISEASE clean up is complete
  
  # foreach my $taxid ($entry->OXs->NCBI_TaxID()->elements()) {
    
  #   if ($taxid -> text() eq '9606') {
     
  #     $isHuman = "yes";
  #   }
  # }

  my $needCheck = 0;

  if ($opt_crosscheck) {
  
    $needCheck = _needCrossCheck($entry, \%xrefs);
  }

  for my $featureType (@featureTypes) {
  
    # treat VARSPLIC features differently to other types of feature now that 
    # CC blocks have been constructed

    if ($commentType{$featureType} == 1) {
    
      # standard code for dealing with variants and conflicts
      
      ($featId, $count) = _processOldCommentBlock($entry, 
                                                  $featureType,
                                                  $count, 
                                                  $ac,
                                                  $featId,
                                                  \%isoformNames, 
                                                  \%isoform2feature,
                                                  \%featId2start,
                                                  \%featId2end,
                                                  \%featId2desc,
                                                  \%featId2featureType,
                                                  \%isoform2note,
                                                  $isHuman);
      
    } elsif ($commentType{$featureType} == 2) {
    
      # special code for dealing with alterntive splicing
    
      ($featId, $count) = _processNewCommentBlock($entry, 
                                                  $featureType, 
                                                  $count, 
                                                  $ac,
                                                  $featId,
                                                  \%isoform2id,
                                                  \%isoformNames,
                                                  \%isoform2feature,
                                                  \%featId2start,
                                                  \%featId2end,
                                                  \%featId2desc,
                                                  \%featId2featureType,
                                                  \%isoform2note); 
    } 
  }
  
  # holding pattern
  
  if (! $opt_variant) {
  
    ${$isoformNames{'VARIANT'}}[0] = 'Displayed';
    ${$isoform2id{'VARIANT'}}[0] = $ac;
  }
  
  if (! $opt_conflict) {
  
    ${$isoformNames{'CONFLICT'}}[0] = 'Displayed';
    ${$isoform2id{'CONFLICT'}}[0] = $ac;
  }
  
  if ((! $opt_varsplic) && ($opt_variant || $opt_conflict)) {
   
    ${$isoformNames{'VAR_SEQ'}}[0] = 'Displayed';
    ${$isoform2id{'VAR_SEQ'}}[0] = $ac;
  }
  
  # are there any new isoforms to work with?
  
  if ($count == 0) {
      
    # if not, don't print anything unless $opt_which has been specified as 
    # "full"
      
    # special output procedures for records with no variants, output files are
    # standard (i.e. no version numbers, warning messages etc.)
     
    if ($opt_which eq "full") {
      
      my $entry2;
      
      if ($addMethionine > 0) {
        
        $entry2 = SWISS::Entry -> new();
        $entry2 -> SQs -> seq("M" . $entry -> SQs -> seq());
        $entry2 -> ACs -> list([${$entry -> ACs -> list()}[0]]);
        $entry2 -> IDs -> list([${$entry -> IDs -> list()}[0]]);
        $entry2 -> DEs -> text(${$entry -> DEs -> elements()}[0] -> text);
        _makeEntry2fromEntry($entry, $entry2);
        
        # CC text: should always adopt this form in this case
        # (hence not generalised in method)
        
        my $CCtext = "This entry " .
                     "represents an alternative sequence " . 
                     "derived from UniProtKB entry " . $ac . 
                     ". Please do not " .
                     "distribute independently of the parent UniProtKB entry";
        my $newCC = SWISS::CC -> new();
        $newCC -> topic("WARNING");
        $newCC -> comment($CCtext);
        $entry2 -> CCs -> add($newCC);
        $entry2 -> update();
      
      } else {
      
        $entry2 = $entry;
      }
      
      if ($opt_pseudo)  {
        
        print PSEUDO $entry2 -> toText();
      }
        
      if ($opt_fasta)  {
      
        _prepareFastas($entry, $entry2, -1);
      }
    }
      
    # increment relevent stats, then move on to next entry
        
    for my $specificationType (@specifications) {
      
      ${$numberOfIsoforms{$specificationType}}[0]++;
    }
    
    if ($needCheck == 1) {
            
      my %name2crc;
      
      $name2crc{$isoformNames{'VAR_SEQ'}[$i]} = 
        $entry -> SQs -> crc(); 
        
      push @{$checksums{$ac}}, \%name2crc;
    }
     
    next ENTRY;
  }
  
  # next, collect all features and sort these
  
  $globalHash = \%featId2start;
  my @featIds = sort _byStartPosition keys %featId2start;
  
  ## Now effect appropriate changes for each isoform
  
  # each isoform can be defined created by a certain combination of VARSPLIC,
  # VARIANT and CONFLICT features
  # but some features of the different types are incompatible
  # and certain information is known/assumptions can be made about the
  # co-occurance of some instances of the same type
    
  # start by storing sequence length
    
  my $oldSize = $entry -> SQs -> length();
  
  # intitialise current isoform
  
  my @currentIsoform;
  my %versionTypes;
  
  for (my $i = 0; $i < scalar @featureTypes; $i ++) {
  
    $currentIsoform[$i] = 0;
    $versionTypes{$featureTypes[$i]} = \$currentIsoform[$i];
  }
  
  # according to user-specified options, do sequence expansion for whatever
  # combination of features is requested
  
  my $i = 0;
  
  do {
   
    
    my $parent = 0;
    my $skip = 0;  
    my $noInitMet = 0;      
        
    # ignore isoforms without details 
    # (currently only possible for VARSPLIC)
    
    if 
   (${${$isoform2feature{$featureTypes[$i]}}[$currentIsoform[$i]]}{"Ignore"}) {
     
      $skip = 1;
    }    
    
    if ($skip == 0) {  
      
      # when feature type is VARSPLIC only and noftids option is not used, 
      # otherwise, don't display master form
      
      #  otherwise, don't display first form
      
      if (($featureTypes{'VAR_SEQ'}) &&
          (! $featureTypes{'VARIANT'}) &&
          (! $featureTypes{'CONFLICT'}) &&
          (${${$isoform2feature{'VAR_SEQ'}}[$currentIsoform[0]]}{"Master"} 
            eq "yes") &&
          (! $opt_noftids)) { 
         
        $parent = 1;  
        
      } else {
      
        $parent = 1;
      
        TYPE: for (my $j = 0; $j < scalar @featureTypes; $j ++) {
        
          if ($currentIsoform[$j] != 0) {
            
            $parent = 0;
            last TYPE;
          }
        }
      } 
        
      # processing parent entry
        
      if ($parent == 1) {
         
        if (! $opt_which) {
          
          if ($needCheck == 1) {
          
             my %name2crc;
             $name2crc{$isoformNames{'VAR_SEQ'}[$i]} = $entry -> SQs -> crc();
             push @{$checksums{$ac}}, \%name2crc;
          }

          $skip = 1;
        }
      }
    }     
       
    if ($skip == 0) {
     
      my $descriptionText;
        
      # $entry2 is initially the same as $entry
      # Changes will be made to $entry2, while $entry will be kept as 
      # reference
        
      my $entry2 = SWISS::Entry -> new();
      my $sq = $entry -> SQs -> seq();
      $entry2 -> SQs -> seq($sq);        
        
      # Reset variables for new form
        
      my @changedCoOrdinates;
      my $sizeChange = 0;
      my $unMatched = 0;
      my $totalAdd = 0;
      my $check = 0;
      my $percent = 0;
      my $startChange = 0;
      my $stopChange = 0;
        
      # for each feature
      #  what type of feature is it
      #  does the current form possess the current feature? 
        
      # we have one big feature list with features of all types on it
          
      # each isoform has a separate list of what features of which type it
      # should contain
        
      # i.e. if the current feature is a VARSPLIC feature, one needs to 
      # look up the VARSPLIC specification of the current isoform and see
      # if it contains this feature
        
      ## implmentation
      
      # Consider each feature at a time
        
      FEATURE: for (my $l = 0; $l < (scalar @featIds); $l++) {
      
        # does the current isoform contain the current feature?
        
        my $thisFeatId = $featIds[$l];
        my $thisFeatureType = $featId2featureType{$thisFeatId};
        my $releventIncrementor = ${$versionTypes{$thisFeatureType}};
         
        if 
   (${${$isoform2feature{$thisFeatureType}}[$releventIncrementor]}{$thisFeatId} 
          eq "yes") {
            
          # data integrity checking
        
          if (($featId2start{$thisFeatId} !~ /^\d+$/) ||
              ($featId2end{$thisFeatId} !~ /^\d+$/)) {
              
            # reject all isoforms containing a feature without clean,
            # numercial co-ordinates
              
            _errorMessage(\@currentIsoform, 
                          $ac, 
                          \%isoformNames,
                          \@featureTypes,
                          \@deName);
            $skip = 1;
          }
          
          if ($skip == 0) {
            
            #  Is this variation compatible with previous variations?
            
            CHECKCORD: for my $changedCoOrdinate (@changedCoOrdinates) {
              
              if ((($featId2start{$thisFeatId} >= 
                    ${$changedCoOrdinate}[0]) &&
                   ($featId2start{$thisFeatId} <= 
                    ${$changedCoOrdinate}[1])) ||
                  (($featId2end{$thisFeatId} >= ${$changedCoOrdinate}[0]) &&
                   ($featId2end{$thisFeatId} <= ${$changedCoOrdinate}[1]))) {
                
                _errorMessage(\@currentIsoform, 
                              $ac, 
                              \%isoformNames,
                              \@featureTypes,
                              \@deName);
                $skip = 1;
                last CHECKCORD;
              }
            }
          }
          
          my ($add, $remove);
          
          if ($skip == 0) {
            
            # add region to be altered now onto list of forbidden regions in
            # future
            
            my @coOrdinatePair = ($featId2start{$thisFeatId},
                                  $featId2end{$thisFeatId});
            push @changedCoOrdinates, \@coOrdinatePair;
              
            # identify what sort of splice variant this is
             
            if ($featId2desc{$thisFeatId} =~ /^Missing/i)  {
              
              # if it is a deletion there is nothing to add
                
              $add = "";
          
            } else {
              
              # if it is a substitution we need to dissect the string which 
              # tells us what to add
              
              # first, must tidy up string
              
              $featId2desc{$thisFeatId}=~ s/[ \n]//g;
              
              unless (($remove, $add) = 
                      ($featId2desc{$thisFeatId} =~ /(\w+)->(\w+)/)) {
                
                _errorMessage(\@currentIsoform, 
                              $ac, 
                              \%isoformNames,
                              \@featureTypes,
                              \@deName);
                $skip = 1;
              }
            }
          }
          
          if ($skip == 0) {
          
            # description of details of proposed changes
              
            if ($opt_describe) {
              
              $descriptionText = $descriptionText . 
                                 $featId2start{$thisFeatId} . "," . 
                                 $featId2end{$thisFeatId} . "," .
                                 (length $add) . "," . 
                  ($featId2end{$thisFeatId} - $featId2start{$thisFeatId} + 1) . 
                                 ";";                
            }
           
            # call routine to delete missing sequence and add new 
            # sequence where appropriate
            
            # get sequence and make appropriate changes to it
              
            my $sq = $entry2 -> SQs -> seq();
            
            if ($featId2start{$thisFeatId} == 1) {
            
              $noInitMet = 1;
            }
            
            my @ret = _seqCut($sq, 
                              $featId2start{$thisFeatId}, 
                              $featId2end{$thisFeatId}, 
                              $add, 
                              $sizeChange);
                               
            if ($sq !~ /\w/) {
             
              _errorMessage(\@currentIsoform, 
                            $ac, 
                            \%isoformNames, 
                            \@featureTypes,
                            \@deName,
                            "Too big!");
              $skip = 1;
            
            } else {
            
              # record extremities of chnages (using co-ordinates of parent
              # sequence)
            
              if ($startChange == 0) {
            
                $startChange = $featId2start{$thisFeatId};
              }
            
              $stopChange = $featId2end{$thisFeatId};
              
              # record any change in size effected by the editing 
              # (subsequently, will need to offset further changes by this
              # ammount from their original sequence co-ordinates) 
          
              # a cleverer approach: implement changes backwards first?
          
              $sizeChange = $ret[2]; 
          
              # keep track of total number of unMatched bases, total number of 
              # bases added
              
              $unMatched += $ret[3];
              $totalAdd += $ret[4]; 
          
              # if we have done a replacement not a cut, we can check that it
              # was correct
              
              if ($add ne "")  {
              
                if ($remove ne $ret[1])  {
                
                  $check = 1;
                }  
              }
            
              ## alter the sequence record
          
              $entry2 -> SQs -> seq($ret[0]);
          
              if ($ret[0] eq '') {
              
                _errorMessage(\@currentIsoform, 
                              $ac, 
                              \%isoformNames,
                              \@featureTypes,
                              \@deName);
              }
            }
          }
        
          # print warning if FT data is not consistent
          # may result from data error or simple incompatability of two
          # expanded features, i.e. it doesn't mean anything is wrong, simply
          # that this so-called isoform can't actually exist and should be
          # ignored
        
          if ($check == 1)  {
          
            _errorMessage(\@currentIsoform, 
                          $ac, 
                          \%isoformNames, 
                          \@featureTypes,
                          \@deName);
            $skip = 1;
          }
        }
      }
      
      if ($skip == 0) {
        
        # correction to sequence indicated by INIT MET annotation
        
        if (($addMethionine > 0) && ($noInitMet == 0)) {
        
          $entry2 -> SQs -> seq("M" . $entry2 -> SQs -> seq());
        }
        
        # syntax check on sequences that don't start with methionine
        
        if ($entry2 -> SQs -> seq !~ /^M/) {
        
          if (!
           _noMethionineOk($entry,
                          ${$isoformNames{'VAR_SEQ'}}[$currentIsoform[0]])) {
        
            # print error message but still produce sequence in this case
        
            _errorMessage(\@currentIsoform, 
                          $ac, 
                          \%isoformNames, 
                          \@featureTypes,
                          \@deName,
                          "No methionine!");
          } 
        }
       
        # have now completed sequence calculations for one form
        # next: display results
        
        # only want to print out record if % change exceeds threshold value
        # this is calculated as follows:
        # a gapped alignment is assumed, aligning all unchanged bases
        # the percentage of unalinged positions in this alignment is then 
        # calculated
        # this is actually done on the fly, now we only need to do a final sum
      
        my $calc = $unMatched/($oldSize + $totalAdd);
        
        if ($calc >= ($opt_threshold/100))  {
        
          my $de = "";
          my $isoformAc = "";
          my $crc = $entry2 -> SQs -> crc();
          
          if ($needCheck == 1) {
            
            my %name2crc;
            $name2crc{${$isoformNames{'VAR_SEQ'}}[$currentIsoform[0]]} = $crc;
            push @{$checksums{$ac}}, \%name2crc;
          }
          
          if ($crcHash{$crc}) {
                 
            if ($varsplicOnly == 1) {
            
              print ERROR "The following variants of AC " . $ac . 
                           " have the same sequences: " . 
                           ${$isoformNames{'VAR_SEQ'}}[$currentIsoform[0]] . " and " . 
                           $crcHash{$crc} . "\n";
                           
            } else {
               
              print ERROR "The following variants of AC " . $ac . 
                          " have the " .
                          "same sequences: ";
                
              my $name = _generateIsoformName(\@featureTypes, 
                                              \@currentIsoform,
                                              \%isoformNames,
                                              \@deName);
              print ERROR $name . " and " . $crcHash{$crc} . "\n";
            }
           
          } elsif ($varsplicOnly == 1) {
          
            $crcHash{$crc}= $isoformNames{'VAR_SEQ'}[$i];
          
          } else {
          
            my $name = _generateIsoformName(\@featureTypes, 
                                            \@currentIsoform,
                                            \%isoformNames,
                                            \@deName);
            $crcHash{$crc}= $name;
          }  
          
          # FASTA file is produced from UniProt file
          # therefore produce UniProt file if either form of output is 
          # necessary
          
          if ($opt_pseudo || $opt_fasta) {
          
            _preparePseudoEntry($entry, 
                                $entry2,
                                \@currentIsoform, 
                                $ac,
                                \%isoform2id, 
                                \%isoformNames,
                                \@featureTypes,
                                \%isoform2note,
                                $varsplicOnly,
                                $startChange,
                                $stopChange);
            
            # finish description (also needs data in UniProt format)
          
            if ($opt_describe) {
        
              $descriptionText = _doDescription($entry2, 
                                                \@currentIsoform, 
                                                \%isoform2id) . 
                                 $descriptionText;
                                 
              print DESCR $descriptionText . "\n";
            }
            
            if ($opt_pseudo) {
          
              print PSEUDO $entry2 -> toText();
            }
          
            if ($opt_fasta)  {
             
              _prepareFastas($entry, 
                             $entry2, 
                             $parent, 
                             $startChange, 
                             $stopChange);
            } 
          }                  
        }        
        
        # prepare statistics for this variant, unless it is a parent!
        # do separate statistics for isoforms, variants, conflicts
        
        # only do statistics for each factor singly, not in combination
        # i.e. only do these statistics if only 1 factor (i.e. varsplic, 
        # variant, conflict) differs from the parent 
        
        # if this is not the parent entry
        
        if ($parent != 1) {
          
          my $procede = 'N';
          my @zeroes = _findZeroes(\@featureTypes, \@currentIsoform);
          my $target; 
               
          if 
  (${${$isoform2feature{'VAR_SEQ'}}[$currentIsoform[0]]}{"Master"} eq "yes") { 
              
            $target = scalar @featureTypes - 2;
              
          } else {
              
            $target = scalar @featureTypes - 1;
          } 
              
          if (scalar @zeroes >= $target) {
            
            $procede = 'Y';
          }
           
          if ($procede eq 'Y') { 
               
            if ($calc == 0) {
           
              my $isoformAc = _getAc(\@currentIsoform, $ac, \%isoform2id);
                print ERROR "WARNING - Apparently unchanged isoform: ", 
                            $isoformAc, "\n";
            }
              
            for my $specificationType (@specifications) {
            
            # if this isoform is not the displayed form in respect of a 
            # given criteria for expansion
            
              my $thisIncrementor = 0;
                
              if ($versionTypes{$specificationType}) {
                
                $thisIncrementor = ${$versionTypes{$specificationType}};
              }  
                
              if (($thisIncrementor > 0) ||
                  (($specificationType eq 'VAR_SEQ') &&
                   ($featureTypes{'VAR_SEQ'} eq "yes") &&
              (${${$isoform2feature{'VAR_SEQ'}}[$currentIsoform[0]]}{"Master"} 
                    ne "yes"))) {
            
                # if this isoform is not the displayed VARSPLIC form,
                # and VARSPLIC is a criteria for expansion, then need to
                # calculate data for calculation of median
                   
                $count{$specificationType}++;   
                $calc *= 100;
                push @{$allCalc{$specificationType}}, $calc;
              
                # assign a category according to the degree of change
              
                ASSIGN: for (my $q = 0; $q < (scalar @range); $q++) {
              
                  my ($from, $to) = split /-/, $range[$q];
                
                  if (($calc > $from) && ($calc <= $to)) {
                
                    ${$noChange{$specificationType}}[$q]++;
                    last ASSIGN;
                  }
                }
              }
            }
          } # end if $procede
        } # end if $parent
      } #end if $skip
    }
    
    if ($currentIsoform[$i] < 
        ((scalar @{$isoformNames{$featureTypes[$i]}}) - 1)) {
    
      # if there are still more isoforms described using the current key, 
      # continue to look at this isoforms caused by this key type
      
      $currentIsoform[$i] ++;
    
    } else {
    
      # or, prepare to increment the next key type (fi extant) instead
    
      for (my $j = 0; $j <= $i; $j ++) {
      
        $currentIsoform[$j] = 0;
      }
        
      my $newI = $i;
      my $oldI = $i; 
      
      while ($i == $oldI) {
      
        if (defined $currentIsoform[$newI + 1]) {
          
          if ($currentIsoform[$newI + 1] < 
              scalar @{$isoformNames{$featureTypes[$newI + 1]}} - 1) {
      
           $i = $newI + 1;
           
          } else {
            
            $newI = $newI + 1;
          }  
          
        } else {
          
          $i = $newI + 1;
        }
      }
        
      for (my $j = 0; $j < $i; $j ++) {
      
        $currentIsoform[$j] = 0;
      }
      
      if ($i < scalar @featureTypes) {
      
        $currentIsoform[$i] ++;
        $i = 0;
      }
        
    }
    
  } until ($i >= scalar @featureTypes);
  
  # have now processed all forms of entry

  # Record maximum number of variants per protein
  
  for my $specificationType (@specifications) {
                 
    if ($count{$specificationType} > $maxNumberOfIsoforms{$specificationType})  {
      
      $maxNumberOfIsoforms{$specificationType} = $count{$specificationType};
    }
    
    # keep count of the number of records with this number of isoforms
    # note: total number of isoforms or additional isoforms;
    
    ${$numberOfIsoforms{$specificationType}}[$count{$specificationType}]++;
  }
}

if ($opt_pseudo)  {

  close PSEUDO;
}

if ($opt_fasta)  {

  close FASTA;
}
  
if ($opt_stats)  {

  _prepareStats();
}
  
if ($opt_crosscheck) {
 
   _crossCheck(\%xrefs, \%checksums);
}
 
if ($opt_filter) {

  system("rm $opt_filter");
}

## subroutines 

sub _byNumber {

  $a <=> $b;
}

sub _byStartPosition {
  
  $$globalHash{$a} <=> $$globalHash{$b};
}

sub _byTypeOrder {

  $typeOrder{$a} <=> $typeOrder{$b}
}

sub _crossCheck {

  my ($xrefs, $checksums);
  my @acs = keys %xrefs;
  
  for my $ac (@acs) {
  
    my %errorText = "";
    my $majorError = "";
    my %theseChecksums;
    my %name2Checksum;
    
    for my $checksumHash (@{$checksums{$ac}}) {
    
      my @names = keys %{$checksumHash};
      
      for my $name (@names) {
      
        my $checksum = $$checksumHash{$name};
        push @{$theseChecksums{$checksum}}, $name;
        $name2Checksum{$name} = $checksum;
      }
    }
    
    my @xrefs = keys %{$xrefs{$ac}};
    
    for my $xref (@xrefs) {
      
      for my $checksumHash (@{$checksums{$xref}}) {
        
        my @names = keys %{$checksumHash};
        
        for my $name (@names) {
        
          my $checksum = $$checksumHash{$name}; 
          
          # check that this entry's account of the checksum of each
          # variant doesn't overlap with another entry
        
          if (($name2Checksum{$name}) && 
              ($name2Checksum{$name} ne $checksum)) {
           
           $majorError = $majorError .
                         "Problem with entry $ac and associates: " .
                         "2 different sequences for isoform $name\n";
             
        
          } else {
        
            $name2Checksum{$name} = $checksum;
          }
          
          if ($majorError != 1) {
          
            if ($theseChecksums{$checksum}) {
        
              # overwitre previous error message for this checksum
              # write a new message now an extra entry has been added to this
              # cluster
        
              $errorText{$checksum} = "Problem with entry " . $ac . 
                                     " and asociates: (Isoforms " . $name;
            
              for my $otherName (@{$theseChecksums{$checksum}}) {
            
                $errorText{$checksum} = $errorText{$checksum} . ", " . 
                                        $otherName;
              }
            
              $errorText{$checksum} = $errorText{$checksum}  . 
                                      " have duplicate sequence, checksum: " .
                                      $checksum . ")\n";
        
            } 
        
            push @{$theseChecksums{$checksum}}, $name;
          }
        }
      }
    }
    
    if ($majorError eq "") {
    
      my @checksums = keys %errorText;
      
      for my $checksum (@checksums) {
    
        print ERROR $errorText{$checksum};
      }
    
    } else {
    
      print ERROR $majorError;
    } 
  }
}

sub _doDescription {
  
  my ($entry2, $currentIsoform, $isoform2id) = @_;
  my @currentIsoform = @$currentIsoform;
  my %isoform2id = %$isoform2id;
  my $ac = ${$entry2 -> ACs -> list}[0];
  my $parentAc = $ac;
  $parentAc =~ s/\-.*//;
  my $text = $parentAc . ":";
       
  if ((! $opt_variant) && (! $opt_conflict)) {
            
    $text = $text . $ac; 
      
  } else {
          
    $text = $text + _getAc(\@currentIsoform, $ac, \%isoform2id);
  }
          
  $text = $text . ":";
  return $text;
}

sub _errorMessage {

  # subroutine to print out a message about incompatible specifications
  
  my ($currentIsoform, $ac, $isoformNames, $featureTypes, $deName, $message) 
    = @_;
    
  my %isoformNames = %$isoformNames;
  my @currentIsoform = @$currentIsoform;
  my @featureTypes = @$featureTypes;
  my @deName = @$deName;
  print ERROR "Data concerning ";
  
  print ERROR _generateIsoformName(\@featureTypes, 
                                   \@currentIsoform, 
                                   \%isoformNames,
                                   \@deName);
                                   
  print ERROR " of " . $ac . " inconsistent.";
  
  if ($message =~ /\w/) {         
   
    print ERROR " " . $message;
  
  } else {
  
    print ERROR " No entry output to file for this variant.";
  }
  
  print ERROR "\n";
   
  my @zeroes = _findZeroes(\@featureTypes, \@currentIsoform, "Y"); 
    
  if (scalar @zeroes == 1) {
            
    $errors{$featureTypes[$zeroes[0]]} ++;
  }
}

sub _findZeroes {
  
  my ($featureTypes, $currentIsoform, $inverse) = @_;
  my @featureTypes = @$featureTypes;
  my @currentIsoform = @$currentIsoform;
  my @zeroes;
              
  for (my $j = 1; $j < scalar @featureTypes; $j ++) {
               
    if (($currentIsoform[$j] == 0) && (! defined $inverse)) {
                 
      push @zeroes, $j;
    
    } elsif (($currentIsoform[$j] != 0) && (defined $inverse)) {  
      
      push @zeroes, $j;
    }
  }
  
  return @zeroes;
}              

sub _getAc {

  # quick method for just getting AC line
  
  my ($currentIsoform, $ac, $isoform2id) = @_;
  my @currentIsoform = @$currentIsoform;
  
  if (($varsplicOnly == 1) && 
      (${${$isoform2id}{'VAR_SEQ'}}[$currentIsoform[0]])) {
    
    return ${${$isoform2id}{'VAR_SEQ'}}[$currentIsoform[0]];

   } else {

     for (my $i = 0; $i < scalar @currentIsoform; $i ++) {

       my $nextBit = $currentIsoform[$i];

       if ($nextBit < 10) {
    
         $nextBit = "0" . $nextBit;
       }
      
       $ac = $ac . "-" . $nextBit;
    }

    return $ac;
  }
}

sub _getDe {

  # quick method for just getting DE line
  
  my ($currentIsoform, $ac, $isoformNames, $featureTypes, $deName, 
      $isoform2note, $varsplicOnly) = @_;
  my @currentIsoform = @$currentIsoform;
  my @featureTypes = @$featureTypes;
  my %isoformNames = %$isoformNames;
  my %isoform2note = %$isoform2note;
  my $name;
  
  if ($varsplicOnly == 1) {
  
     my  $currentName = ${$$isoformNames{'VAR_SEQ'}}[$currentIsoform[0]]; 
     $name = "Isoform " . $currentName. " of " . $deName;
   
     if (${$isoform2note{'VAR_SEQ'}}{$currentName}) {
     
       $name = $name  . " (Note:" . 
               ${$isoform2note{'VAR_SEQ'}}{$currentName} . ")";
     }
  
  } else {
  
     $name = _generateIsoformName(\@featureTypes, 
                                  \@currentIsoform, 
                                  \%isoformNames,
                                  $deName);
     $name = $name . " of " . $ac;
     
  }
  
  return $name;
}

sub _getId {
    
  my ($entry) = @_;
  my @ids = $entry -> IDs -> list();
  my $newid = ${$ids[0]}[0];
  return $newid;
}

sub _generateIsoformName {

  my ($featureTypes, $currentIsoform, $isoformNames, $deName) = @_;
  my @featureTypes = @$featureTypes;
  my @currentIsoform = @$currentIsoform;
  my %isoformNames = %$isoformNames;
  my @deName;
  my $name = "";
  
  # 2 patterns of names: may or may not include a label
  
  for (my $j = 0; $j < scalar @featureTypes; $j ++) {            
     
    $name = $name . $deName[$typeOrder{$featureTypes[$j]}] . " " .
            ${$isoformNames{$featureTypes[$j]}}[$currentIsoform[$j]]; 
         
  }
  
  $name =~ s/\s+$//;
  return $name;               
}

sub _makeEntry2fromEntry {

  my ($entry, $entry2) = @_; 
  my $gn = $entry -> GNs();
  $entry2 -> GNs($gn);
  my @oc = $entry -> OCs -> list();
  $entry2 -> OCs -> list(@oc);
  my @og = $entry -> OGs -> list();
  $entry2 -> OGs -> list(@og);
  my @os = $entry -> OSs -> list();
  $entry2 -> OSs -> list(@os);
  $entry2 -> OXs -> NCBI_TaxID -> list($entry -> OXs ->NCBI_TaxID() -> list());
}      


sub _needCrossCheck {

  my ($entry, $xrefs) = @_;
  my $needCheck = 0;
  my $cc = 
      $entry -> CCs -> filter(&SWISS::CCs::ccTopic("ALTERNATIVE PRODUCTS"));
  
  if (defined ${$cc -> list}[0]) {
      
    my $block = ${$cc -> list}[0];
    my $keyEvent = $block -> keyEvent();
    my @isoformNames = $block -> getFormNames($keyEvent);
     
    for my $isoformName (@isoformNames) {
      
      my @ftIds = $block -> getFeatIds($keyEvent, $isoformName);
        
      if ($ftIds[0] eq "External") {
          
         my @isoIds = 
           $block -> getIsoIds($keyEvent, $isoformName);
         my ($externalAc) = $isoIds[0];
         $externalAc =~ s/\-\d+$//;
         my $ac = ${$entry -> ACs -> list}[0];
         ${$xrefs{$ac}}{$externalAc} = "yes"; 
         $needCheck = 1;
      }
    }
  }
  
  return $needCheck;
}

sub _noMethionineOk {

  my ($entry, $formName) = @_;
  my $FTs = $entry -> FTs;
  my @nonTers = $FTs -> get('NON_TER');
  
  for my $nonTer (@nonTers) {
  
    if (($$nonTer[1] == 1) && ($$nonTer[2] == 1)) {
      
      return 1;
    }
  }   
  
  my @initMets = $FTs -> get('INIT_MET');
  
  for my $initMet (@initMets) {
  
    if (($$initMet[1] == 0) && ($$initMet[2] == 0)) {
    
      return 1;
    }
  } 
  
  my @CCs = $entry -> CCs -> get("ALTERNATIVE PRODUCTS");

  for my $CC (@CCs) {

    my $keyEvent = $CC -> keyEvent();

    if ($CC -> getNote($keyEvent, $formName) =~
        /Incomplete sequence/) {
    
      return 1;
    }
  }
  
  return 0;
}

sub _parseFTs {

  # subroutine to parse FT lines with a stated key

  # $featId is an identifier for all relevent features from this program
  # $ftId is the actual FTId that is found in some (not all) UniProt 
  # features
  
  # later parameters are only needed when parsing VARIANT or CONFLICT lines 

  my ($entry, 
      $ftKey, 
      $featId, 
      $featId2start,
      $featId2end,
      $featId2desc,
      $featId2featureType) = @_;
  
  my %ftId2featId;
  
  # %isoform2feature is only used by VARIANT/CONLICT algorithm
  # VARSPLIC algorithm assigns isoforms in an outer method
  
  # for VARSPLIC we know which isoforms are comprised of which 
  # features from the CC lines, whereas for VARIANT and CONFLICT we have to 
  # get all this information from the FT lines
  
  # for INIT_MET, treat each feature indendently
  
  my %isoform2featId;
  my $isoformNumber = 1;
  
  if (my @fts = $entry->FTs->get($ftKey)) {
  
    FT: for my $ft (@fts)  {
      
      my ($keys, $startPos, $stopPos, $desc, $junk, $ftId) = @$ft;
      
      # tidy and edit
      
      $startPos =~ s/<//;
      $stopPos =~ s/>//;
      
      if ($ftKey eq 'INIT_MET') {
      
        # INIT MET at position 0 indicates alternative sequence of existing
        # isoforms, not the existence of additional isoforms 
      
        if ($startPos == 0) {
        
          next FT;
        }
        
        # treat other INIT_METs (at position n) as deletions from 1 to (n - 1)
        
        $startPos = 1;
        $stopPos = $stopPos - 1;
        $desc = "Missing"; 
      }
      
      $$featId2start{$featId} = $startPos;
      $$featId2end{$featId} = $stopPos;
      $$featId2desc{$featId} = $desc;
      $ftId =~ s/\/FTId\=//;
      $ftId =~ s/\.$//;
      $$featId2featureType{$featId} = $ftKey;
      
      if ($ftKey eq 'VAR_SEQ') {
      
        # set up various hashes
        # use according to method call
        
        $ftId2featId{$ftId} = $featId;
        $$featId2desc{$featId} = $desc;
          
        my ($realDesc, $isoforms) = ($desc =~ /(.+?)\(IN (.+)\)/i);
        my @isoforms  = split /\WAND\W|,\W(?!REF)/i, $isoforms;
          
        for my $isoform (@isoforms) {
       
          if ($isoform !~ /ISOFORMS/i) {
          
            if ($isoform !~ /^isoform/) {
            
              $isoform = "isoform " . $isoform;
              print ERROR $isoform . " of entry ". 
                          ${$entry -> ACs -> list()}[0] . " is referred " .
                          "to w/o 'isoform' prefix in FT lines\n";
            }
        
            push @{$isoform2featId{$isoform}}, $featId;
          }
        }  
        
      } elsif ($ftKey eq 'VARIANT') {
      
        my ($realDesc, $isoforms) = ($desc =~ /(.+?)\(IN (.+)\)/i);
        
        # description of variant is split up into a real description and an 
        # 'isoform' list
        
        if ($isoforms eq  '') {
        
          # if there is no isoform list
        
          # firstly, trim the junk from the real description
          
          ($realDesc) = ($desc  =~ /(.+)\/*/);
          $realDesc =~ s/\(.*//;
          
          # as there is no name for this variant, use FtId if possible, 
          # otherwise create one from the internal isoform count
          
          if ($ftId ne  '') {
          
            $isoforms = $ftId;
          
          } else {
          
            $isoforms = "Unnamed variant isoform ". $isoformNumber;
            $isoformNumber ++;
          }
          
          # next task is to set up an isoform - feature hash
        
          push @{$isoform2featId{$isoforms}}, $featId;
        
        } else {
        
          # if there is an isoform list, we need to parse these
          # each isoform may have a pesudo comment
          # if there is not a specific pseudo comment for a given isoform, 
          # assume that a terminal pseduo comment applies to all isoforms
          
          # the point of all this is that we want to keep/add pseudo-comments 
          # to VARIANT names where these are STRAINs; otherwise not
  
          my ($lastPseudoComment) = ($isoforms =~ /.+; (.*)/);
          $isoforms =~ s/(.+);.*/$1/;
          my @isoforms  = split /\WAND\W|,\W(?!REF)/i, $isoforms;
        
          for my $isoform (@isoforms) {
        
            my ($pseudoComment) = ($isoform =~ /.+?; (.*)/);
            $isoform =~ s/;.*//;
            $isoform = _tidy($isoform);
          
            # specific pseudo comment, if there is none, there may be a general
            # one to replace it with
          
            if ($pseudoComment eq '') {
          
              $pseudoComment = $lastPseudoComment;
            }

            if ($pseudoComment =~ /^STRAIN/) {
          
              $isoform = $isoform . "("  . $pseudoComment . ")";
            }
            
            push @{$isoform2featId{$isoform}}, $featId;
          }
        }
        
        # if there are a list of alternative variants in the one line, just
        # display the first
        
        if ($realDesc =~ / OR .*/i) {
      
          $realDesc =~ s/(.+?) OR .*/$1/i;
        }
        
        $$featId2desc{$featId} = $realDesc;
    
        # only HUMAN variants currently have FTIds
        # also, note that for VARIANTs, we are interested in the data the
        # opposite way around to that which interests us for VARSPLICs
        
        # therefore, use the same hash, but with reverse population
          
        if ($ftId ne "") {
          
          $ftId2featId{$featId} = $ftId;
        
        } 
        
      } elsif ($ftKey eq 'CONFLICT') {
      
        my ($realDesc, $isoforms) = ($desc =~ /(.+?)\(IN (.+)\)/i);
        my @isoforms = split /\WAND\W|,\W(?!REF)/, $isoforms;
       
        if ($realDesc =~ / OR /i) {
               
          $realDesc =~ s/(.+?) OR .*/$1/i;
        }
          
        $$featId2desc{$featId} = $realDesc;
        
        for my $isoform (@isoforms) {
        
          $isoform = _tidy($isoform);
          
          if ($isoform !~ /\bREF\b/i)  {
          
            if ($isoforms[0] =~ /^REF\b/i)  {
            
              $isoform = 'REF. '. $isoform;
            } 
          } 
          
          push @{$isoform2featId{$isoform}}, $featId;  
        }  
      
      } elsif ($ftKey eq 'INIT_MET') {
      
        # no attempt to use correct isoform names for INIT_MET related features
      
        my $isoform = "Initiator at " . ($stopPos + 1);
        push @{$isoform2featId{$isoform}}, $featId;
      }
      
      $featId ++;  
    }
  }
  
  return ($featId, \%ftId2featId, \%isoform2featId);
}

sub _prepareFastas {

  my ($entry, $entry2, $parent, $startChange, $stopChange) = @_;
  my ($newDE);
  my $f = $entry2-> toFasta();
  $f =~ s/PE\=.+?\n/\n/;
  
  # pseudo entries aren't recognised as Swiss-Prot (rather than TrEMBL) by
  # Swissknife FASTA dumper, so need to set data class manually.
  # (oddly, if pseudo entries have been printed, this behaviour disappears...)
  
  if ($entry2 -> IDs -> dataClass eq 'Reviewed') {
  
    $f =~ s/^\>tr/>sp/;
  }
  
  ## linelength
  
  if ($opt_linelength) {
  
    my ($f1, $f2);
    ($f1, $f2) = ($f =~ /^(\>.+? )(.*?)\n/);
    my $length = $opt_linelength - length $f1;
    ($newDE) = ($f2 =~ /(.{$length})/);
     
    if ($newDE eq '') {
              
     $newDE = $f2;
    }
        
    $newDE = $f1 . $newDE . "\n";
    $f =~ s/(^.+?\n)/$newDE/;
  }
  
  $f =~ s/ *$//gm; 
  print FASTA $f;    
}

sub _preparePseudoEntry {
          
  my ($entry, 
      $entry2,
      $currentIsoform,
      $ac, 
      $isoform2id, 
      $isoformNames,
      $featureTypes, 
      $isoform2note,
      $varsplicOnly,
      $startChange,
      $stopChange) = @_;
  my %isoform2id = %$isoform2id;
  my %isoformNames = %$isoformNames;
  my @currentIsoform = @$currentIsoform;
  my @featureTypes = @$featureTypes;
  my %isoform2note = %$isoform2note;
  
  # prepare entry in pseudo-UniProt format
          
  # note, we are doing this even wherethe pseudo option has not been 
  # set, wasteful
        
  # Set some attributes of $entry2 unchanged from $entry
     
  _makeEntry2fromEntry($entry, $entry2);
  
  # output: is determined by which options, if any, have been specified
          
  # case 1: only VARSPLIC has been (explicitly or implicitly) selected
  # case 2: non-VARSPLIC features have been selected
   
  my $isoformAc = _getAc(\@currentIsoform, $ac, \%isoform2id);
  my $newid = _getId($entry);
  
  if ($opt_uniqspids) {
  
      my ($suffix) = ($isoformAc =~ /(-.*)/);
      $newid = $newid . $suffix;
  }
  
  $entry2 -> ACs -> list([$isoformAc]);
  $entry2 -> IDs -> list([$newid]);
  $entry2 -> IDs -> dataClass($entry -> IDs -> dataClass());
  my $de = ${$entry -> DEs -> list}[0];
  my $recName = $de -> text();
  my @des;
  my $de = SWISS::DE -> new();
  $de -> category("IsoName");
  $de -> type("Full");
  
  my $text = _getDe(\@currentIsoform, $ac, \%isoformNames, \@featureTypes,
                    $recName, \%isoform2note, $varsplicOnly);
                    
  $de -> text($text);
  push @des, $de;
  
  if ($opt_locations) {
  
    my $de = SWISS::DE -> new();
    $de -> category("LocName");
    $de -> type("Full");
    $de -> text($startChange . "-" . $stopChange);
    push @des, $de;
  } 
  
  $entry2 -> DEs -> list(\@des);
  
  # change comments
  
  my $CCtext;
       
  if ($varsplicOnly == 1) {
          
    $CCtext = "This entry " .
              "represents splice isoform " .
              ${$isoformNames{'VAR_SEQ'}}[$currentIsoform[0]] .
              " derived from UniProtKB entry " . $ac . ". Please do not " .
              "distribute independently of the parent UniProtKB entry";
          
  } else {
          
    $CCtext = "This entry " .
              "represents an alternative sequence " . 
              "derived from UniProtKB entry " . $ac . ". Please do not " .
              "distribute independently of the parent UniProtKB entry";;
  }
  
  my $newCC = SWISS::CC -> new();
  $newCC -> topic("WARNING");
  $newCC -> comment($CCtext);
  $entry2 -> CCs -> add($newCC);
}
          
sub _prepareStats { 
 
  print STATS "STATISTICS\n";
  print STATS "----------";
  close STATS;
  
  for my $featureType (@featureTypes) {
  
    # initialise working variables
      
    my $totalVar = 0;
    my $totalCalc = 0;
      
    #  Start with stats for records

    open(STATS, ">>$opt_stats")  || die "Can't open $opt_stats";
    print STATS "\n\n", $featureType;
    print STATS "\n\nSTATISTICS FOR ORIGINAL RECORDS\n\n";
    print STATS "Total Number of Original Records: ", $records, "\n\n";
    print STATS 
    "Number of New Forms/Record      Number of Records       % of Records\n";
    my $percent;

    for (my $i = 0; $i <= $maxNumberOfIsoforms{$featureType}; $i++)  {
  
      ${$numberOfIsoforms{$featureType}}[$i] =~ s/^$/0/;
      
      $percent = (${$numberOfIsoforms{$featureType}}[$i] / $records) 
                  * 100;
      
      printf STATS "%13.13s", $i; 
      printf STATS "%28.28s", ${$numberOfIsoforms{$featureType}}[$i];
      printf STATS "%23.2f", $percent;
      print STATS "\n"; 
    
      # count total number of new variants
  
      $totalVar += ($i * ${$numberOfIsoforms{$featureType}}[$i]);
    }

    # Now do stats for new variants
      
    if (($totalVar - $errors{$featureType}) != 0)  {
        
      print STATS "\n\n";
      print STATS 
        "STATISTICS FOR ALTERNATIVE $featureType FORMS\n\n";
      print STATS "Total Number of New Forms: ", $totalVar, "\n\n";
      print STATS "% Change from Original Record   Number of Forms         " .
                  "% of Variants\n";

      for (my $i = 0; $i < (scalar @range); $i++)  {
  
        $percent = (${$noChange{$featureType}}[$i] / $totalVar) * 100;
        printf STATS "%13.13s", $range[$i];
        
        if (${$noChange{$featureType}}[$i] =~ /\d/) {
        
          printf STATS "%28.28s", ${$noChange{$featureType}}[$i];
        
        } else {
        
          printf STATS "%28.28s", 0;
        }
        
        printf STATS "%23.2f", $percent;
        print STATS "\n";
      } 
        
      print STATS "\n";
      print STATS "Number of unexpandable variants: ", 
                     $errors{$featureType}, "\n";
        
      # mean, median
        
      my $totalCalc;
        
      for my $calc (@{$allCalc{$featureType}}) {
        
        $totalCalc += $calc;
      }
      
      my $mean = $totalCalc / ($totalVar - $errors{$featureType});
          
      # find median
        
      my @sortedIsoformList = sort _byNumber @{$allCalc{$featureType}};
      my $median;
      
      if ($totalVar%2 == 1) {
        
        $median = $sortedIsoformList[($totalVar -1)/2];
        
      } else {
        
        $median = ($sortedIsoformList[($totalVar -2)/2] + 
                   $sortedIsoformList[($totalVar)/2])/2;
      }  
        
      print STATS "\n";
      print STATS 'Mean % change per variant: ', $mean, "\n"; 
      print STATS 'Median % change per variant: ', $median, "\n";
      close STATS;
    }
  }
}

sub _processOldCommentBlock {

  # for VARIANT and CONFLICT lines

  my ($entry, 
      $featureType, 
      $count,
      $ac,
      $featId,
      $isoformNames, 
      $isoform2feature,
      $featId2start,
      $featId2end,
      $featId2desc,
      $featId2featureType,
      $isoform2note,
      $isHuman) = @_;
  my %diseases;
  my %diseaseVariantCount;
  
  # set up base forms by default
  
  ${${$isoformNames}{$featureType}}[0] = 'Displayed';
  my $isoformNumber = 1; 
   
  # parse feature table
  
  my ($featId, $featId2ftId, $isoform2featId) = _parseFTs($entry, 
                                                         $featureType, 
                                                         $featId,
                                                         $featId2start,
                                                         $featId2end,
                                                         $featId2desc,
                                                         $featId2featureType);
  
  if ($featureType eq 'VARIANT') {
  
    # get disease information from CC lines
  
    my $cc = $entry -> CCs -> filter(&SWISS::CCs::ccTopic("DISEASE"));
    my @draftIsoformNames = keys %$isoform2featId;
    
    if ($isHuman eq "yes") {
      
      # new format used for human disease lines
    
      for my $disease (@{$cc -> list}) {
  
        my $diseaseText = $disease -> comment();
        my ($diseaseName) = ($diseaseText =~ /cause of (.+)[?=,|;]*/i);
        print $diseaseText, "!!\n";
        print $diseaseName, "\t";
      
        # first sentence only
      
        $diseaseName =~ s/\. .*//; 
        $diseaseName =~ s/ \[MIM:.+?\].*//;
      
        if ($diseaseName =~ /\(.+?\)/) {
      
          $diseaseName =~ s/.+?\((.+?)\).*/$1/;
        }
      
        print $diseaseName, "\n";
        $diseases{$diseaseName} = "yes";
      }
    
    } else {
    
      # old format data (non human entires, can't require the same degree of 
      # syntactic precision)
    
      my @diseases;
    
      for my $disease (@{$cc -> list}) {
  
        push @diseases, $disease -> comment();
      }
  
      DISEASEID: for my $draftIsoformName (@draftIsoformNames) {     
        
        my $testName = quotemeta $draftIsoformName;
     
        for my $disease (@diseases) {
    
          if ($disease =~ /$testName/) {
     
            $diseases{$draftIsoformName} = "yes";
            next DISEASEID;
          }
        }
      }
    } 
  
    # now combine feature information and comment information to complete
    # isoforms.  
   
    # we need to check if each isoform name is actually a disease 
    # name... for diseases, treat each 'variant' as separate, i.e. don't
    # combine all sites of variation from the same disease into 
    # a single form
           
    my @draftIsoformNames = keys %$isoform2featId;
  
    ISOFORM: for my $draftIsoformName (@draftIsoformNames) {     
        
      if ($diseases{$draftIsoformName}) {
     
        for my $featId (@{${$isoform2featId}{$draftIsoformName}}) {
      
          # in the case of diseases, we need to create 1 isoform per feature
          # so the draft name isn't good enough, we need additional names
        
          my $isoformName;
        
          if ($$featId2ftId{$featId} =~ /\w/) {
                  
            $isoformName = $draftIsoformName . "-" . $$featId2ftId{$featId};
        
          } else {
           
            $diseaseVariantCount{$draftIsoformName} ++;
            $isoformName = $draftIsoformName . "-CASE " . 
                           $diseaseVariantCount{$draftIsoformName} ++;
          }
        
          push @{${$isoformNames}{$featureType}}, $isoformName; 
          ${${${$isoform2feature}{$featureType}}[$isoformNumber]}{$featId} = 
            "yes";
          $isoformNumber ++;
          $count ++;
        }
        
        next ISOFORM;
      } 
     
      # standard procedure (can have > 1 feature / isoform, as already 
      # specified)
     
      push @{${$isoformNames}{$featureType}}, $draftIsoformName; 
      
      foreach $featId (@{${$isoform2featId}{$draftIsoformName}}) {
      
        ${${${$isoform2feature}{$featureType}}[$isoformNumber]}{$featId} = 
          "yes";
      }
      
      $isoformNumber ++;
      $count ++;
    }
  
  } elsif ($featureType ne 'VARIANT') {
  
    # conflicts are simpler, always merge all conflicts from 1 reference
    # no need to look at CC block
    
    # if feature Ids have not yet been assigned, can use isoform names to 
    # assign data for VARSPLIC keys as well
  
    my @isoformNames = keys %$isoform2featId;
    
    for my $isoformName (@isoformNames) { 
 
      my $revisedName = $isoformName;
      $revisedName =~ s/^isoform //;
      push @{${$isoformNames}{$featureType}}, $revisedName;      
      
      foreach $featId (@{${$isoform2featId}{$isoformName}}) {
      
        ${${${$isoform2feature}{$featureType}}[$isoformNumber]}{$featId} = 
          "yes";
      }
      
      $isoformNumber ++;
      $count ++;
      
      if (($opt_note) && ($featureType eq 'VAR_SEQ')) {
      
        my $cc = 
          $entry -> CCs -> filter(&SWISS::CCs::ccTopic("ALTERNATIVE PRODUCTS"));
     
        if (defined ${$cc -> list}[0]) {
      
          my $block = ${$cc -> list}[0];
          my $keyEvent = $block -> keyEvent();
          
          ${$$isoform2note{'VAR_SEQ'}}{$isoformName} = 
            $block -> getNote($keyEvent, $isoformName); 
        }
      }
    }
  } 
  
  return ($featId, $count);    
}

sub _processNewCommentBlock {

  # for VARSPLIC lines, to be used in conjunction with new style 
  # CC ALTERNATIVE PRODUCTS blocks

  my ($entry, 
      $featureType, 
      $count, 
      $ac, 
      $featId,
      $isoform2id, 
      $isoformNames, 
      $isoform2feature, 
      $featId2start,
      $featId2end,
      $featId2desc,
      $featId2featureType,
      $isoform2note) = @_;
  my $ftId2featId;
  my %isoform2id = %$isoform2id;
  my %isoformNames = %$isoformNames;
  my %isoform2feature = %$isoform2feature;
  my %isoform2note = %$isoform2note;
  my $isoformNumber = 0;
  my ($ft2featId, $ftId2featId, $ft2featId);
  my $cc = 
      $entry -> CCs -> filter(&SWISS::CCs::ccTopic("ALTERNATIVE PRODUCTS"));
     
  if (defined ${$cc -> list}[0]) {
      
    my $block = ${$cc -> list}[0];
    my $keyEvent = $block -> keyEvent();
    
    @{$isoformNames{$featureType}} = 
      $block -> getFormNames($keyEvent);
      
    # but there might be only alternative initiation in this CC block
     
    if (defined ${$isoformNames{$featureType}}[0]) {
        
      # get a hash mapping FT IDs in UniProt record to program's
      # internal feature IDs
       
      ($featId, $ftId2featId, $ft2featId) = _parseFTs($entry, 
                                                      'VAR_SEQ', 
                                                      $featId,
                                                      $featId2start,
                                                      $featId2end,
                                                      $featId2desc,
                                                      $featId2featureType);
      
      my %ftId2featId = %$ftId2featId;
      
      for my $isoformName (@{$isoformNames{$featureType}}) {
      
        if ($opt_note) {
        
          ${$isoform2note{'VAR_SEQ'}}{$isoformName} = 
            $block -> getNote($keyEvent, $isoformName);
        }
        
        # now map from each internal feature ID back to the forms it is
        # contained in 
             
        # and set up a data structure allowing one to map from a given 
        # isoform specification back to its constituant features
             
        my @ftIds = 
          $block -> getFeatIds($keyEvent, $isoformName);
          
        if (($ftIds[0] eq "") || 
            ($ftIds[0] eq "Not described") || 
            ($ftIds[0] eq "External")) {
            
          
          # isoform not suitable for VARSPLIC expansion in this entry
          # just ignore
          
          ${${$isoform2feature{'VAR_SEQ'}}[$isoformNumber]}{"Ignore"} = "yes";
          push @{$isoform2id{'VAR_SEQ'}}, "-";
          $isoformNumber ++;
          
        } else {
        
          if ($ftIds[0] ne 'Displayed') {
            
            # no need of VARSPLIC expansion for isoforms with 'Displayed'
            # specification
            # first specification of every entry should always be 
            # 'Displayed'
              
            for my $ftId (@ftIds) {
        
              # says that the nth VARSPLIC specification contains this
              # feature
              
              my $thisFeatId = $ftId2featId{$ftId};
              ${${$isoform2feature{'VAR_SEQ'}}[$isoformNumber]}{$thisFeatId} =
               "yes";
            }
            
            # use count variable to count number of alternatives
                
            $count ++;
          
          } else {
          
            ${${$isoform2feature{'VAR_SEQ'}}[$isoformNumber]}{"Master"} 
              = "yes";
          }
           
          $isoformNumber ++;
          my @isoIds = 
            $block -> getIsoIds($keyEvent, $isoformName);
            
          push @{$isoform2id{'VAR_SEQ'}}, $isoIds[0];
        } 
      }
    }
  } 
   
  if ($isoformNumber == 0) {
      
    # no ALTERNATIVE PRODUCTS blocks, can assume no VARSPLIC features
    # either
        
    ${$isoformNames{'VAR_SEQ'}}[0] = 'Displayed';
    ${$isoform2id{'VAR_SEQ'}}[0] = $ac;
    ${${$isoform2feature{'VAR_SEQ'}}[0]}{"Master"} = "yes";
  }
  
  if (($opt_check_vsps) && ($count > 0)) {
  
    # double check of data derived from CC lines with data derived from FT
    # lines
  
    my $isoformNumber = 0;
    
    for my $isoformName (@{$isoformNames{$featureType}}) {
    
      my %ok;
      
      # isform2featId hash is not defined for parental form
      
      my $tempIsoformName = 'isoform ' . $isoformName;
      
      if (defined ${$ft2featId}{$tempIsoformName}) {
      
        my @ftFeatureIds = @{${$ft2featId}{$tempIsoformName}};
    
        for my $ftFeatureId (@ftFeatureIds) {
    
          if 
      (! ${${$isoform2feature{'VAR_SEQ'}}[$isoformNumber]}{$ftFeatureId}) {
      
            print ERROR "Inconsistent VSPs for $tempIsoformName of entry " . 
                        $ac . "\n";
        
          } else {
        
            $ok{$ftFeatureId} = "yes";
          }
        }
      }
      
      my @ccFeatureIds = 
        keys %{${$isoform2feature{'VAR_SEQ'}}[$isoformNumber]};
      
      for my $ccFeatureId (@ccFeatureIds) {
      
        if (($ccFeatureId ne "Master") && 
            ($ccFeatureId ne "Ignore") &&
            (! $ok{$ccFeatureId})) {
        
          my @vsps = keys %$ftId2featId;
        
          for my $vsp (@vsps) {
          
            if (${$ftId2featId}{$vsp} eq $ccFeatureId) {
        
              print ERROR $vsp . " absent from >= 1 isoform description " .
                          "in entry " . $ac . "\n";
            }
          }
        }
      }
      
      $isoformNumber ++;
    }   
  }
  
  return ($featId, $isoformNumber);
}

sub _seqCut {

  ## cuts (and adds) bases to a sequence 

  # name the parameters
  
  my($sq, $start, $stop, $add, $sizechange) = @_;
 
  # declare other variables
  
  my ($pre, $post, $cut);
  my ($st, $sp);
  my ($addsize, $cutsize);
  my $unmatched;
  
  # Use supplied parameters to calculate certain postions in string
  
  $st = $start + $sizechange - 1;
  $sp = $stop + $sizechange - $st;
  
  # can't handle titin correctly, need to prevent program from crashing
  
  if (($st > 32766) || ($sp > 32766)) {
  
    return;
  }
  
  # Cut $sq from $start to $stop: make insertion where appropriate
  
  ($pre, $cut, $post) = ($sq =~ /(\w{$st})(\w{$sp})(\w*)/s);
  $sq = $pre . $add . $post;  
  
  #  Calculate change in size of seqeunce
  
  $addsize = ($add =~ tr/A-Z//);
  $cutsize = ($cut =~ tr/A-Z//);
  
  $sizechange += $addsize - $cutsize;
 
  # $unmatched holds the number of unaligned bases between the orginal protein 
  # and its isoform in gapped alignment (see README for details)
  
  $unmatched = $addsize + $cutsize;
  
  # adjustment for matching first resiude s(annotation artefact)
  
  my ($firstOldAA) = ($add =~ /^(\w)/);
  my ($firstNewAA) = ($cut =~ /^(\w)/);
  
  if ($firstOldAA eq $firstNewAA) {
  
    $unmatched -= 2;
    $addsize -=2;
  }
    
  # return modifed sequence, cut sequence, net sizechange (for co-ordiante
  # calculation), no. of unmatched  bases and number of added bases to gapped 
  # alignment
  
  return ($sq, $cut, $sizechange, $unmatched, $addsize);
}

sub _tidy {

  # used to tidy up variant and conflict entries
  
  my ($var) = @_;
  
  $var =~ s/\bIN //;
  $var =~ s/\bAND\W//;
  $var =~ s/^ //;
  return $var;
}

__END__

 
