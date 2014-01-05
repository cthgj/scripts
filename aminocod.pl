#!/usr/bin/perl -w

use strict;

#
# declaracio de variables
#

my %codamin=("TTT" => "Phe",
          "TTC" => "Phe",
          "TTA" => "Leu",
          "TTG" => "Leu",
          "CTT" => "Leu",
          "CTC" => "Leu",
          "CTA" => "Leu",
          "CTG" => "Leu",
          "ATT" => "Ile",
          "ATC" => "Ile",
          "ATA" => "Ile",
          "ATG" => "Met",
          "GTT" => "Val",
          "GTC" => "Val",
          "GTA" => "Val",
          "GTG" => "Val",
          "TCT" => "Ser",
          "TCC" => "Ser",
          "TCA" => "Ser",
          "TCG" => "Ser",
          "CCT" => "Pro",
          "CCC" => "Pro",
          "CCA" => "Pro",
          "CCG" => "Pro",
          "ACT" => "Thr",
          "ACC" => "Thr",
          "ACA" => "Thr",
          "ACG" => "Thr",
          "GCT" => "Ala",
          "GCC" => "Ala",
          "GCA" => "Ala",
          "GCG" => "Ala",
          "TAT" => "Tyr",
          "TAC" => "Tyr",
          "TAA" => "Stop",
          "TAG" => "Stop",
          "CAT" => "His",
          "CAC" => "His",
          "CAA" => "Gln",
          "CAG" => "Gln",
          "AAT" => "Asn",
          "AAC" => "Asn",
          "AAA" => "Lys",
          "AAG" => "Lys",
          "GAT" => "Asp",
          "GAC" => "Asp",
          "GAA" => "Glu",
          "GAG" => "Glu",
          "TGT" => "Cys",
          "TGC" => "Cys",
          "TGA" => "Stop",
          "TGG" => "Trp",
          "CGT" => "Arg",
          "CGC" => "Arg",
          "CGA" => "Arg",
          "CGG" => "Arg",
          "AGT" => "Ser",
          "AGC" => "Ser",
          "AGA" => "Arg",
          "AGG" => "Arg",
          "GGT" => "Gly",
          "GGC" => "Gly",
          "GGA" => "Gly",
          "GGG" => "Gly");

my %abrev=("A" => "Ala",
           "R" => "Arg",
           "N" => "Asn",
           "D" => "Asp",
           "C" => "Cys",
           "Q" => "Gln",
           "E" => "Glu",
           "G" => "Gly",
           "H" => "His",
           "I" => "Ile",
           "L" => "Leu",
           "K" => "Lys",
           "M" => "Met",
           "F" => "Phe",
           "P" => "Pro",
           "S" => "Ser",
           "X" => "Ter",
           "T" => "Thr",
           "W" => "Trp",
           "Y" => "Tyr",
           "V" => "Val",
           "*" => "Stop",);

my $entrada;                              # variable que enregistra la sequencia
my @codons;                               # vector que enregistra els codons de la sequencia
my $i;                                    # variable per anar recorrent el vector anterior
my %revabr=reverse(%abrev);               # hash per enregistrar l'invers del hash ``abrev''

$entrada = <STDIN>;                       # llegim la sequencia d'entrada pel teclat
chomp($entrada);                          # eliminem el canvi de linia (si existeix)

@codons = ($entrada =~ m/.../g);          # enregistrem els codons utilitzant una exp. reg.


$i=0;
while ($i < scalar(@codons)) {

  print "$revabr{$codamin{$codons[$i]}}"; # imprimim l'abreviatura utilitzant els hash
  $i = $i + 1;                            # de codon-aminoacid primer, i
                                          # aminoacid-abreviatura despres
}

print "\n";                               # canvi de linia per a que no ens surti el
                                          # simbol del shell enganxat a la sortida del
                                          # programa

