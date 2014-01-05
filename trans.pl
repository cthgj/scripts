#!/usr/bin/perl -w

use strict;
use Getopt::Std;

my %opts;
getopts('f:',\%opts)|| die("no opts\n");
my $a = 1;
my $seq = "";
my @dna;
my $name;
my $frame = $opts{f} || 0;
if($frame < 0 || $frame >2){die("frame must be between 0 and 2\n") }
my %aa = ();
while(<>)
{
    
    $a = (( />/) ?  1 : 0);
   if($a == 0)
   {
       chomp;
       $seq = $seq . $_;
#       @dna = ($seq =~ /(...)/g);
              
   }
    else
    {
	&translate($seq) unless $. == 1;;
	$name = $_;
	chomp $name;
	$seq = "";

    }


}
&translate($seq);


sub translate
{
   my @prot;

    $aa{GCA} = "A";    $aa{GCC} = "A";    $aa{GCG} = "A";    $aa{GCT} = "A";    $aa{TGC} = "C";    $aa{TGT} = "C";
    $aa{GAC} = "D";    $aa{GAT} = "D";    $aa{GAA} = "E";    $aa{GAG} = "E";    $aa{TTC} = "F";    $aa{TTT} = "F"; 
    $aa{GGA} = "G";    $aa{GGC} = "G";    $aa{GGG} = "G";    $aa{GGT} = "G";    $aa{CAC} = "H";    $aa{CAT} = "H";  
    $aa{ATA} = "I";    $aa{ATC} = "I";    $aa{ATT} = "I";    $aa{AAA} = "K";    $aa{AAG} = "K";    $aa{TTA} = "L";   
    $aa{TTG} = "L";    $aa{CTA} = "L";    $aa{CTC} = "L";    $aa{CTG} = "L";    $aa{CTT} = "L";    $aa{ATG} = "M";   
    $aa{AAC} = "N";    $aa{AAT} = "N";    $aa{CCA} = "P";    $aa{CCC} = "P";    $aa{CCG} = "P";    $aa{CCT} = "P";  
    $aa{CAA} = "Q";    $aa{CAG} = "Q";    $aa{AGA} = "R";    $aa{AGG} = "R";    $aa{CGA} = "R";    $aa{CGC} = "R";  
    $aa{CGG} = "R";    $aa{CGT} = "R";    $aa{AGC} = "S";    $aa{AGT} = "S";    $aa{TCA} = "S";    $aa{TCC} = "S";  
    $aa{TCG} = "S";    $aa{TCT} = "S";    $aa{ACA} = "T";    $aa{ACC} = "T";    $aa{ACG} = "T";    $aa{ACT} = "T";  
    $aa{GTA} = "V";    $aa{GTC} = "V";    $aa{GTG} = "V";    $aa{GTT} = "V";    $aa{TGG} = "W";    $aa{TAC} = "Y";   
    $aa{TAT} = "Y";    $aa{TAA} = "!";    $aa{TAG} = "#";    $aa{TGA} = "@";    

 
   $aa{gca} = "A";    $aa{gcc} = "A";    $aa{gcg} = "A";    $aa{gct} = "A";    $aa{tgc} = "C";    $aa{tgt} = "C";   
   $aa{gac} = "D";    $aa{gat} = "D";    $aa{gaa} = "E";    $aa{gag} = "E";    $aa{ttc} = "F";    $aa{ttt} = "F";    
   $aa{gga} = "G";    $aa{ggc} = "G";    $aa{ggg} = "G";    $aa{ggt} = "G";    $aa{cac} = "H";    $aa{cat} = "H";     
   $aa{ata} = "I";    $aa{atc} = "I";    $aa{att} = "I";    $aa{aaa} = "K";    $aa{aag} = "K";    $aa{tta} = "L";  
   $aa{ttg} = "L";    $aa{cta} = "L";    $aa{ctc} = "L";    $aa{ctg} = "L";    $aa{ctt} = "L";    $aa{atg} = "M"; 
   $aa{aac} = "N";    $aa{aat} = "N";    $aa{cca} = "P";    $aa{ccc} = "P";    $aa{ccg} = "P";    $aa{cct} = "P"; 
   $aa{caa} = "Q";    $aa{cag} = "Q";    $aa{aga} = "R";    $aa{agg} = "R";    $aa{cga} = "R";    $aa{cgc} = "R";    
   $aa{cgg} = "R";    $aa{cgt} = "R";    $aa{agc} = "S";    $aa{agt} = "S";    $aa{tca} = "S";    $aa{tcc} = "S";    
   $aa{tcg} = "S";    $aa{tct} = "S";    $aa{aca} = "T";    $aa{acc} = "T";    $aa{acg} = "T";    $aa{act} = "T";    
   $aa{gta} = "V";    $aa{gtc} = "V";    $aa{gtg} = "V";    $aa{gtt} = "V";    $aa{tgg} = "W";    $aa{tac} = "Y";    
   $aa{tat} = "Y";    $aa{taa} = "!";    $aa{tag} = "#";    $aa{tga} = "@";



   $seq =~ s/^.{$frame}//;
   @dna= ($seq =~ /(...)/g);
   for (my $i=0; $i<scalar(@dna); $i++)
   {
       
       $aa{$dna[$i]} = "X" unless defined($aa{$dna[$i]});
       $prot[$i] = $aa{$dna[$i]} || die("codon :$dna[$i]\n") ;
   }
   print "$name\n";
   for(my $i=0; $i<scalar(@prot); $i++)
   {
       print "$prot[$i]";
       unless($i==0){print "\n" if $i % 60 == 0;}
       
   }
   
   print "\n";
   
}
