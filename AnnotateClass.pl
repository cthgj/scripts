#!/usr/bin/perl
use strict;



#Majority rule, plus de 30% par défaut, paramétrable
#ouverture du fichier avec la liste des genes et leur association (une ligne par gene)
open (F, $ARGV[0]) or die "cannot open $ARGV[0]\n";
#ouverture du fichier classes
open (G, $ARGV[1]) or die "cannot open $ARGV[1]\n";

my %asso=();
#Lecture du fichier des genes et leurs associations
while (my $line=<F>){
    if ($line=~/^(.*?)\t(.*)$/){
	chomp $line;
	my $geneid=$1;
	my $endline=$2;
	#print "$geneid\n";
	#print "$endline\n";
	$asso{$geneid}=$endline;
    }
    
}
#while ((my $keys, my $values) = each (%asso)) {
#			 	print "key $keys\t value $values\n";
#		 
#		 }

my %countannot=();
#Lecture du fichier des classes
while (my $ligne=<G>){
    if ($ligne=~/.{4,}\t/) {
	my $i =0;
	my %countannot=();
	chomp $ligne;
	my @genelist=split("\t", $ligne);
	#print "Lecture Classe\n";
	#print "@genelist\n";
	foreach my $gene (@genelist){
	    #print "Lecture Gene\n";
	    $i++;
	    #print "$gene\n";
	    if (exists ($asso{$gene})) {
		#print "gene $gene\n";
		#print "asso de gene $asso{$gene}\n";
		my @annot=split("\t", $asso{$gene});
		#print "@annot\n";
		foreach my $annot (@annot){
		    #print "annot $annot\n";
		    if (exists ($countannot{$annot})) {
			$countannot{$annot}++;
		    } else {
			$countannot{$annot}=1;
		    }
		    
		}
	    }
	    
	    
	}
	#print "\n";	
	#print "Etat table comptage\n";
	print "#####taille classe $i#####\n";
	while ((my $keys, my $values) = each (%countannot)) {
	    my $test=$values/$i;
	    
	    if ($test>0.3){
		#print "$keys\t$values\t";
		print "$keys\t$values over $i ==> $test\n";
		#print "$keys\t$values\n";
	    }
	}
    }	
}
