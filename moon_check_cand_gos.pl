#!/usr/bin/perl -w

use strict;
use Switch;

my (@candidates,%gos,%hits,%terms_to_GOs,%seen);
my $have_already_read_terms_file=0;

## Parse tehe output of moon_recip_best.pl
my $parsed_file=shift;
open(A,"$parsed_file");
my $c=0;
while(<A>){
    chomp;
    my @a=split(/\t/);
    map{s/.+_(.+_.+)/$1/}@a;
    $candidates[$c]=\@a;
    map{$hits{$_}++; }@a;
    $c++;
}
close(A);

foreach my $file (@ARGV){
    $file=~/(.+)\.([bimpPcd])\./;
    my $mode=$2;
    my $species=$1;
    
    open(F,"$file");
    while(<F>){
	
	if ($mode eq 'b') {
	    /^(.+?)\s+.+?(GO:\d+)/||die("$mode:$_\n");
	    my ($a,$b)=($1,$2);
	    next unless defined($hits{$a});
	    push @{$gos{$a}{$species}},$b unless $seen{$a}{$species}{$b};
	    $seen{$a}{$species}{$b}++;
	    #print "\@{\$gos{$a}{$species}},$b;\n"
	}
	elsif($mode=~/[cd]/){
	    /^(.+?)\s.+?(GO:\d+).+?\s+:\s+([\w_]+).+?(GO:\d+)/||die("$mode:$_\n");
	    my ($a,$b,$c,$d)=($1,$2,$3,$4);
	    #print "popo  $a : $hits{$a}, $c $hits{$c}\n" if $a =~/LAMB/;
	    unless (defined($hits{$a}) ||  defined($hits{$c})){next}
	    #print "LAMB2_HUMAN- xx -$a-$b-\n$_";die();	    
	    push @{$gos{$a}{$species}},$b unless $seen{$a}{$species}{$b};
	    push @{$gos{$c}{$species}},$d unless $seen{$c}{$species}{$d};
	    #print "\@{\$gos{$a}{$species}},$b;\n"
	    $seen{$a}{$species}{$b}++;
	    $seen{$c}{$species}{$d}++;
	}
	elsif($mode eq 'i'){
	    /^(.+?)\s+.+?\)\s+(.+?)\t/||die("$mode:$_\n");
	    my ($a,$b)=($1,$2);
	    next unless defined($hits{$a});
	    my $t=&terms_to_GOs($b,0);
	    push @{$gos{$a}{$species}},$t unless $seen{$a}{$species}{$t};
	    #print "\@{\$gos{$a}{$species}},$b;\n"
	    $seen{$a}{$species}{$t}++;

	}
	elsif($mode eq 'm'){
	    /^(.+?)\s+\(.+?(GO:\d+)\,/||die("$mode:$_\n");
	    my ($a,$b)=($1,$2);
	    next unless defined($hits{$a});
	    push @{$gos{$a}{$species}},$b unless $seen{$a}{$species}{$b};
	    #print "\@{\$gos{$a}{$species}},$b;\n"
	    $seen{$a}{$species}{$b}++;
	}
	else{die("Unknown case : $file : -$mode-\n");}
	
	
    }
    close(F);
}

for (my $n=0; $n<=$#candidates; $n++){
#foreach my $array (@candidates){
    print "@{$candidates[$n]}\t";
    foreach my $hit (@{$candidates[$n]}){
	my @keys=keys(%{$gos{$hit}});
#	print "KK : @keys\n";
	foreach my $species (keys(%{$gos{$hit}})){
	    print "$species:";
	    map{print "$_ " . &terms_to_GOs($_,1) . "::"}@{$gos{$hit}{$species}};
	}
	
    }
    print "\n";
}


############################################################
sub terms_to_GOs{
#    my $term='transcription, DNA-dependent';#shift;
 #   my $mode=0;#shift; ## 0 will return GO:xxx, 1 will return term name
    my $term=shift;
    my $mode=shift; ## 0 will return GO:xxx, 1 will return term name

#    unless($term eq 'biological_process'){die("crapiola : -$term-\n");}
    $term=~s/_/ /g;
    if($have_already_read_terms_file==0){
	open(T,"/home/terdon/research/GO/GO.terms_alt_ids.txt")|| die("Cannot open terms file : $!\n");
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
	    #print "xx $a[$#a-1]- : $a[0]\nxx $term-\n" if $line=~/GO:0006351/;

	    if($line=~/obs$/){pop(@a);}
	    if($aa[1] =~/GO:/){
		push @terms,split(/\s+/,$aa[1]);
	    }
	    else{
		$terms_to_GOs{TERMS}{$a[$#a-1]}=$a[0];
		map{$terms_to_GOs{GOs}{$_}=$a[$#a-1];}@terms;
	    }
	}
	close(T);
	$have_already_read_terms_file=1;
    }
    print "$term : $terms_to_GOs{TERMS}{$term}  : $terms_to_GOs{GOs}{$term}\n";
#    &debug("term : $term, id:$terms_to_GOs{TERMS}{$term}, id:$terms_to_GOs{TERMS}{$term} " );
#    print STDERR "term : $term, id:$terms_to_GOs{TERMS}{$term} \n";
    #unless( defined($terms_to_GOs{GOs}{$term})|| defined($terms_to_GOs{GOs}{$term})){die("$term\n")}
    $mode==0 ? 
	return($terms_to_GOs{TERMS}{$term}) :
	return($terms_to_GOs{GOs}{$term}) ;
    
}

#GO:0051276
