#!/usr/bin/perl -w

use strict;
use Getopt::Std;

my (%opts,%synonyms,%interactors,%terms_to_GOs,%want);
getopts('hT:N:o:s:O:D:',\%opts);
$synonyms{LOADED}=0;
my $DATADIR=$opts{D}||"/home/terdon/research/testing/new/data";
my $synfile=$opts{s}||"$DATADIR/hs.uni2acc.map";
my $net_file=$opts{N}||"$DATADIR/../HSHQ.gr";
my $cand=$ARGV[0]||die "USAGE:\n\t$0 <cand_prot>";
my $gaf_file=$opts{O}||"$DATADIR/gene_association.human";
my $subonto=$opts{o}||'P';
my $go_terms_file=$opts{T}||"$DATADIR/GO.terms_alt_ids";
my $have_already_read_terms_file=0;
open(N,"$net_file")||die "Could not open $net_file :$!\n";
while(<N>){
    chomp;
    next unless /$cand/;
    my @a=split(/\t/);
    my $n;
    $a[0] eq "$cand" ? ($n=&get_name($a[1],1)) : ($n=&get_name($a[0],1));
    $want{$n}++ ;
}
close(N);
open(G,"$gaf_file")||die "Could not open $gaf_file :$!\n";
while(<G>){
    next if /^!/;
    chomp;
    my @tmparray=split(/\t/,$_);
    next unless $tmparray[8] eq $subonto;
    next if $tmparray[3] !~/^$/;
    my $name=&get_name($tmparray[1],1);
    #print "bnn : $name : $tmparray[1]\n";
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
foreach my $name (keys(%interactors)){
    my @terms=keys(%{$interactors{$name}{TERMS}});
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
