#!/usr/bin/perl -w

use strict;
use Getopt::Std;
my $go_terms_file="$ENV{HOME}/research/testing/new/data/GO.terms_alt_ids";
my $have_already_read_terms_file=0;
my (%terms_to_GOs,%opts,%prob,%pairs);
my $subonto='P';
getopts('p:D:F:l',\%opts);
my $load_entire_file=$opts{F}||undef;
my $stats_dir=$opts{D}||undef;
my $prob_file_read=0;
my $min_prob=$opts{p}||0.09157528;
my $print_line=$opts{l}||undef;


foreach my $file (@ARGV){
    $file=~/(.+)\.([baimpPcd])\./;
    my $mode=$2;
    my $species=$1;
    
    open(F,"$file");
    while(<F>){
	my $name;
	my $pair;
	my $line=$_;
	chomp;
	if ($mode eq 'b') {
	    (my @m) = (/(GO:\d+)/g);
	    $pair=join("_", sort {$b lt $a} @m);
	    /^(.+?)\s+/||die("$mode:$_\n");
	    $name=$1;
	}
	elsif($mode eq'a'){
	    (my @m) = (/(GO:\d+)/g);
	    $pair=join("_", sort {$b lt $a} @m);
	    /^(.+?)\s+/||die("$mode:$_\n");
	    $name=$1;
	}
	elsif($mode=~/[cd]/){
	    (my @m) = (/(GO:\d+)/g);
	    $pair=join("_", sort {$b lt $a} @m);
	    /^(.+?)\s+/||die("$mode:$_\n");
	    $name=$1;
	}
	elsif($mode eq 'i'){
	    /^(.+?)\s+.+?\)\s+(.+?)\t(.+)\s*\t/||die("$mode:$_\n");
	    my $i=&terms_to_GOs($2,0);
	    my $ii=&terms_to_GOs($3,0);
	    $pair=join("_", sort {$b lt $a} ($i,$ii));
	    $line=~/^(.+?)\s/||die("$mode:$_\n");
	    $name=$1;
	}
	elsif($mode eq 'm'){
	    (my @m) = (/(GO:\d+)/g);
	    $pair=join("_", sort {$b lt $a} @m);
	    $line=~/^(.+?)\s.+:\s+(.+?)\s/||die("$mode:$_\n");
	    $name=$1 . "\n" . $2;
	    
	}
	else{die("Unknown case : $file : -$mode-\n");}
	$pairs{$pair}++;
	my $P=&interactome_probability($pair);
	if($print_line){
	    print "$line" if $P <= $min_prob;
	}
	else{
	    print "$name\n" if $P <= $min_prob;
	    }
    }
    close(F);
}



############################################################
sub interactome_probability{
    my $gopair=shift;
    unless($prob_file_read==1){
	if($load_entire_file){
	    open(A,"$load_entire_file");
	    print STDERR "Reading";
	    while(<A>){
		print STDERR "." if $. % 100000 ==0;
		chomp;
		my @a=split(/\t/);
		$prob{$a[0]}=$a[6];
	    }
	    $prob_file_read=1;
	}
	else{die("Haven't done this yet\n")}
	print STDERR "Done\n";

    }
    defined($prob{$gopair}) ? 
	return($prob{$gopair}):
	return(1);
}
############################################################
sub terms_to_GOs{
    my $term=shift;
    $term=~s/^\s*//;
    $term=~s/\s*$//;
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
