#!/usr/bin/perl -w
use strict;
use Getopt::Std;
use Switch;
my %opts;
getopts('s:f:',\%opts);
my $synfile=$opts{s};
my $flat=$opts{f};
my $go_terms_file='/home/terdon/research/testing/new/data/GO.terms_alt_ids';
my %synonyms;
$synonyms{LOADED}=0;
my @names;
my $from;
if($synfile =~/hs.uni2acc.map/) {$from='AC'}
elsif($synfile =~/dro.uni2fb.map/) {$from='FlyBase'}
elsif($synfile =~/mus.uni2mgi.map/) {$from='MGI'}
elsif($synfile =~/scc.uni2sgd.map/) {$from='SGD'}
elsif($synfile =~/ele.uni2wb.map/) {$from='WormBase'}
else{die("Bad synfile : $synfile\n");}


open(A,"$ARGV[0]");
while(<A>){
    chomp;
    push @names,&get_name($_);

}

#print STDERR "uniprot_parse.pl -Ff ID -n \"@names\" $flat\n";
system("uniprot_parse.pl -vFf $from -n \"@names\" $flat");


sub get_name{
    my $name=$_[0];
    my $old_name=$name;
    my $called_by=$_[1]||"null";
    my $reverse=$_[2]||undef;
    if ($synonyms{LOADED}==0){
	if(-e $synfile){
	    open(S,$synfile);
	    while(<S>){
		if(/^(.+)\t(.+)/){
		    $synonyms{NAME}{$1}=$2;
		}	
                
	    }
	    close(S);
	}
	else{die("Could not open $synfile\n" );
	 }	
 	$synonyms{LOADED}=1;
    }
    if(defined($synonyms{NAME}{$name})){
	$name=$synonyms{NAME}{$name} ;
    }
    
    return $name;
    
}
