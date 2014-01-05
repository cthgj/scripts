#!/usr/bin/perl 
use strict;
use Getopt::Std;
## This script will look for cases where P(ancestor{go1}-go2)<P(go1-go2)



my %opts;
getopts('hvo:g:',\%opts) ;
my $geneology_file=$opts{g}||"/home/cchapple/research/GO/biological_process.geneology";
my %prob;
my $verbose=$opts{v}||undef;


my %ancestors=&load_ancestors($geneology_file);
open(A,$ARGV[0]);
while(<A>){   #37483939
    next if /GO:0008150/;
    printf STDERR ("$. of 37483939\r");
    chomp;
    s/_/\t/;
    my @a=split(/\t/);
#    $prob{$a[0]}{$a[1]}=$prob{$a[1]}{$a[0]}=$a[$#a];
# 
    if(defined($ancestors{$a[0]}{$a[1]})){
	print STDOUT "$a[0]_$a[1]\t$a[$#a]\n";
    }
    elsif(defined($ancestors{$a[1]}{$a[0]})){
	print STDOUT "$a[1]_$a[0]\t$a[$#a]\n";

    }
    else{next}
    
}

print STDERR "\n--------------------------\n";



# ## Check if ancestor:child is among the underrepresented
# foreach my $g1 (keys(%prob)){
#     foreach my $papa (@{$ancestors{$g1}}){
# 	if (defined($prob{$papa}{$g1})){
# 	    print "$g1_$papa\t$prob{$papa}{$g1} \n";
# 	}
# 	elsif(defined($prob{$g1}{$papa})){
# 	    print "$g1_$papa\t$prob{$g1}{$papa} \n";
# 	}
# 	else{next};


#     }
# }



sub load_ancestors{
    print STDERR "Loading ancestors...\n" if $verbose;
    open (ANC,"$_[0]")|| die("cannot open $_[0]:$!\n");
    my (%ancestors,%offspring);
    while(<ANC>){
	chomp;
	my @terms=split(/\t/);
	my %seen;
	@terms = grep { ! $seen{$_} ++ } @terms;
	my $child=shift(@terms);
	map{$ancestors{$child}{$_}++}@terms;
    }
    close(ANC);
    return(%ancestors)
}
