#!/usr/bin/perl -w

## Add extra names to fasta file from synonyms.pl output


my %synonyms;
my $synonyms_file=$ARGV[0] || die("need a synonyms file");
my $fasta_file=$ARGV[1] || die("need a fasta file");
#Q9BWF2
open(S, "$synonyms_file");
while(<S>){
    s/://;
    s/\s+/ /g;
    s/^\s//;
    my @a=/(.+?)\s/g;
    my %k;
    map{$k{$_}++}@a;
    @a=keys(%k);
    map{$synonyms{$_}=\@a}@a;
}
open(F, "$fasta_file");
while(<F>){
    if(/>/){
	chomp;
	my $idline=$_;
	my @a=split(/[\W\.]/);
	name:foreach my $name (@a){
	    next if length($name)<1;
	    if (defined($synonyms{$name})){
		chomp;
		$idline .= " @{$synonyms{$name}}";
		next name;
	    }
	   
	}
	print "$idline\n";
    }
    else{
	print;
    }

}
