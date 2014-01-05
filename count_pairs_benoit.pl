#!/usr/bin/perl -w

my $geneology_file=$ARGV[0];
my $gaf_file=$ARGV[1];


my %ancestors=&read_geneology($geneology_file);

my(%gos,%go_count);
open(A,"$gaf_file")||die("cannot open $gaf_file: $!\n");
while(<A>){
    next if /^!/;
    chomp;
    my @tmparray=split(/\t/);
    $gos{$tmparray[1]}{$tmparray[4]}++; ## $tmparray[1]== OMIM, $tmparray[4] == HP
    $go_count{$tmparray[4]}{$tmparray[1]}++;
    #print STDERR "1: \$go_count{$tmparray[4]}{$tmparray[1]}\n";
    # map{$gos{$tmparray[5]}{$_}++; $go_count{$_}{$tmparray[5]}++;}@{$ancestors{ANC}{$tmparray[4]}

}
close(A);
## Number of diseases with >1 DIRECT phenotype
my $total_gt1;
map{$total_gt1++ if (keys(%{$gos{$_}})) > 1}keys(%gos);

foreach my $omim (keys(%gos)){
    foreach my $hp(keys(%{$gos{$omim}})){
	map{
	    ## Inherit ancestor terms
	    $gos{$omim}{$_}++;
            ##count them
	    $go_count{$_}{$omim}++;
	}@{$ancestors{ANC}{$hp}};
    }
#    map{print "$_\t";}keys(%{$gos{$omim}});
}
my @GOs=keys(%go_count);
foreach my $hp(@GOs){
    map{$go_count{tot}{$hp}++;}keys(%{$go_count{$hp}});
}
for (my $n=0; $n<=$#GOs;$n++){
    for (my $k=$n+1; $k<scalar(@GOs);$k++){
	print STDERR "$n:$k of $#GOs\r";
	my @mal1=keys(%{$go_count{$GOs[$n]}});
#	my @mal2=keys(%{$go_count{$GOs[$k]}});
	my $both=0;
	map{
	    if (defined($gos{$_}{$GOs[$k]})){
		$both++;
	
	    }
	}@mal1;

	print "$GOs[$n]_$GOs[$k]\t$total_gt1\t$go_count{tot}{$GOs[$n]}\t$go_count{tot}{$GOs[$k]}\t$both\n";
    }
}

sub read_geneology{
## Read geneology file
    my $file=shift;
    open(A,"$file")||die("cannot open $file: $!\n");
    while(<A>){
	chomp;
	s/\|/\t/g;
	my @terms=split(/\t/);
	my $child=shift(@terms);
	@{$ancestors{ANC}{$child}}=@terms;
	map{push @{$ancestors{KIDs}{$_}},$child}@terms;
    }
    close(A);
    return(%ancestors)
}
