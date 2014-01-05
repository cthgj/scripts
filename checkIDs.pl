#!/usr/bin/perl 

### The idea is that this script will take list of fish-tetra hits and print out the known SP
### each hit againts. Doesn't work too well...



use strict;

my $tetras = '/home/ug/cchapple/research/selenoproteins/tetraodon/SECISearch/sergi_cDNAs/cDNA.known.tetraodon.SP_table';
my $fishs = '/home/ug/cchapple/research/selenoproteins/tetraodon/SECISearch/sergi_cDNAs/cDNA.known.fish.SP_table';
my %fishsel = ();
my %tetrasel = ();
my $ids = '/home/ug/cchapple/research/selenoproteins/tetraodon/SECISearch/fishes/blast/fish_tetra.std.out.pairs';

my %tetrapair = ();
my %fishpair = ();


open(ID,"$ids")|| die "shit :$!\n";
while(<ID>) 
{
    /^(.*?)\s(.*?)$/;

    $fishpair{$2} = $1;  # {fishname} = tetra
    $tetrapair{$1} = $2; # $tetrapair{tetraneme} = fishname
        
}
my @tetra = keys(%tetrapair);
my @fishes = keys(%fishpair);

close(ID);

open(FISH,"$fishs");
while(<FISH>){
    foreach my $fi (@fishes)    {
	if (/$fi/)	{
	    /^(.*?)\s/; 
	    $fishsel{$fi} = $1;
	}
    }
}
open(TETRA,"$tetras");
while(<TETRA>){
    foreach my $te(@tetra)    { 
	if (/$te/)	{
	    /^(.*?)\s/; 
	    $tetrasel{$te} = $1;
	}
    }
}
close(TETRA);
foreach my $te (@tetra){
    chomp $te;
    if ($tetrasel{$te}){
	print STDOUT "$te\t$tetrapair{$te} \t$fishsel{$tetrapair{$te}}(fish) \t$tetrasel{$te} \n" if $tetrasel{$te} eq $fishsel{$tetrapair{$te}}

    }
}

print STDOUT "********************\n";
foreach my $te (@tetra){
    chomp $te;
    if ($tetrasel{$te}){
	print STDOUT "$te\t$tetrapair{$te} \t$fishsel{$tetrapair{$te}}(fish) \t$tetrasel{$te} \n" unless $tetrasel{$te} eq $fishsel{$tetrapair{$te}}

    }
}
print STDOUT "********************\n";
foreach my $te (@tetra){
    chomp $te;
    unless ($tetrasel{$te}){
	print STDOUT "$te\n";#\t$tetrapair{$te} \t$fishsel{$tetrapair{$te}}(fish) \t$tetrasel{$te} \n" unless $tetrasel{$te} eq $fishsel{$tetrapair{$te}}

    }
}


my  @nontet = qw(FD0AHC72AC02.contig FD0ADA11BH09.contig FD0AHC14BF04.contig FD0AEB2AB05.contig FD0AHA15DH09.contig FD0ADA11CH04.contig FD0ADB6AG03.contig FD0AHA49AA12BBP1 FD0AHA53CC01AAP1 FD0AHC61CD02BBP1 FD0AHC39DC07BBP1 FD0AFA6BF04BBM1 FD0AHA3DA12.contig FD0AHA37BC07BBP1 FD0AHC41DF04.contig FD0ACA38CA07BBP1 FD0AFB9CG03.contig FD0ADA44DA02BBP1 FD0CCA14CE08BBP1 FD0AFA4BC08BBM1 FD0AHA46CD01AAP1 FD0AHA44CA10BBP1 FD0ABA12CD07.contig FD0ADA6BF10.contig FD0ADA8AD06.contig FD0AFB1BD06.contig FD0ADA9CH09BBP1 FD0ACA26CE04BBP1 FD0AFA3CD03BBM1 FD0AHC30AF02.contig FD0AHA53AF07BBP1 FD0AHA45DF10BBP1 FD0AHC22BC11BBM1 FD0AHA14BC06.contig FD0AEB17BF06.contig FD0AHC42AC10AAP1 FD0AHC73BA10AAP1);


print STDOUT "********************************\n";

foreach my $non (@nontet)
{
    print STDOUT ".$non\n";
}
