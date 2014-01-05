#!/usr/bin/perl -w

# removes repeated sequences from a tbl file
# repeated : diff name same seq


my %hash = ();
my %ids = ();
my @id;

while(<>)
{
    /^(.*?)\s+([^\s]*)$/g;
    $hash{$2} = $1;
    push @{$ids{$2}}, $1;
}
my @keys = keys(%hash);
map
{
    for(my $i=0; $i<(scalar(@{$ids{$_}})); $i++){
	print STDOUT "${$ids{$_}}[$i]";
	print STDOUT "------" if ${$ids{$_}}[$i+1];
	}
    print STDOUT "\t$_\n";
} @keys;
