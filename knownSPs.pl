#!/usr/bin/perl -w


my %ids  = ();
my %sels = ();

open(FILE,"$ARGV[0]");
while(<FILE>)
{
    my $a =  $_;
    chomp $a;
    $ids{$a}++

}
close(FILE);

open(SEL, "$ARGV[1]");
while(<SEL>)
{
    my $b = $_;
    chomp $b;
    $sel{$b}++
}

my @keys = keys(%ids);

map {print STDOUT "$_\n" unless $sel{$_}}@keys;
