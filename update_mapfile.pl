#!/usr/bin/perl -w

use strict;
use Getopt::Std;

print STDERR "DOESN'T WORK, CHECK SOURCE\n";
die();

unless($ARGV[0]){
    print STDERR "This script will take 2 map files (old_name\tnew_name\taccession) as input and replace any names in the second that are different in the first. Intended to be used in conjunction with uniprot_fix_ids.pl. \nUSAGE: uniprot_fix_ids.pl names | $0 -m old_mapfile\n ";
    exit(0);
}
my %opts;
getopts('m:',\%opts);
my $map=$opts{m}||die("Need a map file to update (-m)\n");

my(%newname,%oldname);


while(my $line=<>){
   chomp($line);
   my @a=split(/\t/,$line);
   my $newline=$a[2] . "\t"
    ## $accs{old_acc}=new_acc
   map{$newname{$_}=$line; print STDERR "\$newname{$_}=$line\n"}@a;
}


open(MAP,"$map")|| die("Could not open $map :$!\n");
while(<MAP>){
    chomp;
    my @a=split(/\t/,$_);
    ## $accs{old_acc}=new_acc
    my $c='bob';
    map{$c=$_ if defined($newname{$_})}@a;
    $c eq 'bob' ?
	print "$_\n" : 
	print "nn $_ : $newname{$c}\n" ;
	
}
close(MAP);

