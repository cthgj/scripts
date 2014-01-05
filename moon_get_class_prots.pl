#!/usr/bin/perl -w

my (%prots,%classes);

unless($ARGV[0]){
    print STDERR "This script will take a list of class names and a class file as input and will return the proteins in each of the input classes. If a protein belongs to multiple classes, it will be printed ONCE\n\n";
    exit();
}

if(-e $ARGV[0]){
    open(A,"$ARGV[0]")||die("Cannot open $ARGV[0] : $!\n");
    while(<A>){
	chomp;
	$classes{$_}++;
    }
}
else{
    $classes{$ARGV[0]}++;
}


close(A);

open(B,"$ARGV[1]")||die("Cannot open $ARGV[1] : $!\n");
my $class;
while(<B>){
    next if /^\#/;
    chomp;
    
    if(/^\[CLASS:\s*(\d+)/){$class=$1}
    # if(/^\P#\s*(\d+)/){
	
    # }
    if(/^PN\s*(.+)/){
	my @a=split(/,\s/,$1);
	map{$prots{$_}++}@a if defined($classes{$class});
    }
}
close(B);

map{print "$_\n"}keys(%prots);
