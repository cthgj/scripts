#!/usr/bin/perl -w

my %want;
open(A,"$ARGV[0]");
while(<A>){
    chomp;
     $want{$_}++;
}


my $class;
while(<>){
    if(/\[CLASS:\s*(\d+)/){$class=$1; }
    elsif(/^CA\s+(.+)/){
	print "$class : $1\n" if defined($want{$class});
    }
}
