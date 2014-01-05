#!/usr/bin/perl 

# This script will take

%T = (); 
open(FILE, "$ARGV[0]");

while (<FILE>)
{
    chomp;
    @F = split /\s+/og, $_;  ### @F[0] is now TOG# and all other @F[n] species name
    
    for ($n=1; $n < scalar(@F); $n++) ### Cycle through @F, creating anonymous array references...
    {
	defined($T{$F[$n]}) || (@{ $T{$F[$n]}} = ()); push @{ $T{$F[$n]} }, $F[0]; 

    };
    print STDERR ".";  
    
}
 	foreach $k (sort keys %T) 
	{ 
	    print STDOUT "$k: @{ $T{$k} }\n"; 
	}; 
    
