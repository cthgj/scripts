#!/usr/bin/perl


use strict;

my $name = 'ha';
my %sel  = ();
my @list;
my $sel;

while(<>)
{
    /^(.*?)\s(.*?)$/og;
    
    if ($name eq 'ha'){
	$sel = $1;
	$name = $2; 
	push @{$sel{$sel}}, $name; 
    }
    elsif ($sel eq $1){

	push @{$sel{$sel}}, $2; 
    }
   
    else
    {
	print STDERR "\$sel : $sel\n\@list : @list\n";
	
	$sel = $1;
	$name = $2;
	push @{$sel{$sel}}, $name; 
    }
}

my @keys = sort(keys(%sel));


map {print STDOUT "$_ : @{$sel{$_}}\n"} @keys;
