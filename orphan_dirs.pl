#!/usr/bin/perl -w
## This script will find all subdirectories in the current dir
## and delete those dirs that have less than 3 files


use strict;
my @dirs;

system("find . -type d > listofdirs");
open(F,"listofdirs");
while(<F>){
    if(/^\.+$/){next}
    elsif(/\/\./){next}
    chomp; 
    unshift @dirs, $_;
}

print "###### GOT DIRS ########\n";
 
map{
    opendir(A,"$_"); 
    my @a=readdir(A); 
    close(A);
    
    if($#a-1<3){ 
	
#	system( "rm -rf $_");
#	system( "rm \"$_\"/*; rmdir \"$_\"");
	print "rm $_/*; rmdir \"$_\"\n";
	print "Deleting $_ (",$#a-1," files)\n";
    }
}@dirs;
close(F);
#unlink("listofdirs");


