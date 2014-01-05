#!/usr/bin/perl -w

# This script will take an exp sequence information file (genbank ncbi style) 
# and convert it to FASTA

while (<>)
{
    if(/^ID\s*(.*?)$/){
	print ">$1\n"
	} 
    elsif (/^\s/){
	s/\n/zpoz/g; s/\s//g; 
	s/\d//g; s/zpoz/\n/g;  
	print;
    }
    else {
	next;
    }
}
        
