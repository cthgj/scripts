#!/usr/bin/perl -w

##########################################################
# This script will take an OCG class file and convert it #
# to Clust&See format. 					 #
##########################################################
usage() unless $ARGV[1];



my $name=$ARGV[0];
print "#ModClust analysis export\n#Network:$name\n#Scope:Network\n\n";
open(A,"$ARGV[1]") || die "Could not open file $ARGV[1] : $!\n";

while (<A>) {
    chomp;
    next if /^\#/;
    /^\[.+?(\d+)/&& do{print ">ClusterID:$1||\n"};
    if (/^CA/ ||/^CM/) {
	/^..\s+(.+)/;
	my $kk=$1;
	my @annot=split(/\s+/,$kk);
#	print "@annot|\n";
    }
    elsif (/^PN\s+(.+)/) {
	my @p=split(/, /,$1);
	local $"="\n";
	print "@p\n\n";
    }
    
    
}






sub usage{
    $0=~/.+?([^\/]+)$/;
    my $p_name=$1;
    print <<EndOfHelp;
USAGE: 
\t$p_name <network name> <class file>

EndOfHelp
exit;
}





