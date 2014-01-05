#!/usr/bin/perl 
use warnings;
use strict;

my @bins=(1e-2,1e-10,1e-20,1e-30,1e-40,1e-50,1e-60,1e-70,1e-80,1e-90,1e-100,1e-110,1e-120,1e-130,1e-140,1e-150,1e-160,1e-170,1e-180,1e-190,1e-200,1e-210,1e-220,1e-230,1e-240,1e-250,1e-260,1e-270,1e-280,1e-290,1e-300,1e-310,1e-320,1e-330,1e-340,1e-350,1e-360,1,1e-10,1e-160,1e-210,1e-260,1e-310,0);

while(<>){
    chomp; 
    for(my $i=1; $i<=$#bins; $i++){ 
	if($bins[$i]<$_){
#	    print "$bins[$i-1] ($_,$bins[$i],$bins[$i-1],$i,$.)\n";
	    print "$bins[$i-1]\n";
	    $i=1020;
	}
	# elsif($bins[$i]>$_){

	# }
	# else{

	# }
    }
}
