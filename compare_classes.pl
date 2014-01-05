#!/usr/bin/perl

use strict;
use Term::ANSIColor; 

my $file1=$ARGV[0]||die("Need at least two filenames!\n");
my $file2=$ARGV[1]||die("Need at least two filenames!\n");

my %classes;
my ($f1_classes, $f2_classes)=(0,0);
#########################
# Read first class file #
#########################
open(my $fh1, "<", $file1)||die "Could not open first file ($file2) : $!\n";
while (<$fh1>) {
    chomp;
    my @pp=split(/\s/);
    $classes{$.}{file1}=\@pp;
    $f1_classes++;
}

##########################
# Read second class file #
##########################
open(my $fh2, "<", $file2)||die "Could not open second file ($file2) : $!\n";
while (<$fh2>) {
    chomp;
    my @pp=split(/\s/);
    $classes{$.}{file2}=\@pp;
    $f2_classes++;    
}

print "## File 1 Classes : $f1_classes\n";
print "## File 2 Classes : $f2_classes\n";

########################################################
# In order to cycle through all the classes, I need to #
# check which file has the greater number of classes.  #
########################################################
my $class_num;
if ($f1_classes > $f2_classes){
    $class_num=$f1_classes;
}
else {
    $class_num=$f2_classes;
}   


###############################################
# Count the number of proteins shared between #
# each pair of classes.			      #
###############################################
my %shared;
my (%best1, %best2);
for (my $i=1; $i<=$f1_classes; $i++) {
    my @class1=@{$classes{$i}{file1}};
    my $num1=scalar(@class1);
    for (my $k=1; $k<=$f2_classes; $k++) {
	my @class2=@{$classes{$k}{file2}};
	my $num2=scalar(@class2);
	my $same=0;
	## For each protein in class1
	foreach my $p1 (@class1) {
	    ## For each protein in class2
	    foreach my $p2 (@class2) {
		$same++ if $p1 eq $p2;
	    }
	}
	####################################################
        # If we already have a value for this class	   #
	# update it if the current one is better.	   #
        ####################################################
	if (defined($best1{$i})) {
	    $best1{$i}{CLASS}=$k if $same > $best1{$i}{SHARED};
	    $best1{$i}{SHARED}=$same if $same > $best1{$i}{SHARED};
	}
	else {
	    $best1{$i}{CLASS}=$k;
	    $best1{$i}{SHARED}=$same;
  	}
	###########################################################
        # And the second class (maybe we dont always have 	  #
	# reciprocal best hits.					  #
        ###########################################################
	if (defined($best1{$k})) {
	    $best2{$k}{CLASS}=$i if $same > $best1{$k}{SHARED};
	    $best2{$k}{SHARED}=$same if $same > $best1{$k}{SHARED};
	}
	else {
	    $best2{$k}{CLASS}=$i;
	    $best2{$k}{SHARED}=$same;
  	}

	##Sort class names to get a uniqe pair name
	my $pair=join("_",sort($i,$k));
	$shared{$i}{$k}=$same;
    }
}


for (my $c1=1; $c1<=$class_num;$c1++) {

#foreach my $c1 (sort { $a <=> $b } (keys(%best1))){
    ##################################
    # Is this a reciprocal best hit? #
    ##################################
    if ($best1{$c1}{CLASS} eq $best2{$c1}{CLASS}) {
	print "YES $c1 : $best1{$c1}{CLASS}($best1{$c1}{SHARED})\n";
    }
    else {
	print "$c1 : $best1{$c1}{CLASS}($best1{$c1}{SHARED}) : $best2{$c1}{CLASS}($best2{$c1}{SHARED})\n";
    }
} 






############################################
# Find the pairs of classes that share the #
# most proteins.			   #
############################################
for (my $i=1; $i<=$f1_classes; $i++) {
         $best1{$i}=0;

     
 }
# foreach my $class_pair (keys(%shared)) {
#     print "$class_pair : $shared{$class_pair}\n";
# }    



# for (my $i=1; $i<=$class_num; $i++) {
#     ####################################################
#     # The number of elements in this class from file 1 #
#     ####################################################
#     my $num1=0;
#     my @class1;
#     defined($classes{$i}{file1}) && do {
# 	$num1=scalar(@{$classes{$i}{file1}});
# 	@class1=@{$classes{$i}{file1}};
#     };
#     ####################################################
#     # The number of elements in this class from file 2 #
#     ####################################################
#     my $num2=0;
#     my @class2;
#     defined($classes{$i}{file2}) && do {
# 	$num2=scalar(@{$classes{$i}{file2}});
# 	@class2=@{$classes{$i}{file2}};
#     };
#     #################################################
#     # Count how many proteins the two classes share #
#     #################################################
#     my $same=0;
#     foreach my $p1 (@class1) {
# 	foreach my $p2 (@class2) {
# 	    $same++ if $p1 eq $p2;
# 	}
#     }
#     if ($same == $num1) {
# 	print color("green")."Class $i : $num1\t$num2\t$same\n".color("reset");
#     }
#     else {
# 	print "Class $i : $num1\t$num2\t$same\n";	
#     }

    
# }    
