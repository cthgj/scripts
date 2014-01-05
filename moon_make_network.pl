#!/usr/bin/perl -w

use strict;
use Getopt::Std;

my (%opts);
getopts('MvTgncnhSl:s:N:w:f:C:m:p:r:R:',\%opts);
my $DATADIR=$opts{D}||"/home/cchapple/research/testing/new/data";
my $net_file=$opts{N}||"$DATADIR/../HSHQ.gr";
my $cand=$ARGV[0];

$net_file=~/.+\/(.+?)\..{2,3}/;
my $net_name=$1;
my $subnet_file=$cand . "_" . $net_name .  ".sif";
die("$subnet_file exists") if -e $subnet_file;
`grep -w $cand $net_file | gawk 'BEGIN{OFS="\t"}{print \$1,"pp",\$2}' > $subnet_file`;

my $attr_file=$cand  . "_" . $net_name . ".class.attr";
my $attr_file=$cand  . "_" . $net_name . ".class.attr";
my $attr_file1=$cand  . "_" . $net_name . ".true_class.attr";
my $attr_file2=$cand  . "_" . $net_name . ".names.attr";
open(C,"> $attr_file")||die("Could not open $attr_file for writing:$!\n");
open(C1,"> $attr_file1")||die("Could not open $attr_file1 for writing:$!\n");
open(C2,"> $attr_file2")||die("Could not open $attr_file2 for writing:$!\n");
print C "Class (class=String)\n$cand" . "_$cand = 0\n";
print C1 "True_Class (class=String)\n";
print C2 "Label\n$cand" . "_$cand = $cand\n";
print C1 "$cand" . "_$cand = ($class1" . "::$class2)\n" if $mode eq 'i';
foreach my $name (keys(%names)){
    my $NetName=$cand . "_" . $name;
    print C2 "$NetName = $name\n";
    if (defined($names{$name}{$class1}) && defined($names{$name}{$class2})){
	print C "$NetName = AB\n";
	print C1 "$NetName = ($class1" . "::$class2)\n";
    }
    elsif (defined($names{$name}{$class1})){
	print C "$NetName = A\n";
	print C1 "$NetName = $class1\n";
    }
    elsif(defined($names{$name}{$class2})){
	print C "$NetName = B\n";
	print C1 "$NetName = $class2\n";
    }
    else{die("Huh!!??\n$NetName,$class1,$class2,$class1,$class2\n");}
}
close(C);








sub make_networks{
 #    my $cand=shift;
#     my $mode=shift;
#     my $i=shift;
# #    my $cand_attr=$cand . "." . $mode . ".cand.attr";
#  #   open(A,">$cand_attr");
#     if(($mode eq 'i')||($mode eq 'B') || ($mode eq 'b') || ($mode eq 'a')){
# 	my $class1=shift;
# 	my $class2;
# 	$mode eq 'a' ? ($class2='none') : ($class2=shift);
# 	my %names;
# 	open(A,"moon_get_class_prots.pl $class1 $class_file |");
# 	while(<A>){
# 	    chomp;
# 	    next if $cand eq $_;
# 	    $names{$_}{$class1}++;
# 	}
# 	close(A);
# 	unless ($mode eq 'a'){
# 	    open(B,"moon_get_class_prots.pl $class2 $class_file |");
# 	    while(<B>){
# 		chomp;
# 		next if $cand eq $_;
# 		$names{$_}{$class2}++;
# 	    }
# 	    close(B);
# 	}
# 	my $attr_file=$cand . "." . $mode . ".$i" . ".class.attr";
#         my $attr_file1=$cand . "." . $mode . ".$i" . ".true_class.attr";
#         my $attr_file2=$cand . "." . $mode . ".$i" . ".names.attr";
# 	open(C,"> $attr_file")||die("Could not open $attr_file for writing:$!\n");
# 	open(C1,"> $attr_file1")||die("Could not open $attr_file1 for writing:$!\n");
# 	open(C2,"> $attr_file2")||die("Could not open $attr_file2 for writing:$!\n");
# 	print C "Class (class=String)\n$cand" . "_$cand = 0\n";
# 	print C1 "True_Class (class=String)\n";
# 	print C2 "Label\n$cand" . "_$cand = $cand\n";
# 	print C1 "$cand" . "_$cand = ($class1" . "::$class2)\n" if $mode eq 'i';
# 	foreach my $name (keys(%names)){
# 	    my $NetName=$cand . "_" . $name;
# 	    print C2 "$NetName = $name\n";
# 	    if (defined($names{$name}{$class1}) && defined($names{$name}{$class2})){
# 		print C "$NetName = AB\n";
# 		print C1 "$NetName = ($class1" . "::$class2)\n";
# 	    }
# 	    elsif (defined($names{$name}{$class1})){
# 		print C "$NetName = A\n";
# 		print C1 "$NetName = $class1\n";
# 	    }
# 	    elsif(defined($names{$name}{$class2})){
# 		print C "$NetName = B\n";
# 		print C1 "$NetName = $class2\n";
# 	    }
# 	    else{die("Huh!!??\n$NetName,$class1,$class2,$class1,$class2\n");}
# 	}
# 	close(C);
# 	my $subnet_file=$cand .  $mode . ".$i.sif";
# 	open(S,">$subnet_file")||die "Could not open $subnet_file:$!\n";
# 	open(N,"$network_file")||die("Could not open $network_file:$!\n");
# 	while(<N>){
# 	    next if /^\d+$/;
# 	    next if /^\s+$/;
# 	    chomp;
# 	    my @a=split(/\t/);
# 	    if(((($a[0] eq $cand) || ($a[1] eq $cand)) && (defined($names{$a[0]}) ||defined($names{$a[1]})))||
# 	       (defined($names{$a[0]}) && defined($names{$a[1]}))){
# 		print S "$cand" . "_$a[0] pp $cand" . "_$a[1]\n";
# 	    }
# 	}
# 	close(S);
# 	close(N);
	
# 	my $cyto_script=$cand . ".$i." . $mode . ".cyto";
# 	open(C3,">$cyto_script")||die "Could not open $cyto_script:$!\n";
# 	my $cur_dir=`pwd`;
# 	chomp($cur_dir);
# 	print C3 "session open file=\"/home/cchapple/research/testing/new/results/candidates.cys\"\n";
# 	print C3 "network import file=\"$cur_dir/$subnet_file\"\n";
# 	print C3 "node import attributes file=\"$cur_dir/$attr_file\"\n";
# 	print C3 "node import attributes file=\"$cur_dir/$attr_file1\"\n";
# 	print C3 "node import attributes file=\"$cur_dir/$attr_file2\"\n";

#     } 
}


##
