#!/usr/bin/perl -w

use Getopt::Std;
my (%opts,%prots,%classes,%want);
getopts("nc:h", \%opts);

my $no_zero=$opts{n}||undef;
my %k; 
my $n=$ARGV[0]; 
my $file=$ARGV[1];
open(A,"$file")||die("Cannot open $file : $!\n");
while(<A>){
    chomp;
    if(/^(.+?)\s+=\s+\((.+)\)/){
    my ($foo,$bar)=($1,$2); 
    my @a=split(/::/,$bar); 
    map{$k{$foo}{$_}++; }@a; 
    }
}
close(A);
print "MultiClasses\n";
my @lp=keys(%{$k{$n}});
foreach my $p (keys(%k)){
    my $line= "$p = ("; 
    my %c;
    my @oo=keys(%{$k{$p}});
    foreach my $class (@oo){
	if(defined($k{$n}{$class})){
	    $c{$class}++ ;
	}
	else{$c{0}++;}
    }
    my @cc=keys(%c);
    my @classes=();
    for (my $ii=0; $ii<=$#cc; $ii++){
	if($no_zero){next if $cc[$ii] eq 0;}
	push @classes,$cc[$ii];
    }
    if($#classes==0){$line.="$classes[0])";}
#    elsif($#classes==1){$line.="$classes[0])"; print "cc : @classes\n"; die();}
    else{
	for (my $ii=0; $ii<=$#classes; $ii++){
	    $line.="$classes[$ii]";
	    $line.="::" unless $ii==$#classes;
	}
	 $line.=")";
    }
    print "$line\n" unless $line=~/\(\)/;

    # if($no_zero){
    # 	    $line.=  "$cc[0]" unless $cc[0] eq 0;
    # 	}
    # 	else{
    # 	    $line.=   "$cc[0]";
    # 	}
    # for (my $ii=1; $ii<$#cc; $ii++){print "ii $ii : $#cc\n";
    # 	if($no_zero){
    # 	    $line.=  "::$cc[$ii]" unless $cc[$ii] eq 0;
    # 	}
    # 	else{
    # 	    $line.=   "::$cc[$ii]";
    # 	}
    # }
    # if($#cc>0){
    # 	if($no_zero){
    # 	    if($cc[$#cc] eq 0){
    # 		$line.=  ")" ;
    # 	    }
    # 	    else{
    # 		$line.=  "::$cc[$#cc])";
    # 	    }
    # 	}
    # 	else{
    # 	    $line.=  "::$cc[$#cc])";
    # 	}
    # }
    # else{   $line.=  ")";}
}

