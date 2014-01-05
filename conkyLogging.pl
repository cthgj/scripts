#!/usr/bin/perl
my $lim=32;
my @a;
while(<>){/.*?\]\s*(.+)$/;
 push @a,$1;
}
my $k;
for($n=scalar(@a)-1;$n>=0; $n--){
    $_=$a[$n];
    @c=split(/[\s]+/);
    @b=split(//);
    my $k=0;
    my $kk=0;
    print "    ";
    for($i=0;$i<=$#b; $i++){
	$_=$b[$i];
	if (/^\s+$/){
	    $k+=length($c[$kk+1]);
	    $kk++ ;
	}

	if($k>$lim){#print "\nKK($k>$lim,$kk) :: $c[$kk] :: $c[$kk+1]\n";
	    s/([=,\-\s])/$1\n\t  /;
	    $k=0;
	}
	
	print STDOUT "$_";
    }
#    if($k<$lim) {print "\n";}
    print STDOUT "\n";
}
#print "\n";    if($k<$lim) {}
