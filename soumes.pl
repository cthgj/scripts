#!/usr/bin/perl -w

while(<>)
{
    /^(.+?)\s+(.+?)\s+(.+)/;
    my $n=$1;
    my $p=$3;
    my @a=split(/\,/, $2);
#    print "$1:$2:$3 :  @a\n";die();
    my $k;
    map{$k+=$_}@a;
    my $o=$k-$p;
    print "$n : $o ($k-$p)\n";
    
	     

}
