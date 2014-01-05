#!/usr/bin/env perl

my $target=$ARGV[0] || die("Need a target as $ARGV[0]\n");
my $netstats;#=$ARGV[1] || die("Need a string as $ARGV[1]\n");
my $netdir="/var/www_8080/PoGo/data";
my $urldir="data";
my $netstats='<ul id="netstats">';



foreach (qw(human fly yeast worm mouse)) {
    my $prots=`gawk '{print \$1\"\\n\"\$2}' $netdir/$_.gr | sort | uniq | wc -l`;
    chomp($prots);
    my $ints=`cat $netdir/$_.gr | wc -l`;
    chomp($ints);
    $netstats.="<li><a href=\"$urldir/$_.gr\">$_</a>, $ints interactions between $prots proteins</li>";
}
$netstats.="</ul>";

open(A,"$target")||die("Could not open file $target : $!\n");
while (<A>) {
    if(/id=.netstats/){
	print "$netstats";
    }
    else{print}
}
